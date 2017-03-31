/*
  Copyright 2017 SINTEF ICT, Applied Mathematics.
  Copyright 2017 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <opm/utility/ECLPropTable.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <exception>
#include <iterator>
#include <stdexcept>
#include <utility>

Opm::ECLPropTable1D::SingleTable::
SingleTable(ElmIt               xBegin,
            ElmIt               xEnd,
            std::vector<ElmIt>& colIt)
{
    // There must be at least one dependent variable/result variable.
    assert (colIt.size() >= 1);

    const auto nRows = std::distance(xBegin, xEnd);

    this->x_.reserve(nRows);
    this->y_.reserve(nRows * colIt.size());

    auto keyValid = [](const double xi)
    {
        // Indep. variable values <= -1.0e20 or >= 1.0e20 signal "unused"
        // table nodes (rows).  These nodes are in the table to fill out the
        // allocated size if one particular sub-table does not use all
        // nodes.  The magic value 1.0e20 is documented in the Fileformats
        // Reference Manual.
        return std::abs(xi) < 1.0e20;
    };

    while (xBegin != xEnd) {
        // Extract relevant portion of the table.  Preallocated rows that
        // are not actually part of the result set (i.e., those that are set
        // to a sentinel value) are discarded.
        if (keyValid(*xBegin)) {
            this->x_.push_back(*xBegin);

            for (auto ci : colIt) {
                // Store 'y_' with column index cycling most rapidly.
                this->y_.push_back(*ci);
            }
        }

        // -------------------------------------------------------------
        // Advance iterators.

        // 1) Independent variable.
        ++xBegin;

        // 2) Dependent/result/columns.
        for (auto& ci : colIt) {
            ++ci;
        }
    }

    // Dispose of any excess capacity.
    if (this->x_.size() < static_cast<decltype(this->x_.size())>(nRows)) {
        this->x_.shrink_to_fit();
        this->y_.shrink_to_fit();
    }

    if (this->x_.size() < 2) {
        // Table has no interval that supports interpolation.  Either just a
        // single node or no nodes at all.  We can't do anything useful
        // here, so don't pretend that this is okay.

        throw std::invalid_argument {
            "No Interpolation Intervals of Non-Zero Size"
        };
    }
}

double
Opm::ECLPropTable1D::SingleTable::
y(const ECLPropTableRawData::SizeType nCols,
  const ECLPropTableRawData::SizeType row,
  const ResultColumn&                 c) const
{
    assert (row * nCols < this->y_.size());
    assert (c.i < nCols);

    // Recall: 'y_' stored with column index cycling the most rapidly (row
    // major ordering).
    return this->y_[row*nCols + c.i];
}

std::vector<double>
Opm::ECLPropTable1D::SingleTable::
interpolate(const ECLPropTableRawData::SizeType nCols,
            const ResultColumn&                 c,
            const std::vector<double>&          x) const
{
    auto y = std::vector<double>{};  y.reserve(x.size());

    auto yval = [nCols, c, this]
        (const ECLPropTableRawData::SizeType i)
    {
        return this->y(nCols, i, c);
    };

    const auto first = ECLPropTableRawData::SizeType{ 0 };
    const auto last  = ECLPropTableRawData::SizeType{ this->x_.size() - 1 };

    for (const auto& xi : x) {
        y.push_back(0.0);
        auto& yi = y.back();

        if (! (xi > this->x_.front())) {
            // Constant extrapolation to the left of range.
            yi = yval(first);
        }
        else if (! (xi < this->x_.back())) {
            // Constant extrapolation to the right of range.
            yi = yval(last);
        }
        else {
            // Somewhere in [min(x_), max(x_)].  Primary key (indep. var) is
            // sorted range.  Recall: lower_bound() returns insertion point,
            // which translates to the *upper* (right-hand) end-point of the
            // interval in this context.
            auto b = std::begin(this->x_);
            auto p = std::lower_bound(b, std::end(this->x_), xi);

            assert ((p != b) && "Logic Error Left End-Point");
            assert ((p != std::end(this->x_)) &&
                    "Logic Error Right End-Point");

            const auto i  = p - b;
            const auto xl = this->x_[i - 1];
            const auto yl = yval(i - 1);
            const auto yr = yval(i + 0);

            const auto t = (xi - xl) / (this->x_[i + 0] - xl);

            yi = (1.0 - t)*yl + t*yr;
        }
    }

    return y;
}

// =====================================================================

Opm::ECLPropTable1D::ECLPropTable1D(const ECLPropTableRawData& raw)
    : nResCols_(raw.numCols - 1)
{
    if (raw.numCols < 2) {
        throw std::invalid_argument {
            "Malformed Property Table"
        };
    }

    this->table_.reserve(raw.numTables);

    // Table format: numRows*numTables values of first column (indep. var)
    // followed by numCols-1 dependent variable (function value result)
    // columns of numRows*numTables values each, one column at a time.
    const auto colStride = raw.numRows * raw.numTables;

    // Position
    auto xBegin = std::begin(raw.data);
    auto colIt  = std::vector<decltype(xBegin)>{ xBegin + colStride };
    for (auto col = 0*raw.numCols + 1; col < raw.numCols - 1; ++col) {
        colIt.push_back(colIt.back() + colStride);
    }

    for (auto t = 0*raw.numTables;
              t <   raw.numTables;
         ++t, xBegin += raw.numRows)
    {
        auto xEnd = xBegin + raw.numRows;

        // Note: The SingleTable ctor advances each 'colIt' across numRows
        // entries.  That is a bit of a layering violation, but helps in the
        // implementation of this loop.
        this->table_.push_back(SingleTable(xBegin, xEnd, colIt));
    }
}

std::vector<double>
Opm::ECLPropTable1D::interpolate(const InTable&             t,
                                 const ResultColumn&        c,
                                 const std::vector<double>& x) const
{
    if (t.i >= this->table_.size()) {
        throw std::invalid_argument {
            "Invalid Table ID"
        };
    }

    if (c.i >= this->nResCols_) {
        throw std::invalid_argument {
            "Invalid Result Column ID"
        };
    }

    return this->table_[t.i].interpolate(this->nResCols_, c, x);
}
