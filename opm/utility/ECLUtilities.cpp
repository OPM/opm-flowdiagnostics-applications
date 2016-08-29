/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

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

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <opm/utility/ECLUtilities.hpp>


namespace Opm
{

// ======================================================================

std::vector<double>
ECL::getPVolVector(const ecl_grid_type* G,
                   const ecl_file_type* init)
{
    auto make_szt = [](const int i)
    {
        return static_cast<std::vector<double>::size_type>(i);
    };

    const auto nx = make_szt(ecl_grid_get_nx(G));
    const auto ny = make_szt(ecl_grid_get_ny(G));
    const auto nz = make_szt(ecl_grid_get_nz(G));

    const auto nglob = nx * ny * nz;

    auto pvol = std::vector<double>(nglob, 1.0);

    if (ecl_file_has_kw(init, "PORV")) {
        auto porv =
            ecl_file_iget_named_kw(init, "PORV", 0);

        assert ((make_szt(ecl_kw_get_size(porv)) == nglob)
                && "Pore-volume must be provided for all global cells");

        ecl_kw_get_data_as_double(porv, pvol.data());
    }

    return pvol;
}

ECL::GridPtr
ECL::loadCase(const boost::filesystem::path& grid)
{
    auto G = GridPtr{
        ecl_grid_load_case(grid.generic_string().c_str())
    };

    if (! G) {
        std::ostringstream os;

        os << "Failed to load ECL Grid from '"
           << grid.generic_string() << '\'';

        throw std::invalid_argument(os.str());
    }

    return G;
}

ECL::FilePtr
ECL::loadFile(const boost::filesystem::path& file)
{
    // Read-only, keep open between requests
    const auto open_flags = 0;

    auto F = FilePtr{
        ecl_file_open(file.generic_string().c_str(), open_flags)
    };

    if (! F) {
        std::ostringstream os;

        os << "Failed to load ECL File object from '"
           << file.generic_string() << '\'';

        throw std::invalid_argument(os.str());
    }

    return F;
}

std::vector<ecl_nnc_type>
ECL::loadNNC(const ecl_grid_type* G,
             const ecl_file_type* init)
{
    auto make_szt = [](const int n)
    {
        return static_cast<std::vector<ecl_nnc_type>::size_type>(n);
    };

    auto nncData = std::vector<ecl_nnc_type>{};

    const auto numNNC = make_szt(ecl_nnc_export_get_size(G));

    if (numNNC > 0) {
        nncData.resize(numNNC);

        ecl_nnc_export(G, init, nncData.data());
    }

    std::sort(nncData.begin(), nncData.end(),
              [](const ecl_nnc_type& nd1, const ecl_nnc_type& nd2) {
                  return nd1.input_index < nd2.input_index;
              });

    return nncData;
}

// ======================================================================

ECL::ScatterMap::ScatterMap(const ecl_grid_type*       G,
                            const std::vector<double>& pvol)
{
    auto make_szt = [](const int i)
    {
        return static_cast<std::size_t>(i);
    };

    this->nact_ = make_szt(ecl_grid_get_nactive(G));

    {
        const auto nx = make_szt(ecl_grid_get_nx(G));
        const auto ny = make_szt(ecl_grid_get_ny(G));
        const auto nz = make_szt(ecl_grid_get_nz(G));

        this->nglob_ = nx * ny * nz;
    }

    subset_.clear();
    subset_.reserve(this->nact_);

    if (pvol.empty()) {
        for (decltype(ecl_grid_get_nactive(G))
                 act = 0, nact = ecl_grid_get_nactive(G);
             act < nact; ++act)
        {
            const auto glob =
                make_szt(ecl_grid_get_global_index1A(G, act));

            this->subset_.push_back(ID{ make_szt(act), glob });
        }
    }
    else {
        for (decltype(ecl_grid_get_nactive(G))
                 act = 0, nact = ecl_grid_get_nactive(G);
             act < nact; ++act)
        {
            const auto glob =
                make_szt(ecl_grid_get_global_index1A(G, act));

            if (pvol[glob] > 0.0) {
                this->subset_.push_back(ID{ make_szt(act), glob });
            }
        }
    }
}

std::vector<std::size_t>
ECL::ScatterMap::activeGlobalCells() const
{
    auto active = std::vector<std::size_t>{};
    active.reserve(this->subset_.size());

    for (const auto& id : this->subset_) {
        active.push_back(id.glob);
    }

    return active;
}

template <typename T>
std::vector<T>
ECL::ScatterMap::scatterToGlobal(const std::vector<T>& x) const
{
    // Assume that input vector 'x' is either defined on explicit
    // notion of active cells (ACTNUM != 0) or on all global cells or
    // some other contiguous index set (e.g., the NNCs).

    if (x.size() != static_cast<decltype(x.size())>(this->nact_)) {
        // Input not defined on explictly active cells.  Let caller
        // deal with it.  This typically corresponds to the set of
        // global cells or the list of NNCs.
        return x;
    }

    auto y = std::vector<T>(this->nglob_);

    for (const auto& i : this->subset_) {
        y[i.glob] = x[i.act];
    }

    return y;
}

// ======================================================================

ECL::GlobalCellData::GlobalCellData(const ecl_grid_type*       grid,
                                    const std::vector<double>& pvol)
    : map_(grid, pvol)
{
}

void
ECL::GlobalCellData::assignDataSource(const boost::filesystem::path& src)
{
    this->src_ = ECL::loadFile(src);
}

void
ECL::GlobalCellData::assignDataSource(ECL::FilePtr src)
{
    this->src_ = std::move(src);
}

const ecl_file_type*
ECL::GlobalCellData::getDataSource() const
{
    return this->src_.get();
}

std::vector<std::size_t>
ECL::GlobalCellData::activeGlobalCells() const
{
    return this->map_.activeGlobalCells();
}

bool
ECL::GlobalCellData::haveVector(const std::string& vector) const
{
    return ecl_file_has_kw(this->src_.get(), vector.c_str());
}

const ECL::GlobalCellData::DataVector&
ECL::GlobalCellData::getVector(const std::string& vector,
                               const int          occurrence)
{
    auto& p = this->coll_[vector][occurrence];

    if (! p) {
        auto x = this->map_.scatterToGlobal(load(vector, occurrence));

        p = DataHandle{ new DataVector{ std::move(x) } };
    }

    return *p;
}

ECL::GlobalCellData::DataVector
ECL::GlobalCellData::load(const std::string& vector,
                          const int          occurrence)
{
    auto x = ecl_file_iget_named_kw(this->src_.get(),
                                    vector.c_str(),
                                    occurrence);

    auto y = DataVector(ecl_kw_get_size(x));

    ecl_kw_get_data_as_double(x, y.data());

    return y;
}

// ======================================================================

CellCollection::CellCollection(const ecl_grid_type*       G,
                               std::vector<std::size_t>&& active_glob_cells)
    : nx_       (static_cast<std::size_t>(ecl_grid_get_nx(G)))
    , ny_       (static_cast<std::size_t>(ecl_grid_get_ny(G)))
    , nz_       (static_cast<std::size_t>(ecl_grid_get_nz(G)))
    , glob_cell_(std::move(active_glob_cells))
    , active_id_(nx_ * ny_ * nz_, -1)
{
    auto act = 0;

    for (const auto& id : glob_cell_) {
        this->active_id_[static_cast<std::size_t>(id)] = act++;
    }
}

const std::vector<std::size_t>&
CellCollection::activeGlobalCells() const
{
    return this->glob_cell_;
}

std::size_t
CellCollection::numActive() const
{
    return this->glob_cell_.size();
}

std::size_t
CellCollection::numGlobalCells() const
{
    return this->active_id_.size();
}

int
CellCollection::activeCell(const std::size_t globalCell) const
{
    if (globalCell >= numGlobalCells()) { return -1; }

    return this->active_id_[globalCell];
}

int
CellCollection::cartesianNeighbour(const std::size_t globalCell,
                                   const Direction   d) const
{
    if (globalCell >= numGlobalCells()) { return -1; }

    auto ijk = ind2sub(globalCell);

    if      (d == Direction::I) { ijk[0] += 1; }
    else if (d == Direction::J) { ijk[1] += 1; }
    else if (d == Direction::K) { ijk[2] += 1; }
    else {
        return -1;
    }

    const auto globNeigh = globIdx(ijk);

    if (globNeigh >= numGlobalCells()) { return -1; }

    return this->active_id_[globNeigh];
}

std::size_t
CellCollection::globIdx(const IJKTuple& ijk) const
{
    if (ijk[0] >= this->nx_) { return -1; }
    if (ijk[1] >= this->ny_) { return -1; }
    if (ijk[2] >= this->nz_) { return -1; }

    return ijk[0] + nx_*(ijk[1] + ny_*ijk[2]);
}

CellCollection::IJKTuple
CellCollection::ind2sub(const std::size_t globalCell) const
{
    assert (globalCell < numGlobalCells());

    auto ijk = IJKTuple{};
    auto g   = globalCell;

    ijk[0] = g % this->nx_;  g /= this->nx_;
    ijk[1] = g % this->ny_;
    ijk[2] = g / this->ny_;  assert (ijk[2] < this->nz_);

    assert (globIdx(ijk) == globalCell);

    return ijk;
}

} // namespace Opm
