/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil ASA.

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

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <opm/flowdiagnostics-applications/ECLGraph.hpp>

#include <memory>
#include <utility>

#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_grid.h>
#include <ert/ecl/ecl_kw_magic.h>

class Opm::Utility::ECLGraph::Impl
{
public:
    Impl(const Path& grid, const Path& init);

private:

    template <class T, void (*release)(T* obj)>
    struct DestroyECL
    {
        void operator()(T* obj)
        {
            release(obj);
        }
    };

    template <class T, void (*release)(T* obj)>
    using ECLPtr  = std::unique_ptr<T, DestroyECL<T, release>>;

    using GridPtr = ECLPtr<ecl_grid_type, ecl_grid_free>;
    using FilePtr = ECLPtr<ecl_file_type, ecl_file_close>;
};

// ======================================================================

Opm::Utility::ECLGraph::ECLGraph(ECLGraph&& rhs)
    : pImpl_(std::move(rhs.pImpl_))
{}

Opm::Utility::ECLGraph::~ECLGraph()
{}

Opm::Utility::ECLGraph&
Opm::Utility::ECLGraph::operator=(ECLGraph&& rhs)
{
    pImpl_ = std::move(rhs.pImpl_);

    return *this;
}

Opm::Utility::ECLGraph
Opm::Utility::ECLGraph::load(const Path& grid, const Path& init)
{
    auto pImpl = ImplPtr{new Impl(grid, init)};

    return { std::move(pImpl) };
}

std::size_t
Opm::Utility::ECLGraph::numCells() const
{
    return pImpl_->numCells();
}

const std::vector<int>&
Opm::Utility::ECLGraph::neighbours() const
{
    return pImpl_->neighbours();
}
