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

#ifndef OPM_ECLGRAPH_HEADER_INCLUDED
#define OPM_ECLGRAPH_HEADER_INCLUDED

#include <cstddef>
#include <memory>
#include <vector>

#include <boost/filesystem.hpp>

namespace Opm { namespace Utility {

    class ECLGraph
    {
    public:
        using Path = boost::filesystem::path;

        ECLGraph() = delete;

        ~ECLGraph();

        ECLGraph(ECLGraph&& rhs);
        ECLGraph(const ECLGraph& rhs) = delete;

        ECLGraph& operator=(ECLGraph&& rhs);
        ECLGraph& operator=(const ECLGraph& rhs) = delete;

        static ECLGraph
        load(const Path& grid, const Path& init);

        std::size_t             numCells()  const;
        const std::vector<int>& neigbours() const;

    private:
        class Impl;

        using ImplPtr = std::unique_ptr<Impl>;

        ECLGraph(ImplPtr pImpl);
        
        ImplPtr pImpl_;
    };

}} // namespace Opm::Utility

#endif // OPM_ECLGRAPH_HEADER_INCLUDED
