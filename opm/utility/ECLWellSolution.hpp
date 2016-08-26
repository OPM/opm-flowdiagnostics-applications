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

#ifndef OPM_ECLWELLSOLUTION_HEADER_INCLUDED
#define OPM_ECLWELLSOLUTION_HEADER_INCLUDED

#include <ert/ecl/ecl_file.h>
#include <ert/util/ert_unique_ptr.hpp>
#include <boost/filesystem.hpp>
#include <array>
#include <string>
#include <utility>
#include <vector>

namespace Opm
{

    class ECLWellSolution
    {
    public:
        using Path = boost::filesystem::path;

        /// Construct with path to restart file.
        explicit ECLWellSolution(const Path& restart_filename);

        /// Contains the well data extracted from the restart file.
        struct WellData
        {
            std::string name;
            struct Completion
            {
                std::array<int, 3> ijk;
                double reservoir_rate;
            };
            std::vector<Completion> completions;
        };

        /// Return well solution for given occurrence (time step).
        std::vector<WellData> solution(const int occurrence);

    private:
        using FilePtr = ERT::ert_unique_ptr<ecl_file_type, ecl_file_close>;
        FilePtr restart_;
        std::vector<double> loadDoubleField(const std::string& fieldname);
        std::vector<int> loadIntField(const std::string& fieldname);
    };


} // namespace Opm

#endif // OPM_ECLWELLSOLUTION_HEADER_INCLUDED
