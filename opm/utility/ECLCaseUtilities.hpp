/*
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

#ifndef OPM_ECLCASEUTILITIES_HEADER_INCLUDED
#define OPM_ECLCASEUTILITIES_HEADER_INCLUDED

#include <vector>

#include <boost/filesystem/path.hpp>

namespace Opm { namespace ECLCaseUtilities {

    class ResultSet
    {
    public:
        using Path = boost::filesystem::path;

        explicit ResultSet(const Path& casename);

        Path gridFile() const;
        Path initFile() const;
        Path restartFile(const int reportStepID) const;

        bool isUnifiedRestart() const;

        std::vector<int> reportStepIDs() const;

    private:
        /// Case prefix.  Typical form is "path/to/case", but might be
        /// "path/to/case.01" too.
        boost::filesystem::path prefix_;

        /// Whether or not this result set has a unified restart file
        /// (.UNRST, .FUNRST).
        bool isUnified_;
    };

}} // Opm::ECLCaseUtilities

#endif // OPM_ECLCASEUTILITIES_HEADER_INCLUDED
