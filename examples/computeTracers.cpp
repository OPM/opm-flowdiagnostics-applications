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

#include "exampleSetup.hpp"


int main(int argc, char** argv)
try {
    auto fdTool = example::setup(argc, argv);

    // Solve for tracers.
    Opm::FlowDiagnostics::CellSetID id("Example start set ID");
    Opm::FlowDiagnostics::CellSet start;
    start.identify(id);
    start.insert(123);
    auto sol = fdTool.computeInjectionDiagnostics({start});
    const auto& conc_sol = sol.fd.concentration(id);

    // Write it to standard out.
    std::cout.precision(16);
    const int num_c = conc_sol.cellValueCount();
    for (int ii = 0; ii < num_c; ++ii) {
        auto p = conc_sol.cellValue(ii);
        std::cout << p.first << "\t" << p.second << '\n';
    }
}
catch (const std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << '\n';
}
