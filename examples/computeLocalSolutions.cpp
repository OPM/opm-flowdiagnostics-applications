/*
  Copyright 2017 SINTEF ICT, Applied Mathematics.
  Copyright 2017 Statoil ASA.

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

#include "exampleSetup.hpp"
#include <opm/flowdiagnostics/CellSet.hpp>

#include <fstream>
#include <ios>
#include <type_traits>

namespace {
    using StartSets = std::vector< ::Opm::FlowDiagnostics::CellSet>;

    template <class WellData>
    std::vector<int>
    activeCompletions(const ::Opm::ECLGraph& G,
                      const WellData&        well)
    {
        auto completion_cells = std::vector<int>{};

        completion_cells.reserve(well.completions.size());

        for (const auto& completion : well.completions) {
            const auto cell_index =
                G.activeCell(completion.ijk, completion.gridName);

            if (cell_index >= 0) {
                completion_cells.push_back(cell_index);
            }
        }

        return completion_cells;
    }

    template <class Select>
    StartSets
    getStartPointsFromWells(const example::Setup& setup,
                            Select&&              pick)
    {
        auto start = StartSets{};

        for (const auto& well : setup.well_fluxes) {
            if (pick(well)) {
                start.emplace_back(Opm::FlowDiagnostics::CellSetID(well.name),
                                   activeCompletions(setup.graph, well));
            }
        }

        return start;
    }

    StartSets injectors(const example::Setup& setup)
    {
        using WData =
            std::decay<decltype(setup.well_fluxes[0])>::type;

        return getStartPointsFromWells(setup, [](const WData& well)
                                       { return well.is_injector_well; });
    }

    StartSets producers(const example::Setup& setup)
    {
        using WData =
            std::decay<decltype(setup.well_fluxes[0])>::type;

        return getStartPointsFromWells(setup, [](const WData& well)
                                       { return ! well.is_injector_well; });
    }

    void printSolution(const ::Opm::FlowDiagnostics::CellSetValues& x,
                       const ::Opm::FlowDiagnostics::CellSetID&     id,
                       const std::string&                           varname)
    {
        const auto filename =
            varname + '-' + id.to_string() + ".out";

        std::ofstream os(filename);

        if (os) {
            os.precision(16);
            os.setf(std::ios_base::scientific);

            for (const auto& item : x) {
                os << item.first << ' ' << item.second << '\n';
            }
        }
    }

    std::vector<int>
    extractCellIDs(const ::Opm::FlowDiagnostics::CellSetValues& x)
    {
        auto i = std::vector<int>{};
        i.reserve(x.size());

        for (const auto& xi : x) { i.push_back(xi.first); }

        std::sort(std::begin(i), std::end(i));

        return i;
    }

    bool
    sameReachability(const ::Opm::FlowDiagnostics::CellSetValues& tof,
                     const ::Opm::FlowDiagnostics::CellSetValues& conc)
    {
        if (tof.size() != conc.size()) {
            return false;
        }

        const auto tof_id  = extractCellIDs(tof);
        const auto conc_id = extractCellIDs(conc);

        return tof_id == conc_id;
    }

    void runAnalysis(const StartSets&                        start,
                     const ::Opm::FlowDiagnostics::Solution& sol)
    {
        auto ok = std::vector<bool>{};
        ok.reserve(start.size());

        for (const auto& pt : start) {
            const auto& id = pt.id();

            const auto& tof  = sol.timeOfFlight (id);
            const auto& conc = sol.concentration(id);

            printSolution(tof , id, "tof" );
            printSolution(conc, id, "conc");

            if (! sameReachability(tof, conc)) {
                std::cout << id.to_string() << ": FAIL\n";
            }
        }
    }
}

// Syntax (typical):
//   computeLocalSolutions case=<ecl_case_prefix> step=<report_number>
int main(int argc, char* argv[])
try {
    example::Setup setup(argc, argv);
    auto& fdTool = setup.toolbox;

    {
        const auto inj = injectors(setup);
        const auto fwd = fdTool.computeInjectionDiagnostics(inj);

        runAnalysis(inj, fwd.fd);
    }

    {
        const auto prod = producers(setup);
        const auto rev  = fdTool.computeProductionDiagnostics(prod);

        runAnalysis(prod, rev.fd);
    }
}
catch (const std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << '\n';
}
