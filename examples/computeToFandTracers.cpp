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

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/flowdiagnostics/ConnectivityGraph.hpp>
#include <opm/flowdiagnostics/ConnectionValues.hpp>
#include <opm/flowdiagnostics/Toolbox.hpp>

#include <opm/utility/ECLGraph.hpp>

#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/filesystem.hpp>

namespace {
    std::vector<std::string>
    commandArguments(int argc, char* argv[])
    {
        // Initialise vector<> from range implied by pair of (begin,end)
        // iterators.  We don't need POSIX guarantee that
        //
        //     argv[argc] == nullptr
        return { argv, argv + argc };
    }

    bool isFile(const boost::filesystem::path& p)
    {
        namespace fs = boost::filesystem;

        auto is_regular_file = [](const fs::path& pth)
        {
            return fs::exists(pth) && fs::is_regular_file(pth);
        };

        return is_regular_file(p)
            || (fs::is_symlink(p) &&
                is_regular_file(fs::read_symlink(p)));
    }

    boost::filesystem::path
    deriveFileName(boost::filesystem::path         file,
                   const std::vector<std::string>& extensions)
    {
        for (const auto& ext : extensions) {
            file.replace_extension(ext);

            if (isFile(file)) {
                return file;
            }
        }

        const auto prefix = file.parent_path() / file.stem();

        std::ostringstream os;

        os << "Unable to derive valid filename from model prefix "
           << prefix.generic_string();

        throw std::invalid_argument(os.str());
    }

    boost::filesystem::path
    gridFile(const std::vector<std::string>& cmd_args)
    {
        if (cmd_args.size() < 2) {
            throw std::invalid_argument("Must supply at least prefix "
                                        "of ECLIPSE output dataset");
        }

        return cmd_args[1];
    }

    boost::filesystem::path
    initFile(const std::vector<std::string>& cmd_args)
    {
        if (cmd_args.size() < 3) {
            return deriveFileName(gridFile(cmd_args),
                                  { "INIT", "FINIT" });
        }

        return cmd_args[2];
    }

    boost::filesystem::path
    restartFile(const std::vector<std::string>& cmd_args)
    {
        if (cmd_args.size() < 4) {
            return deriveFileName(gridFile(cmd_args),
                                  { "UNRST", "FUNRST" });
        }

        return cmd_args[3];
    }

    int
    stepNumber(const std::vector<std::string>& cmd_args)
    {
        if (cmd_args.size() < 5) {
            return 0;
        }

        return std::stoi(cmd_args[4]);
    }

    Opm::FlowDiagnostics::ConnectionValues
    extractFluxField(/* mutable */ Opm::ECLGraph& G,
                     const int                    step)
    {
        using ConnVals = Opm::FlowDiagnostics::ConnectionValues;

        using NConn = ConnVals::NumConnections;
        using NPhas = ConnVals::NumPhases;

        const auto nconn = NConn{G.numConnections()};
        const auto nphas = NPhas{3};

        auto flux = ConnVals(nconn, nphas);

        auto phas = ConnVals::PhaseID{0};

        for (const auto& p : { Opm::BlackoilPhases::Aqua   ,
                               Opm::BlackoilPhases::Liquid ,
                               Opm::BlackoilPhases::Vapour })
        {
            const auto pflux = G.flux(p, step);

            if (! pflux.empty()) {
                assert (pflux.size() == nconn.total);

                auto conn = ConnVals::ConnID{0};

                for (const auto& v : pflux) {
                    flux(conn, phas) = v;

                    conn.id += 1;
                }
            }

            phas.id += 1;
        }

        return flux;
    }

    Opm::FlowDiagnostics::Toolbox
    initialiseFlowDiagnostics(/* mutable */ Opm::ECLGraph& G, const int step)
    {
        const auto connGraph = Opm::FlowDiagnostics::
            ConnectivityGraph{ static_cast<int>(G.numCells()),
                               G.neighbours() };

        using FDT  = Opm::FlowDiagnostics::Toolbox;
        using PVol = FDT::PoreVolume;
        using Flux = FDT::ConnectionFlux;

        auto tool = FDT{ connGraph };

        tool.assign(PVol{ G.poreVolume() })
            .assign(Flux{ extractFluxField(G, step) });

        return tool;
    }
} // Anonymous namespace

// Syntax:
//   computeToFandTracers ecl_case_prefix

int main(int argc, char* argv[])
try {
    const auto cmd_args = commandArguments(argc, argv);

    auto graph = Opm::ECLGraph::load(gridFile(cmd_args),
                                     initFile(cmd_args));

    graph.assignFluxDataSource(restartFile(cmd_args));

    auto fdTool = initialiseFlowDiagnostics(graph, stepNumber(cmd_args));

    // Solve for time of flight.
    using FDT = Opm::FlowDiagnostics::Toolbox;
    std::vector<Opm::FlowDiagnostics::CellSet> start;
    auto sol = fdTool.computeInjectionDiagnostics(FDT::StartCells{start});
    const auto& tof = sol.fd.timeOfFlight();

    // Write it to standard out.
    for (double t : tof) {
        std::cout << t << '\n';
    }
}
catch (const std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << '\n';
}
