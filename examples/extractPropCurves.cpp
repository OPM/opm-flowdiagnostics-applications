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

#include <examples/exampleSetup.hpp>

#include <opm/utility/ECLCaseUtilities.hpp>
#include <opm/utility/ECLPhaseIndex.hpp>
#include <opm/utility/ECLPropertyUnitConversion.hpp>
#include <opm/utility/ECLPvtCommon.hpp>
#include <opm/utility/ECLPvtCurveCollection.hpp>
#include <opm/utility/ECLResultData.hpp>
#include <opm/utility/ECLSaturationFunc.hpp>

#include <array>
#include <cstddef>
#include <exception>
#include <initializer_list>
#include <iomanip>
#include <ios>
#include <iostream>
#include <regex>
#include <vector>

#include <boost/filesystem.hpp>

namespace {
    template <class OStream>
    void printGraph(OStream&                                        os,
                    const std::string&                              indep,
                    const std::string&                              dep,
                    const std::vector<Opm::FlowDiagnostics::Graph>& graphs)
    {
        const auto oprec  = os.precision(16);
        const auto oflags = os.setf(std::ios_base::scientific);

        for (const auto& graph : graphs) {
            const auto& x = graph.first;
            const auto& y = graph.second;

            os << indep << ' ' << dep << '\n';

            for (auto n = x.size(), i = 0*n; i < n; ++i) {
                os << x[i] << ' ' << y[i] << '\n';
            }

            os << "\n\n";
        }

        os.setf(oflags);
        os.precision(oprec);
    }

    template <class OStream>
    void printGraph(OStream&                                  os,
                    const std::string&                        pressure,
                    const std::string&                        mixRat,
                    const std::string&                        function,
                    const std::vector<Opm::ECLPVT::PVTGraph>& graphs)
    {
        const auto oprec  = os.precision(16);
        const auto oflags = os.setf(std::ios_base::scientific);

        for (const auto& graph : graphs) {
            const auto& p = graph.press;
            const auto& R = graph.mixRat;
            const auto& f = graph.value;

            os << pressure << ' ' << mixRat << ' ' << function << '\n';

            for (auto n = p.size(), i = 0*n; i < n; ++i) {
                os << p[i] << ' ' << R[i] << ' ' << f[i] << '\n';
            }

            os << "\n\n";
        }

        os.setf(oflags);
        os.precision(oprec);
    }

    // -----------------------------------------------------------------
    // Relative permeability

    using CurveSet = Opm::ECLSaturationFunc::RawCurve::CurveSet;

    std::string curveName(const CurveSet     curveSet,
                          const std::string& baseName)
    {
        if (curveSet == CurveSet::Imbibition) {
            return "I" + baseName;
        }
        else {
            return baseName;
        }
    }

    void krg(const Opm::ECLSaturationFunc&                 sfunc,
             const Opm::ECLInitFileData&                   init,
             const std::string&                            gridID,
             const int                                     activeCell,
             const Opm::ECLSaturationFunc::SatFuncScaling& scaling,
             const CurveSet                                curveSet)
    {
        using RC = Opm::ECLSaturationFunc::RawCurve;

        auto func = std::vector<RC>{};
        func.reserve(1);

        // Request krg (gas rel-perm in oil-gas system)
        func.push_back(RC{
            RC::Function::RelPerm,
            RC::SubSystem::OilGas,
            Opm::ECLPhaseIndex::Vapour,
            curveSet,
        });

        const auto graph = sfunc
            .getSatFuncCurve(func, init, gridID, activeCell, scaling);

        printGraph(std::cout, "Sg", curveName(curveSet, "Krg"), graph);
    }

    void krog(const Opm::ECLSaturationFunc&                 sfunc,
              const Opm::ECLInitFileData&                   init,
              const std::string&                            gridID,
              const int                                     activeCell,
              const Opm::ECLSaturationFunc::SatFuncScaling& scaling,
              const CurveSet                                curveSet)
    {
        using RC = Opm::ECLSaturationFunc::RawCurve;

        auto func = std::vector<RC>{};
        func.reserve(1);

        // Request krog (oil rel-perm in oil-gas system)
        func.push_back(RC{
            RC::Function::RelPerm,
            RC::SubSystem::OilGas,
            Opm::ECLPhaseIndex::Liquid,
            curveSet,
        });

        const auto graph = sfunc
            .getSatFuncCurve(func, init, gridID, activeCell, scaling);

        printGraph(std::cout, "Sg", curveName(curveSet, "Krog"), graph);
    }

    void krow(const Opm::ECLSaturationFunc&                 sfunc,
              const Opm::ECLInitFileData&                   init,
              const std::string&                            gridID,
              const int                                     activeCell,
              const Opm::ECLSaturationFunc::SatFuncScaling& scaling,
              const CurveSet                                curveSet)
    {
        using RC = Opm::ECLSaturationFunc::RawCurve;

        auto func = std::vector<RC>{};
        func.reserve(1);

        // Request krow (oil rel-perm in oil-water system)
        func.push_back(RC{
            RC::Function::RelPerm,
            RC::SubSystem::OilWater,
            Opm::ECLPhaseIndex::Liquid,
            curveSet,
        });

        const auto graph = sfunc
            .getSatFuncCurve(func, init, gridID, activeCell, scaling);

        printGraph(std::cout, "Sw", curveName(curveSet, "Krow"), graph);
    }

    void krw(const Opm::ECLSaturationFunc&                 sfunc,
             const Opm::ECLInitFileData&                   init,
             const std::string&                            gridID,
             const int                                     activeCell,
             const Opm::ECLSaturationFunc::SatFuncScaling& scaling,
             const CurveSet                                curveSet)
    {
        using RC = Opm::ECLSaturationFunc::RawCurve;

        auto func = std::vector<RC>{};
        func.reserve(1);

        // Request krw (water rel-perm in oil-water system)
        func.push_back(RC{
            RC::Function::RelPerm,
            RC::SubSystem::OilWater,
            Opm::ECLPhaseIndex::Aqua,
            curveSet,
        });

        const auto graph = sfunc
            .getSatFuncCurve(func, init, gridID, activeCell, scaling);

        printGraph(std::cout, "Sw", curveName(curveSet, "Krw"), graph);
    }

    // -----------------------------------------------------------------
    // Capillary pressure

    void pcgo(const Opm::ECLSaturationFunc&                 sfunc,
              const Opm::ECLInitFileData&                   init,
              const std::string&                            gridID,
              const int                                     activeCell,
              const Opm::ECLSaturationFunc::SatFuncScaling& scaling,
              const CurveSet                                curveSet)
    {
        using RC = Opm::ECLSaturationFunc::RawCurve;

        auto func = std::vector<RC>{};
        func.reserve(1);

        // Request pcgo (gas/oil capillary pressure (Pg-Po) in G/O system)
        func.push_back(RC{
            RC::Function::CapPress,
            RC::SubSystem::OilGas,
            Opm::ECLPhaseIndex::Vapour,
            curveSet,
        });

        const auto graph = sfunc
            .getSatFuncCurve(func, init, gridID, activeCell, scaling);

        printGraph(std::cout, "Sg", curveName(curveSet, "Pcgo"), graph);
    }

    void pcgw(const Opm::ECLSaturationFunc&                 sfunc,
              const Opm::ECLInitFileData&                   init,
              const std::string&                            gridID,
              const int                                     activeCell,
              const Opm::ECLSaturationFunc::SatFuncScaling& scaling,
              const CurveSet                                curveSet)
    {
        using RC = Opm::ECLSaturationFunc::RawCurve;

        auto func = std::vector<RC>{};
        func.reserve(1);

        // Request pcgw (gas/water capillary pressure (Pg-Pw) in G/W system)
        func.push_back(RC{
            RC::Function::CapPress,
            RC::SubSystem::GasWater,
            Opm::ECLPhaseIndex::Vapour,
            curveSet,
        });

        const auto graph = sfunc
            .getSatFuncCurve(func, init, gridID, activeCell, scaling);

        printGraph(std::cout, "Sg", curveName(curveSet, "Pcgw"), graph);
    }

    void pcow(const Opm::ECLSaturationFunc&                 sfunc,
              const Opm::ECLInitFileData&                   init,
              const std::string&                            gridID,
              const int                                     activeCell,
              const Opm::ECLSaturationFunc::SatFuncScaling& scaling,
              const CurveSet                                curveSet)
    {
        using RC = Opm::ECLSaturationFunc::RawCurve;

        auto func = std::vector<RC>{};
        func.reserve(1);

        // Request pcow (oil/water capillary pressure (Po-Pw) in O/W system)
        func.push_back(RC{
            RC::Function::CapPress,
            RC::SubSystem::OilWater,
            Opm::ECLPhaseIndex::Aqua,
            curveSet,
        });

        const auto graph = sfunc
            .getSatFuncCurve(func, init, gridID, activeCell, scaling);

        printGraph(std::cout, "Sw", curveName(curveSet, "Pcow"), graph);
    }

    // -----------------------------------------------------------------
    // PVT Curves

    void Bg(const Opm::ECLPVT::ECLPvtCurveCollection& pvtCurves,
            const int                                 pvtNum)
    {
        using RC = Opm::ECLPVT::RawCurve;

        const auto graph = pvtCurves
            .getPvtCurve(RC::FVF, Opm::ECLPhaseIndex::Vapour, pvtNum);

        printGraph(std::cout, "Pg", "Rv", "Bg", graph);
    }

    void mu_g(const Opm::ECLPVT::ECLPvtCurveCollection& pvtCurves,
              const int                                 pvtNum)
    {
        using RC = Opm::ECLPVT::RawCurve;

        const auto graph = pvtCurves
            .getPvtCurve(RC::Viscosity, Opm::ECLPhaseIndex::Vapour, pvtNum);

        printGraph(std::cout, "Pg", "Rv", "mu_g", graph);
    }

    void Bo(const Opm::ECLPVT::ECLPvtCurveCollection& pvtCurves,
            const int                                 pvtNum)
    {
        using RC = Opm::ECLPVT::RawCurve;

        const auto graph = pvtCurves
            .getPvtCurve(RC::FVF, Opm::ECLPhaseIndex::Liquid, pvtNum);

        printGraph(std::cout, "Po", "Rs", "Bo", graph);
    }

    void mu_o(const Opm::ECLPVT::ECLPvtCurveCollection& pvtCurves,
              const int                                 pvtNum)
    {
        using RC = Opm::ECLPVT::RawCurve;

        const auto graph = pvtCurves
            .getPvtCurve(RC::Viscosity, Opm::ECLPhaseIndex::Liquid, pvtNum);

        printGraph(std::cout, "Po", "Rs", "mu_o", graph);
    }

    // -----------------------------------------------------------------
    // Saturated states (RvSat(Pg) and RsSat(Po))

    void rvSat(const Opm::ECLPVT::ECLPvtCurveCollection& pvtCurves,
               const int                                 pvtNum)
    {
        using RC = Opm::ECLPVT::RawCurve;
        using PI = Opm::ECLPhaseIndex;

        const auto graph = pvtCurves
            .getPvtCurve(RC::SaturatedState, PI::Vapour, pvtNum);

        printGraph(std::cout, "Pg", "Rv", "rvSat", graph);
    }

    void rsSat(const Opm::ECLPVT::ECLPvtCurveCollection& pvtCurves,
               const int                                 pvtNum)
    {
        using RC = Opm::ECLPVT::RawCurve;
        using PI = Opm::ECLPhaseIndex;

        const auto graph = pvtCurves
            .getPvtCurve(RC::SaturatedState, PI::Liquid, pvtNum);

        printGraph(std::cout, "Po", "Rs", "rsSat", graph);
    }

    // -----------------------------------------------------------------

    std::unique_ptr<const Opm::ECLUnits::UnitSystem>
    makeUnits(const std::string&          unit,
              const Opm::ECLInitFileData& init)
    {
        if ((unit == "si") || (unit == "SI") || (unit == "internal")) {
            return {};          // No conversion needed.
        }

        if ((unit == "metric") || (unit == "Metric") || (unit == "METRIC")) {
            return Opm::ECLUnits::metricUnitConventions();
        }

        if ((unit == "field") || (unit == "Field") || (unit == "FIELD")) {
            return Opm::ECLUnits::fieldUnitConventions();
        }

        if ((unit == "lab") || (unit == "Lab") || (unit == "LAB")) {
            return Opm::ECLUnits::labUnitConventions();
        }

        if ((unit == "pvt-m") || (unit == "PVT-M") || (unit == "PVTM")) {
            return Opm::ECLUnits::pvtmUnitConventions();
        }

        std::cerr << "Unit convention '" << unit << "' not recognized\n"
                  << "Using 'native' (input/serialised) conventions.\n";

        return Opm::ECLUnits::serialisedUnitConventions(init);
    }

    int getActiveCell(const Opm::ECLInitFileData&             init,
                      const Opm::ECLCaseUtilities::ResultSet& rset,
                      const Opm::ParameterGroup&              prm)
    {
        if (prm.has("cell")) { return prm.get<int>("cell"); }

        const auto graph = Opm::ECLGraph::load(rset.gridFile(), init);
        const auto grid = prm.getDefault("gridName", std::string{});

        auto ijk = std::array<int,3>{ { 0, 0, 0 } };

        if      (prm.has("I")) { ijk[0] = prm.get<int>("I") - 1; } // 1-based
        else if (prm.has("i")) { ijk[0] = prm.get<int>("i") - 0; } // 0-based

        if      (prm.has("J")) { ijk[1] = prm.get<int>("J") - 1; } // 1-based
        else if (prm.has("j")) { ijk[1] = prm.get<int>("j") - 0; } // 0-based

        if      (prm.has("K")) { ijk[2] = prm.get<int>("K") - 1; } // 1-based
        else if (prm.has("k")) { ijk[2] = prm.get<int>("k") - 0; } // 0-based

        const auto acell = graph.localActiveCell(ijk, grid);

        if (acell < 0) {
            std::cerr << "Cell ("
                      << (ijk[0] + 1) << ", "
                      << (ijk[1] + 1) << ", "
                      << (ijk[2] + 1) << ") "
                      << "in "
                      << (grid.empty()
                          ? std::string{"Main Grid"}
                          : "Local Grid: " + grid)
                      << " is inactive.\n"
                      << "Using Cell 0 in Main Grid\n";

            return 0;
        }

        return acell;
    }

    Opm::ECLSaturationFunc::SatFuncScaling
    saturationFuncScaling(const Opm::ParameterGroup& prm)
    {
        using T = Opm::ECLSaturationFunc::SatFuncScaling::Type;

        const auto opt =
            std::regex_constants::icase |
            std::regex_constants::ECMAScript;

        const auto horiz  = std::regex { "horizontal", opt };
        const auto vert   = std::regex { "vertical"  , opt };
        const auto useEPS =
            prm.getDefault("useEPS", std::string{"off"});

        auto scaling = Opm::ECLSaturationFunc::SatFuncScaling{};
        scaling.enable = static_cast<unsigned char>(0);

        if (std::regex_search(useEPS, horiz)) {
            scaling.enable |= T::Horizontal;
        }

        if (std::regex_search(useEPS, vert)) {
            scaling.enable |= T::Vertical;
        }

        return scaling;
    }
} // namespace Anonymous

int main(int argc, char* argv[])
try {
    const auto prm     = example::initParam(argc, argv);
    const auto scaling = saturationFuncScaling(prm);

    const auto rset  = example::identifyResultSet(prm);
    const auto init  = Opm::ECLInitFileData(rset.initFile());

    const auto gridID = prm.getDefault("gridName", std::string{});
    const auto cellID = getActiveCell(init, rset, prm);

    const auto pvtNum = init.haveKeywordData("PVTNUM")
        ? init .keywordData<int>("PVTNUM", gridID)[cellID] : 1;
    const auto haveHyst = init.haveHysteresis();

    auto sfunc = Opm::ECLSaturationFunc(init);
    auto pvtCC = Opm::ECLPVT::ECLPvtCurveCollection(init);

    if (prm.has("unit")) {
        auto units = makeUnits(prm.get<std::string>("unit"), init);
        
        sfunc.setOutputUnits(units->clone());
        pvtCC.setOutputUnits(std::move(units));
    }

    // -----------------------------------------------------------------
    // Relative permeability

    if (prm.getDefault("krg", false)) {
        krg(sfunc, init, gridID, cellID, scaling, CurveSet::Drainage);
    }

    if (prm.getDefault("ikrg", false)) {
        if (haveHyst) {
            krg(sfunc, init, gridID, cellID, scaling, CurveSet::Imbibition);
        }
        else {
            std::cerr << "IKrg not available in this result set\n";
        }
    }

    if (prm.getDefault("krog", false)) {
        krog(sfunc, init, gridID, cellID, scaling, CurveSet::Drainage);
    }

    if (prm.getDefault("ikrog", false)) {
        if (haveHyst) {
            krog(sfunc, init, gridID, cellID, scaling, CurveSet::Imbibition);
        }
        else {
            std::cerr << "IKrog not available in this result set\n";
        }
    }

    if (prm.getDefault("krow", false)) {
        krow(sfunc, init, gridID, cellID, scaling, CurveSet::Drainage);
    }

    if (prm.getDefault("ikrow", false)) {
        if (haveHyst) {
            krow(sfunc, init, gridID, cellID, scaling, CurveSet::Imbibition);
        }
        else {
            std::cerr << "IKrow not available in this result set\n";
        }
    }

    if (prm.getDefault("krw", false)) {
        krw(sfunc, init, gridID, cellID, scaling, CurveSet::Drainage);
    }

    if (prm.getDefault("ikrw", false)) {
        if (haveHyst) {
            krw(sfunc, init, gridID, cellID, scaling, CurveSet::Imbibition);
        }
        else {
            std::cerr << "IKrw not available in this result set\n";
        }
    }

    // -----------------------------------------------------------------
    // Capillary pressure
    if (prm.getDefault("pcog", false) || // Alias pcog -> pcgo
        prm.getDefault("pcgo", false))
    {
        pcgo(sfunc, init, gridID, cellID, scaling, CurveSet::Drainage);
    }

    if (prm.getDefault("ipcog", false) || // Alias pcog -> pcgo
        prm.getDefault("ipcgo", false))
    {
        if (haveHyst) {
            pcgo(sfunc, init, gridID, cellID, scaling, CurveSet::Imbibition);
        }
        else {
            std::cerr << "IPcog not available in this result set\n";
        }
    }

    if (prm.getDefault("pcgw", false)) {
        pcgw(sfunc, init, gridID, cellID, scaling, CurveSet::Drainage);
    }

    if (prm.getDefault("ipcgw", false)) {
        if (haveHyst) {
            pcgw(sfunc, init, gridID, cellID, scaling, CurveSet::Imbibition);
        }
        else {
            std::cerr << "IPcgw not available in this result set\n";
        }
    }

    if (prm.getDefault("pcow", false)) {
        pcow(sfunc, init, gridID, cellID, scaling, CurveSet::Drainage);
    }

    if (prm.getDefault("ipcow", false)) {
        if (haveHyst) {
            pcow(sfunc, init, gridID, cellID, scaling, CurveSet::Imbibition);
        }
        else {
            std::cerr << "IPcow not available in this result set\n";
        }
    }

    // -----------------------------------------------------------------
    // PVT Curves

    if (prm.getDefault("Bg"  , false)) { Bg  (pvtCC, pvtNum); }
    if (prm.getDefault("mu_g", false)) { mu_g(pvtCC, pvtNum); }
    if (prm.getDefault("Bo"  , false)) { Bo  (pvtCC, pvtNum); }
    if (prm.getDefault("mu_o", false)) { mu_o(pvtCC, pvtNum); }

    if (prm.getDefault("rvSat", false)) { rvSat(pvtCC, pvtNum); }
    if (prm.getDefault("rsSat", false)) { rsSat(pvtCC, pvtNum); }
}
catch (const std::exception& e) {
    std::cerr << "Caught Exception: " << e.what() << '\n';

    return EXIT_FAILURE;
}
