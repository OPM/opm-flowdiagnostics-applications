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

#include <opm/utility/ECLPvtCurveCollection.hpp>

#include <opm/utility/ECLResultData.hpp>

#include <cassert>
#include <initializer_list>
#include <vector>

#include <opm/parser/eclipse/Units/Units.hpp>

#include <ert/ecl/ecl_kw_magic.h>

namespace {
    std::vector<int>
    pvtnumVector(const ::Opm::ECLGraph&        G,
                 const ::Opm::ECLInitFileData& init)
    {
        auto pvtnum = G.rawLinearisedCellData<int>(init, "PVTNUM");

        if (pvtnum.empty()) {
            // PVTNUM missing in one or more of the grids managed by 'G'.
            // Put all cells in PVTNUM region 1.
            pvtnum.assign(G.numCells(), 1);
        }

        return pvtnum;
    }

    std::unique_ptr<const Opm::ECLUnits::UnitSystem>
    createUnitSystem(const ::Opm::ECLInitFileData& init)
    {
        const auto& ih = init.keywordData<int>(INTEHEAD_KW);

        return ::Opm::ECLUnits::createUnitSystem(ih[ INTEHEAD_UNIT_INDEX ]);
    }

    std::vector<Opm::FlowDiagnostics::Graph> emptyFDGraph()
    {
        return { Opm::FlowDiagnostics::Graph{} };
    }

    template <class PvtPtr>
    std::vector<Opm::FlowDiagnostics::Graph>
    rawPvtCurve(PvtPtr&                          pvt, // ref to shared_ptr<>
                const Opm::ECLPVT::RawCurve      curve,
                const int                        regID,
                const Opm::ECLUnits::UnitSystem& usys)
    {
        if (pvt != nullptr) {
            return pvt->getPvtCurve(curve, regID, usys);
        }

        // Result set does not provide requisite tabulated properties.
        // Return empty.
        return emptyFDGraph();
    }

    std::vector<double>
    fromSI(std::vector<double>&& x, const double x_scale)
    {
        for (auto& xi : x) {
            xi = Opm::unit::convert::to(xi, x_scale);
        }

        return x;
    }

    template <class PVTInterp, class Pressure, class MixRatio>
    std::vector<double>
    formationVolumeFactor(const PVTInterp& pvt,
                          const int        regID,
                          const Pressure&  press,
                          const MixRatio&  R,
                          const double     fvfScale)
    {
        return fromSI(pvt.formationVolumeFactor(regID, R, press), fvfScale);
    }

    template <class PVTInterp, class Pressure, class MixRatio>
    std::vector<double>
    viscosity(const PVTInterp& pvt,
              const int        regID,
              const Pressure&  press,
              const MixRatio&  R,
              const double     muScale)
    {
        return fromSI(pvt.viscosity(regID, R, press), muScale);
    }

    std::vector<double>
    gasProperty(const std::shared_ptr<Opm::ECLPVT::Gas>& pvt,
                const Opm::ECLPVT::RawCurve              property,
                const int                                regID,
                const std::vector<double>&               Pg,
                const std::vector<double>&               Rv,
                const Opm::ECLUnits::UnitSystem&         usys)
    {
        if (pvt == nullptr) {
            // Nu such property interpolant.  Return empty.
            return {};
        }

        assert ((property == Opm::ECLPVT::RawCurve::FVF) ||
                (property == Opm::ECLPVT::RawCurve::Viscosity));

        const auto pg = Opm::ECLPVT::Gas::GasPressure { Pg };
        auto       rv = Opm::ECLPVT::Gas::VaporizedOil{ Rv };
        if (rv.data.empty()) {
            rv.data.assign(pg.data.size(), 0.0);
        }

        if (property == Opm::ECLPVT::RawCurve::FVF) {
            const auto fvfScale = usys.reservoirVolume()
                / usys.surfaceVolumeGas();

            return formationVolumeFactor(*pvt, regID, pg, rv, fvfScale);
        }

        return viscosity(*pvt, regID, pg, rv, usys.viscosity());
    }

    std::vector<double>
    oilProperty(const std::shared_ptr<Opm::ECLPVT::Oil>& pvt,
                const Opm::ECLPVT::RawCurve              property,
                const int                                regID,
                const std::vector<double>&               Po,
                const std::vector<double>&               Rs,
                const Opm::ECLUnits::UnitSystem&         usys)
    {
        if (pvt == nullptr) {
            // Nu such property interpolant.  Return empty.
            return {};
        }

        assert ((property == Opm::ECLPVT::RawCurve::FVF) ||
                (property == Opm::ECLPVT::RawCurve::Viscosity));

        const auto po = Opm::ECLPVT::Oil::OilPressure { Po };
        auto       rs = Opm::ECLPVT::Oil::DissolvedGas{ Rs };
        if (rs.data.empty()) {
            rs.data.assign(po.data.size(), 0.0);
        }

        if (property == Opm::ECLPVT::RawCurve::FVF) {
            const auto fvfScale = usys.reservoirVolume()
                / usys.surfaceVolumeLiquid();

            return formationVolumeFactor(*pvt, regID, po, rs, fvfScale);
        }

        return viscosity(*pvt, regID, po, rs, usys.viscosity());
    }
}

Opm::ECLPVT::ECLPvtCurveCollection::
ECLPvtCurveCollection(const ECLGraph&        G,
                      const ECLInitFileData& init)
    : pvtnum_(pvtnumVector(G, init))
    , gas_   (CreateGasPVTInterpolant::fromECLOutput(init)) // u_p<> -> s_p<>
    , oil_   (CreateOilPVTInterpolant::fromECLOutput(init)) // u_p<> -> s_p<>
    , usys_  (createUnitSystem(init))                       // u_p<> -> s_p<>
{}

std::vector<Opm::FlowDiagnostics::Graph>
Opm::ECLPVT::ECLPvtCurveCollection::
getPvtCurve(const RawCurve      curve,
            const ECLPhaseIndex phase,
            const int           activeCell) const
{
    if (phase == ECLPhaseIndex::Aqua) {
        // Not supported at this time.
        // Return empty.
        return emptyFDGraph();
    }

    if (static_cast<decltype(this->pvtnum_.size())>(activeCell)
        >= this->pvtnum_.size())
    {
        // Active cell index out of bounds.  Return empty.
        return emptyFDGraph();
    }

    // PVTNUM is traditional one-based region identifier.  Subtract one to
    // form valid index into std::vector<>s.
    const auto regID = this->pvtnum_[activeCell] - 1;

    if (phase == ECLPhaseIndex::Liquid) {
        // Caller requests oil properties.
        return rawPvtCurve(this->oil_, curve, regID, *this->usys_);
    }

    // Call requests gas properties.
    return rawPvtCurve(this->gas_, curve, regID, *this->usys_);
}

std::vector<double>
Opm::ECLPVT::ECLPvtCurveCollection::
getDynamicProperty(const RawCurve             property,
                   const ECLPhaseIndex        phase,
                   const int                  activeCell,
                   const std::vector<double>& phasePress,
                   const std::vector<double>& mixRatio) const
{
    if (phase == ECLPhaseIndex::Aqua) {
        // Not supported at this time.
        // Return empty.
        return {};
    }

    if (property == RawCurve::SaturatedState) {
        // Not a supported request.  Return empty.
        return {};
    }

    if (static_cast<decltype(this->pvtnum_.size())>(activeCell)
        >= this->pvtnum_.size())
    {
        // Active cell index out of bounds.  Return empty.
        return {};
    }

    // PVTNUM is traditional one-based region identifier.  Subtract one to
    // form valid index into std::vector<>s.
    const auto regID = this->pvtnum_[activeCell] - 1;

    if (phase == ECLPhaseIndex::Liquid) {
        // Caller requests oil properties.
        return oilProperty(this->oil_, property, regID,
                           phasePress, mixRatio, *this->usys_);
    }

    // Call requests gas properties.
    return gasProperty(this->gas_, property, regID,
                       phasePress, mixRatio, *this->usys_);
}
