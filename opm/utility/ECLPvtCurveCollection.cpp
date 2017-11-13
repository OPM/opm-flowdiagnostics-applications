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

#include <vector>

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
}

Opm::ECLPVT::ECLPvtCurveCollection::
ECLPvtCurveCollection(const ECLGraph&        G,
                      const ECLInitFileData& init)
    : pvtnum_(pvtnumVector(G, init))
    , gas_   (CreateGasPVTInterpolant::fromECLOutput(init)) // u_p<> -> s_p<>
    , oil_   (CreateOilPVTInterpolant::fromECLOutput(init)) // u_p<> -> s_p<>
{}

Opm::FlowDiagnostics::Graph
Opm::ECLPVT::ECLPvtCurveCollection::
getPvtCurve(const RawCurve      curve,
            const ECLPhaseIndex phase,
            const int           activeCell) const
{
    if (phase == ECLPhaseIndex::Aqua) {
        // Not supported at this time.
        // Return empty.
        return FlowDiagnostics::Graph{};
    }

    if (static_cast<decltype(this->pvtnum_.size())>(activeCell)
        >= this->pvtnum_.size())
    {
        // Active cell index out of bounds.  Return empty.
        return FlowDiagnostics::Graph{};
    }

    // PVTNUM is traditional one-based region identifier.  Subtract one to
    // form valid index into std::vector<>s.
    const auto regID = this->pvtnum_[activeCell] - 1;

    if (phase == ECLPhaseIndex::Liquid) {
        // return this->oil_->getPvtCurve(curve, regID);
        return FlowDiagnostics::Graph{};
    }

    if (this->gas_) {
        return this->gas_->getPvtCurve(curve, regID);
    }
    else {
        // Result set does not provide tabulated gas properties.  Return
        // empty.
        return FlowDiagnostics::Graph{};
    }
}
