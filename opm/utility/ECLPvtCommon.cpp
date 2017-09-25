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

#include <opm/utility/ECLPvtCommon.hpp>

#include <opm/utility/ECLUnitHandling.hpp>

#include <opm/parser/eclipse/Units/Units.hpp>

#include <functional>

namespace {
    double fvfScale(const ::Opm::ECLUnits::UnitSystem& usys)
    {
        // B = [rVolume / sVolume(Liquid)]
        return usys.reservoirVolume()
            /  usys.surfaceVolumeLiquid();
    }

    double fvfGasScale(const ::Opm::ECLUnits::UnitSystem& usys)
    {
        // B = [rVolume / sVolume(Gas)]
        return usys.reservoirVolume()
            /  usys.surfaceVolumeGas();
    }

    double rsScale(const ::Opm::ECLUnits::UnitSystem& usys)
    {
        // Rs = [sVolume(Gas) / sVolume(Liquid)]
        return usys.surfaceVolumeGas()
            /  usys.surfaceVolumeLiquid();
    }

    double rvScale(const ::Opm::ECLUnits::UnitSystem& usys)
    {
        // Rv = [sVolume(Liq) / sVolume(Gas)]
        return usys.surfaceVolumeLiquid()
            /  usys.surfaceVolumeGas();
    }

    ::Opm::ECLPVT::ConvertUnits::Converter
    createConverterToSI(const double uscale)
    {
        return ::Opm::ECLPVT::ConvertUnits::Converter {
            [uscale](const double q) -> double
            {
                return ::Opm::unit::convert::from(q, uscale);
            }
        };
    }
}

Opm::ECLPVT::ConvertUnits::Converter
Opm::ECLPVT::CreateUnitConverter::ToSI::
pressure(const ::Opm::ECLUnits::UnitSystem& usys)
{
    return createConverterToSI(usys.pressure());
}

::Opm::ECLPVT::ConvertUnits::Converter
Opm::ECLPVT::CreateUnitConverter::ToSI::
compressibility(const ::Opm::ECLUnits::UnitSystem& usys)
{
    return createConverterToSI(1.0 / usys.pressure());
}

Opm::ECLPVT::ConvertUnits::Converter
Opm::ECLPVT::CreateUnitConverter::ToSI::
disGas(const ::Opm::ECLUnits::UnitSystem& usys)
{
    return createConverterToSI(rsScale(usys));
}

Opm::ECLPVT::ConvertUnits::Converter
Opm::ECLPVT::CreateUnitConverter::ToSI::
vapOil(const ::Opm::ECLUnits::UnitSystem& usys)
{
    return createConverterToSI(rvScale(usys));
}

Opm::ECLPVT::ConvertUnits::Converter
Opm::ECLPVT::CreateUnitConverter::ToSI::
recipFvf(const ::Opm::ECLUnits::UnitSystem& usys)
{
    return createConverterToSI(1.0 / fvfScale(usys));
}

Opm::ECLPVT::ConvertUnits::Converter
Opm::ECLPVT::CreateUnitConverter::ToSI::
recipFvfDerivPress(const ::Opm::ECLUnits::UnitSystem& usys)
{
    // d(1/B)/dp
    const auto B_scale = fvfScale(usys);
    const auto P_scale = usys.pressure();

    return createConverterToSI(1.0 / (B_scale * P_scale));
}

Opm::ECLPVT::ConvertUnits::Converter
Opm::ECLPVT::CreateUnitConverter::ToSI::
recipFvfDerivVapOil(const ::Opm::ECLUnits::UnitSystem& usys)
{
    // d(1/B)/dRv
    const auto B_scale  = fvfScale(usys);
    const auto Rv_scale = rvScale(usys);

    return createConverterToSI(1.0 / (B_scale * Rv_scale));
}

Opm::ECLPVT::ConvertUnits::Converter
Opm::ECLPVT::CreateUnitConverter::ToSI::
recipFvfVisc(const ::Opm::ECLUnits::UnitSystem& usys)
{
    const auto Bscale     = fvfScale(usys);
    const auto visc_scale = usys.viscosity();

    return createConverterToSI(1.0 / (Bscale * visc_scale));
}

Opm::ECLPVT::ConvertUnits::Converter
Opm::ECLPVT::CreateUnitConverter::ToSI::
recipFvfViscDerivPress(const ::Opm::ECLUnits::UnitSystem& usys)
{
    // d(1/(B*mu))/dp
    const auto B_scale  = fvfScale(usys);
    const auto P_scale  = usys.pressure();
    const auto mu_scale = usys.viscosity();

    return createConverterToSI(1.0 / (B_scale * mu_scale * P_scale));
}

::Opm::ECLPVT::ConvertUnits::Converter
Opm::ECLPVT::CreateUnitConverter::ToSI::
recipFvfViscDerivVapOil(const ::Opm::ECLUnits::UnitSystem& usys)
{
    // d(1/(B*mu))/dRv
    const auto B_scale  = fvfScale(usys);
    const auto mu_scale = usys.viscosity();
    const auto Rv_scale = rvScale(usys);

    return createConverterToSI(1.0 / (B_scale * mu_scale * Rv_scale));
}

Opm::ECLPVT::ConvertUnits::Converter
Opm::ECLPVT::CreateUnitConverter::ToSI::
recipFvfGas(const ::Opm::ECLUnits::UnitSystem& usys)
{
    return createConverterToSI(1.0 / fvfGasScale(usys));
}

Opm::ECLPVT::ConvertUnits::Converter
Opm::ECLPVT::CreateUnitConverter::ToSI::
recipFvfGasDerivPress(const ::Opm::ECLUnits::UnitSystem& usys)
{
    // d(1/B)/dp
    const auto B_scale = fvfGasScale(usys);
    const auto P_scale = usys.pressure();

    return createConverterToSI(1.0 / (B_scale * P_scale));
}

Opm::ECLPVT::ConvertUnits::Converter
Opm::ECLPVT::CreateUnitConverter::ToSI::
recipFvfGasDerivVapOil(const ::Opm::ECLUnits::UnitSystem& usys)
{
    // d(1/B)/dRv
    const auto B_scale  = fvfGasScale(usys);
    const auto Rv_scale = rvScale(usys);

    return createConverterToSI(1.0 / (B_scale * Rv_scale));
}

Opm::ECLPVT::ConvertUnits::Converter
Opm::ECLPVT::CreateUnitConverter::ToSI::
recipFvfGasVisc(const ::Opm::ECLUnits::UnitSystem& usys)
{
    const auto Bscale     = fvfGasScale(usys);
    const auto visc_scale = usys.viscosity();

    return createConverterToSI(1.0 / (Bscale * visc_scale));
}

Opm::ECLPVT::ConvertUnits::Converter
Opm::ECLPVT::CreateUnitConverter::ToSI::
recipFvfGasViscDerivPress(const ::Opm::ECLUnits::UnitSystem& usys)
{
    // d(1/(B*mu))/dp
    const auto B_scale  = fvfGasScale(usys);
    const auto P_scale  = usys.pressure();
    const auto mu_scale = usys.viscosity();

    return createConverterToSI(1.0 / (B_scale * mu_scale * P_scale));
}

::Opm::ECLPVT::ConvertUnits::Converter
Opm::ECLPVT::CreateUnitConverter::ToSI::
recipFvfGasViscDerivVapOil(const ::Opm::ECLUnits::UnitSystem& usys)
{
    // d(1/(B*mu))/dRv
    const auto B_scale  = fvfGasScale(usys);
    const auto mu_scale = usys.viscosity();
    const auto Rv_scale = rvScale(usys);

    return createConverterToSI(1.0 / (B_scale * mu_scale * Rv_scale));
}

// =====================================================================

Opm::ECLPVT::PVDx::PVDx(ElemIt               xBegin,
                        ElemIt               xEnd,
                        const ConvertUnits&  convert,
                        std::vector<ElemIt>& colIt)
    : interp_(Extrap{}, xBegin, xEnd, colIt,
              convert.indep, convert.column)
{}

std::vector<double>
Opm::ECLPVT::PVDx::formationVolumeFactor(const std::vector<double>& p) const
{
    return this->computeQuantity(p,
        [this](const EvalPt& pt) -> double
    {
        // 1 / (1 / B)
        return 1.0 / this->fvf_recip(pt);
    });
}

std::vector<double>
Opm::ECLPVT::PVDx::viscosity(const std::vector<double>& p) const
{
    return this->computeQuantity(p,
        [this](const EvalPt& pt) -> double
    {
        // (1 / B) / (1 / (B * mu)
        return this->fvf_recip(pt) / this->fvf_mu_recip(pt);
    });
}
