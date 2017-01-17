/*
  Copyright 2017 SINTEF ICT, Applied Mathematics.
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

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <opm/utility/ECLUnitHandling.hpp>

#include <opm/parser/eclipse/Units/Units.hpp>

#include <exception>
#include <stdexcept>
#include <string>

#include <ert/ecl/ecl_util.h>

namespace Opm { namespace ECLUnits {
    namespace Impl {
        ert_ecl_unit_enum getUnitConvention(const int usys);

        template <ert_ecl_unit_enum convention>
        class USys;

        template <>
        class USys<ECL_METRIC_UNITS> : public ::Opm::ECLUnits::UnitSystem
        {
        public:
            virtual double pressure() const override
            {
                return Metric::Pressure;
            }

            virtual double reservoirRate() const override
            {
                return Metric::ReservoirVolume / Metric::Time;
            }

            virtual double reservoirVolume() const override
            {
                return Metric::ReservoirVolume;
            }

            virtual double time() const override
            {
                return Metric::Time;
            }

            virtual double transmissibility() const override
            {
                return Metric::Transmissibility;
            }
        };

        template <>
        class USys<ECL_FIELD_UNITS> : public ::Opm::ECLUnits::UnitSystem
        {
        public:
            virtual double pressure() const override
            {
                return Field::Pressure;
            }

            virtual double reservoirRate() const override
            {
                return Field::ReservoirVolume / Field::Time;
            }

            virtual double reservoirVolume() const override
            {
                return Field::ReservoirVolume;
            }

            virtual double time() const override
            {
                return Field::Time;
            }

            virtual double transmissibility() const override
            {
                return Field::Transmissibility;
            }
        };

        template <>
        class USys<ECL_LAB_UNITS> : public ::Opm::ECLUnits::UnitSystem
        {
        public:
            virtual double pressure() const override
            {
                return Lab::Pressure;
            }

            virtual double reservoirRate() const override
            {
                return Lab::ReservoirVolume / Lab::Time;
            }

            virtual double reservoirVolume() const override
            {
                return Lab::ReservoirVolume;
            }

            virtual double time() const override
            {
                return Lab::Time;
            }

            virtual double transmissibility() const override
            {
                return Lab::Transmissibility;
            }
        };

        template <>
        class USys<ECL_PVT_M_UNITS> : public ::Opm::ECLUnits::UnitSystem
        {
        public:
            virtual double pressure() const override
            {
                return unit::atm;
            }

            virtual double reservoirRate() const override
            {
                using namespace prefix;
                using namespace unit;

                return cubic(meter) / day;
            }

            virtual double reservoirVolume() const override
            {
                using namespace prefix;
                using namespace unit;

                return cubic(meter);
            }

            virtual double time() const override
            {
                return unit::day;
            }

            virtual double transmissibility() const override
            {
                using namespace prefix;
                using namespace unit;

                return centi*Poise * cubic(meter) / (day * atm);
            }
        };
    } // namespace Impl
}} // namespace Opm::ECLUnits

ert_ecl_unit_enum
Opm::ECLUnits::Impl::getUnitConvention(const int usys)
{
    switch (usys) {
    case 1: return ECL_METRIC_UNITS;
    case 2: return ECL_FIELD_UNITS;
    case 3: return ECL_LAB_UNITS;
    case 4: return ECL_PVT_M_UNITS;
    }

    throw std::runtime_error("Unsupported Unit Convention: "
                             + std::to_string(usys));
}

std::unique_ptr<const ::Opm::ECLUnits::UnitSystem>
Opm::ECLUnits::createUnitSystem(const int usys)
{
    using UPtr =
        std::unique_ptr<const ::Opm::ECLUnits::UnitSystem>;

    switch (Impl::getUnitConvention(usys)) {
    case ECL_METRIC_UNITS:
        return UPtr{ new Impl::USys<ECL_METRIC_UNITS>{} };

    case ECL_FIELD_UNITS:
        return UPtr{ new Impl::USys<ECL_FIELD_UNITS>{} };

    case ECL_LAB_UNITS:
        return UPtr{ new Impl::USys<ECL_LAB_UNITS>{} };

    case ECL_PVT_M_UNITS:
        return UPtr{ new Impl::USys<ECL_PVT_M_UNITS>{} };
    }

    throw std::runtime_error("Unsupported Unit Convention: "
                             + std::to_string(usys));
}
