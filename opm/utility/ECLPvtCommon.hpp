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

#ifndef OPM_ECLPVTCOMMON_HEADER_INCLUDED
#define OPM_ECLPVTCOMMON_HEADER_INCLUDED

#include <opm/utility/ECLPiecewiseLinearInterpolant.hpp>
#include <opm/utility/ECLPropTable.hpp>
#include <opm/utility/ECLTableInterpolation1D.hpp>
#include <opm/utility/ECLUnitHandling.hpp>

#include <functional>
#include <memory>
#include <utility>
#include <type_traits>

/// \file
///
/// Facility for evaluating pressure-dependent fluid properties (formation
/// volume factor, viscosities &c) for oil or gas based on tabulated
/// descriptions as represented in an ECL result set (INIT file 'TAB'
/// vector).

namespace Opm { namespace ECLPVT {

    /// Protocol for converting raw table input data to strict SI unit
    /// conventions.
    struct ConvertUnits
    {
        /// Convenience type alias for a value transformation.
        using Converter = std::function<double(const double)>;

        /// How to convert the independent variate (1st column)
        Converter indep;

        /// How to convert the dependent variates (2nd... columns).
        std::vector<Converter> column;
    };

    /// Collection of unit converters for PVT quantities tabulated in a
    /// result set's INIT file.
    struct CreateUnitConverter {
        /// Convert quantities from native representations to strict SI
        /// units of measure.
        struct ToSI {
            /// Convert quantities of type pressure to strict SI units of
            /// measure (i.e., to Pascal units).
            ///
            /// \param[in] usys Native unit system for particular result
            ///    set.
            ///
            /// \return Value transformation function affecting requisite
            ///    unit conversion.
            static ConvertUnits::Converter
            pressure(const ::Opm::ECLUnits::UnitSystem& usys);

            /// Convert quantities of type dissolved gas-oil ratio (Rs) to
            /// strict SI units of measure (i.e., to Sm^3/Sm^3).
            ///
            /// \param[in] usys Native unit system for particular result
            ///    set.
            ///
            /// \return Value transformation function affecting requisite
            ///    unit conversion.
            static ConvertUnits::Converter
            disGas(const ::Opm::ECLUnits::UnitSystem& usys);

            /// Convert quantities of type vaporised oil-gas ratio (Rv) to
            /// strict SI units of measure (i.e., to Sm^3/Sm^3).
            ///
            /// \param[in] usys Native unit system for particular result
            ///    set.
            ///
            /// \return Value transformation function affecting requisite
            ///    unit conversion.
            static ConvertUnits::Converter
            vapOil(const ::Opm::ECLUnits::UnitSystem& usys);

            /// Convert quantities of type reciprocal formation volume
            /// factor (1/B) to strict SI units of measure (i.e., to
            /// Sm^3/Rm^3).
            ///
            /// \param[in] usys Native unit system for particular result
            ///    set.
            ///
            /// \return Value transformation function affecting requisite
            ///    unit conversion.
            static ConvertUnits::Converter
            recipFvf(const ::Opm::ECLUnits::UnitSystem& usys);

            /// Convert derivatives of quantities of type reciprocal
            /// formation volume factor (1/B) with respect to fluid pressure
            /// to strict SI units of measure (i.e., to Sm^3/(Rm^3 * Pa)).
            ///
            /// \param[in] usys Native unit system for particular result
            ///    set.
            ///
            /// \return Value transformation function affecting requisite
            ///    unit conversion.
            static ConvertUnits::Converter
            recipFvfDerivPress(const ::Opm::ECLUnits::UnitSystem& usys);

            /// Convert derivatives of quantities of type reciprocal
            /// formation volume factor (1/B) with respect to vaporised
            /// oil-gas ratio to strict SI units of measure (i.e., to
            /// Sm^3/(Rm^3 * (Sm^3 / Sm^3))).
            ///
            /// \param[in] usys Native unit system for particular result
            ///    set.
            ///
            /// \return Value transformation function affecting requisite
            ///    unit conversion.
            static ConvertUnits::Converter
            recipFvfDerivVapOil(const ::Opm::ECLUnits::UnitSystem& usys);

            /// Convert quantities of type reciprocal product of formation
            /// volume factor and phase viscosity (1/(B * mu)) to strict SI
            /// units of measure (i.e., to Sm^3/(Rm^3 * Pa*s)).
            ///
            /// \param[in] usys Native unit system for particular result
            ///    set.
            ///
            /// \return Value transformation function affecting requisite
            ///    unit conversion.
            static ConvertUnits::Converter
            recipFvfVisc(const ::Opm::ECLUnits::UnitSystem& usys);

            /// Convert derivatives of quantities of type reciprocal product
            /// of formation volume factor and phase viscosity (1/(B * mu))
            /// with respect to fluid pressure to strict SI units of measure
            /// (i.e., to Sm^3/(Rm^3 * Pa*s * Pa)).
            ///
            /// \param[in] usys Native unit system for particular result
            ///    set.
            ///
            /// \return Value transformation function affecting requisite
            ///    unit conversion.
            static ConvertUnits::Converter
            recipFvfViscDerivPress(const ::Opm::ECLUnits::UnitSystem& usys);

            /// Convert derivatives of quantities of type reciprocal product
            /// of formation volume factor and phase viscosity (1/(B * mu))
            /// with respect to vaporised oil-gas ratio to strict SI units
            /// of measure (i.e., to Sm^3/(Rm^3 * Pa*s * (Sm^3/Sm^3))).
            ///
            /// \param[in] usys Native unit system for particular result
            ///    set.
            ///
            /// \return Value transformation function affecting requisite
            ///    unit conversion.
            static ConvertUnits::Converter
            recipFvfViscDerivVapOil(const ::Opm::ECLUnits::UnitSystem& usys);

            /// Convert quantities of type reciprocal formation volume
            /// factor (1/B) to strict SI units of measure (i.e., to
            /// Sm^3/Rm^3).
            ///
            /// Specialisation for Gas.
            ///
            /// \param[in] usys Native unit system for particular result
            ///    set.
            ///
            /// \return Value transformation function affecting requisite
            ///    unit conversion.
            static ConvertUnits::Converter
            recipFvfGas(const ::Opm::ECLUnits::UnitSystem& usys);

            /// Convert derivatives of quantities of type reciprocal
            /// formation volume factor (1/B) with respect to fluid pressure
            /// to strict SI units of measure (i.e., to Sm^3/(Rm^3 * Pa)).
            ///
            /// Specialisation for Gas.
            ///
            /// \param[in] usys Native unit system for particular result
            ///    set.
            ///
            /// \return Value transformation function affecting requisite
            ///    unit conversion.
            static ConvertUnits::Converter
            recipFvfGasDerivPress(const ::Opm::ECLUnits::UnitSystem& usys);

            /// Convert derivatives of quantities of type reciprocal
            /// formation volume factor (1/B) with respect to vaporised
            /// oil-gas ratio to strict SI units of measure (i.e., to
            /// Sm^3/(Rm^3 * (Sm^3 / Sm^3))).
            ///
            /// Specialisation for Gas.
            ///
            /// \param[in] usys Native unit system for particular result
            ///    set.
            ///
            /// \return Value transformation function affecting requisite
            ///    unit conversion.
            static ConvertUnits::Converter
            recipFvfGasDerivVapOil(const ::Opm::ECLUnits::UnitSystem& usys);

            /// Convert quantities of type reciprocal product of formation
            /// volume factor and phase viscosity (1/(B * mu)) to strict SI
            /// units of measure (i.e., to Sm^3/(Rm^3 * Pa*s)).
            ///
            /// Specialisation for Gas.
            ///
            /// \param[in] usys Native unit system for particular result
            ///    set.
            ///
            /// \return Value transformation function affecting requisite
            ///    unit conversion.
            static ConvertUnits::Converter
            recipFvfGasVisc(const ::Opm::ECLUnits::UnitSystem& usys);

            /// Convert derivatives of quantities of type reciprocal product
            /// of formation volume factor and phase viscosity (1/(B * mu))
            /// with respect to fluid pressure to strict SI units of measure
            /// (i.e., to Sm^3/(Rm^3 * Pa*s * Pa)).
            ///
            /// Specialisation for Gas.
            ///
            /// \param[in] usys Native unit system for particular result
            ///    set.
            ///
            /// \return Value transformation function affecting requisite
            ///    unit conversion.
            static ConvertUnits::Converter
            recipFvfGasViscDerivPress(const ::Opm::ECLUnits::UnitSystem& usys);

            /// Convert derivatives of quantities of type reciprocal product
            /// of formation volume factor and phase viscosity (1/(B * mu))
            /// with respect to vaporised oil-gas ratio to strict SI units
            /// of measure (i.e., to Sm^3/(Rm^3 * Pa*s * (Sm^3/Sm^3))).
            ///
            /// Specialisation for Gas.
            ///
            /// \param[in] usys Native unit system for particular result
            ///    set.
            ///
            /// \return Value transformation function affecting requisite
            ///    unit conversion.
            static ConvertUnits::Converter
            recipFvfGasViscDerivVapOil(const ::Opm::ECLUnits::UnitSystem& usys);
        };
    };

    /// Evaluate pressure-dependent properties (formation volume factor,
    /// viscosity &c) for dead oil (PVDO) or dry gas (PVDG) from tabulated
    /// functions as represented in an ECL result set (ECLInitData).
    class PVDx
    {
    public:
        /// Convenience type alias.
        using ElemIt = ECLPropTableRawData::ElementIterator;

        /// Constructor.
        ///
        /// \param[in] xBegin Starting position of range of independent
        ///    variable (phase pressure).
        ///
        /// \param[in] xEnd One past the end of range of independent
        ///    variable.  Must be reachable from \p xBegin.
        ///
        /// \param[in,out] colIt Column iterators that reference the
        ///    dependent variables (reciprocal FVF &c).  Should be size 4 to
        ///    represent the FVF, the viscosity and the derivatives with
        ///    respect to phase pressure.  On input, positioned at the
        ///    beginning of a single table's dependent variables columns.
        ///    On output, advanced across \code std::distance(xBegin, xEnd)
        ///    \endcode rows/entries.
        PVDx(ElemIt               xBegin,
             ElemIt               xEnd,
             const ConvertUnits&  convert,
             std::vector<ElemIt>& colIt);

        /// Evaluate the phase FVF in selection of pressure points.
        ///
        /// \param[in] p Set of phase pressure points.
        ///
        /// \return Phase Formation volume factors for each pressure point.
        std::vector<double>
        formationVolumeFactor(const std::vector<double>& p) const;

        /// Evaluate the phase viscosity in selection of pressure points.
        ///
        /// \param[in] p Set of phase pressure points.
        ///
        /// \return Phase viscosity for each pressure point.
        std::vector<double>
        viscosity(const std::vector<double>& p) const;

    private:
        /// Extrapolation policy for property evaluator/interpolant.
        using Extrap = ::Opm::Interp1D::PiecewisePolynomial::
            ExtrapolationPolicy::LinearlyWithDerivatives;

        /// Type of fundamental table interpolant.
        using Backend = ::Opm::Interp1D::PiecewisePolynomial::Linear<Extrap>;

        /// Convenience type alias representing the interpolant's evaluation
        /// point conventions.
        using EvalPt = std::decay<
            decltype(std::declval<Backend>().classifyPoint(0.0))
        >::type;

        /// Interpolant for table of
        ///
        ///    [1/B, 1/(B*mu), d(1/B)/dp, d(1/(B*mu))/dp]
        ///
        /// versus "pressure" p--typically phase pressure for oil (Po) or
        /// gas (Pg).
        Backend interp_;

        /// Translate pressure value to evaluation point and identify
        /// relevant extrapolation case if needed.
        EvalPt getInterpPoint(const double p) const
        {
            return this->interp_.classifyPoint(p);
        }

        /// Interpolate reciprocal FVF at evaluation point.
        double fvf_recip(const EvalPt& pt) const
        {
            const auto col = std::size_t{0};

            return this->interp_.evaluate(col, pt);
        }

        /// Interpolate reciprocal product of FVF and viscosity at
        /// evaluation point.
        double fvf_mu_recip(const EvalPt& pt) const
        {
            const auto col = std::size_t{1};

            return this->interp_.evaluate(col, pt);
        }

        /// Compute a dynamic quantity such as the FVF or viscosity for a
        /// selection of pressure values.
        template <class EvalDynamicQuant>
        std::vector<double>
        computeQuantity(const std::vector<double>& p,
                        EvalDynamicQuant&&         eval) const
        {
            auto result = std::vector<double>{};
            result.reserve(p.size());

            for (const auto& pi : p) {
                result.push_back(eval(this->getInterpPoint(pi)));
            }

            return result;
        }
    };

    /// Evaluate pressure-dependent properties (formation volume factor,
    /// viscosity &c) for live oil (PVTO) or wet gas (PVTG) from tabulated
    /// functions as represented in an ECL result set (ECLInitData).
    template <class SubtableInterpolant>
    class PVTx
    {
    public:
        PVTx(std::vector<double>              key,
             std::vector<SubtableInterpolant> propInterp)
            : key_       (std::move(key))
            , propInterp_(std::move(propInterp))
        {
            if (this->key_.size() != this->propInterp_.size()) {
                throw std::invalid_argument {
                    "Size of Key Table Does Not Match "
                    "Number of Sub-Table Interpolants"
                };
            }

            if (this->key_.size() < 2) {
                throw std::invalid_argument {
                    "Mixing-Dependent Property Evaluator "
                    "Must Have At Least Two Inner Tables"
                };
            }
        }

        struct PrimaryKey
        {
            const std::vector<double>& data;
        };

        struct InnerVariate
        {
            const std::vector<double>& data;
        };

        std::vector<double>
        formationVolumeFactor(const PrimaryKey&   key,
                              const InnerVariate& x) const
        {
            const auto col = std::size_t{0};

            return this->computeQuantity(key, x,
                [this, col](const std::size_t     curve,
                            const InnerEvalPoint& pt) -> double
            {
                // 1 / (1 / B).
                return 1 / this->propInterp_[curve].evaluate(col, pt);
            });
        }

        std::vector<double>
        viscosity(const PrimaryKey& key,
                  const InnerVariate& x) const
        {
            return this->computeQuantity(key, x,
                [this](const std::size_t     curve,
                       const InnerEvalPoint& pt) -> double
            {
                const auto& I = this->propInterp_[curve];

                const auto fvf_recip    = I.evaluate(0, pt);
                const auto fvf_mu_recip = I.evaluate(1, pt);

                // (1 / B) / (1 / (B*mu))
                return fvf_recip / fvf_mu_recip;
            });
        }

    private:
        using InnerEvalPoint = typename std::decay<
            decltype(std::declval<SubtableInterpolant>().classifyPoint(0.0))
        >::type;

        using OuterInterpPoint =
            ::Opm::Interp1D::PiecewisePolynomial::LocalInterpPoint;

        struct InnerInterpPoint {
            InnerEvalPoint left;
            InnerEvalPoint right;
        };

        std::vector<double> key_;
        std::vector<SubtableInterpolant> propInterp_;

        InnerInterpPoint getInterpPoint(const std::size_t i,
                                        const double      x) const
        {
            assert ((i + 1) < this->propInterp_.size());

            return {
                this->propInterp_[i + 0].classifyPoint(x),
                this->propInterp_[i + 1].classifyPoint(x)
            };
        }

        template <class Function>
        double interpolate(Function&&              func,
                           const OuterInterpPoint& outer,
                           const double            x) const
        {
            assert (outer.cat == ::Opm::Interp1D::PointCategory::InRange);

            const auto pt =
                this->getInterpPoint(outer.interval, x);

            const auto yLeft  = func(outer.interval + 0, pt.left);
            const auto yRight = func(outer.interval + 1, pt.right);

            const auto t = outer.t / (this->key_[outer.interval + 1] -
                                      this->key_[outer.interval + 0]);

            // t == 0 => yLeft, t == 1 => yRight
            return t*yRight + (1.0 - t)*yLeft;
        }

        template <class Function>
        double extrapLeft(Function&&              func,
                          const OuterInterpPoint& outer,
                          const double            x) const
        {
            assert (outer.cat == ::Opm::Interp1D::PointCategory::LeftOfRange);
            assert (outer.interval == 0*this->key_.size());

            const auto pt =
                this->getInterpPoint(outer.interval, x);

            const auto yLeft  = func(0, pt.left);
            const auto yRight = func(1, pt.right);

            const auto dydk =
                (yRight - yLeft) / (this->key_[1] - this->key_[0]);

            return yLeft + outer.t*dydk;
        }

        template <class Function>
        double extrapRight(Function&&              func,
                           const OuterInterpPoint& outer,
                           const double            x) const
        {
            const auto nIntervals = this->key_.size() - 1;

            assert (outer.cat == ::Opm::Interp1D::PointCategory::RightOfRange);
            assert (outer.interval == nIntervals);

            const auto pt = this->getInterpPoint(nIntervals - 1, x);

            const auto yLeft  = func(nIntervals - 1, pt.left);
            const auto yRight = func(nIntervals - 0, pt.right);

            const auto dydk =
                (yRight - yLeft) / (this->key_[nIntervals - 0] -
                                    this->key_[nIntervals - 1]);

            return yRight + outer.t*dydk;
        }

        template <class Function>
        double evaluate(const double key,
                        const double x,
                        Function&&   func) const
        {
            const auto outer = ::Opm::Interp1D::PiecewisePolynomial::
                LocalInterpPoint::identify(this->key_, key);

            switch (outer.cat) {
            case ::Opm::Interp1D::PointCategory::InRange:
                return this->interpolate(std::forward<Function>(func),
                                         outer, x);

            case ::Opm::Interp1D::PointCategory::LeftOfRange:
                return this->extrapLeft(std::forward<Function>(func),
                                        outer, x);

            case ::Opm::Interp1D::PointCategory::RightOfRange:
                return this->extrapRight(std::forward<Function>(func),
                                         outer, x);
            }

            throw std::logic_error {
                "Outer/Primary Key Cannot Be Classified"
            };
        }

        template <class Function>
        std::vector<double>
        computeQuantity(const PrimaryKey&   key,
                        const InnerVariate& x,
                        Function            func) const
        {
            auto result = std::vector<double>{};

            const auto nVals = key.data.size();

            if (x.data.size() != nVals) {
                throw std::invalid_argument {
                    "Number of Inner Sampling Points Does Not Match "
                    "Number of Outer Sampling Points"
                };
            }

            result.reserve(nVals);

            for (auto i = 0*nVals; i < nVals; ++i) {
                const auto q =
                    this->evaluate(key.data[i], x.data[i], func);

                result.push_back(q);
            }

            return result;
        }
    };

}} // Opm::ECLPVT

#endif // OPM_ECLPVTCOMMON_HEADER_INCLUDED
