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

#ifndef OPM_ECLSATURATIONFUNC_HEADER_INCLUDED
#define OPM_ECLSATURATIONFUNC_HEADER_INCLUDED

#include <opm/utility/ECLPhaseIndex.hpp>

#include <memory>
#include <vector>

/// \file
///
/// Public interface to relative permeability evaluation machinery.  The
/// back-end is aware of ECLIPSE's standard three-phase model for relative
/// permeability of oil and the two- and three-point saturation end-point
/// scaling options.  Vertical scaling of relative permeability is not
/// supported at present.

namespace Opm {
    class ECLGraph;
    class ECLRestartData;
    class ECLInitFileData;

    /// Gateway to engine for computing relative permeability values based
    /// on tabulated saturation functions in ECL output.
    class ECLSaturationFunc
    {
    public:
        /// Constructor
        ///
        /// \param[in] G Connected topology of current model's active cells.
        ///    Needed to linearise region mapping (e.g., SATNUM) that is
        ///    distributed on local grids to all of the model's active cells
        ///    (\code member function G.rawLinearisedCellData() \endcode).
        ///
        /// \param[in] init Container of tabulated saturation functions and
        ///    saturation table end points, if applicable, for all active
        ///    cells in the model \p G.
        ///
        /// \param[in] useEPS Whether or not to include effects of
        ///    saturation end-point scaling.  No effect if the INIT result
        ///    set does not actually include saturation end-point scaling
        ///    data.  Otherwise, enables turning EPS off even if associate
        ///    data is present in the INIT result set.
        ///
        ///    Default value (\c true) means that effects of EPS are
        ///    included if requisite data is present in the INIT result.
        ECLSaturationFunc(const ECLGraph&        G,
                          const ECLInitFileData& init,
                          const bool             useEPS = true);

        /// Destructor.
        ~ECLSaturationFunc();

        /// Move constructor.
        ///
        /// Subsumes the implementation of an existing object.
        ///
        /// \param[in] rhs Existing engine for saturation function
        ///    evaluation.  Does not have a valid implementation when the
        ///    constructor completes.
        ECLSaturationFunc(ECLSaturationFunc&& rhs);

        /// Copy constructor.
        ///
        /// \param[in] rhs Existing engine for saturation function
        ///    evaluation.
        ECLSaturationFunc(const ECLSaturationFunc& rhs);

        /// Move assignment operator.
        ///
        /// Subsumes the implementation of an existing object.
        ///
        /// \param[in] rhs Existing engine for saturation function
        ///    evaluation.  Does not have a valid implementation when the
        ///    constructor completes.
        ///
        /// \return \code *this \endcode.
        ECLSaturationFunc& operator=(ECLSaturationFunc&& rhs);

        /// Assignment operator.
        ///
        /// \param[in] rhs Existing engine for saturation function
        ///    evaluation.
        ///
        /// \return \code *this \endcode.
        ECLSaturationFunc& operator=(const ECLSaturationFunc& rhs);

        /// Compute relative permeability values in all active cells for a
        /// single phase.
        ///
        /// \param[in] G Connected topology of current model's active cells.
        ///    Needed to linearise phase saturations (e.g., SOIL) that are
        ///    distributed on local grids to all of the model's active cells
        ///    (\code member function G.rawLinearisedCellData() \endcode).
        ///
        /// \param[in] rstrt ECLIPSE restart vectors.  Result set view
        ///    assumed to be positioned at a particular report step of
        ///    interest.
        ///
        /// \param[in] p Phase for which to compute relative permeability
        ///    values.
        ///
        /// \return Derived relative permeability values of active phase \p
        ///    p for all active cells in model \p G.  Empty if phase \p p is
        ///    not actually active in the current result set.
        std::vector<double>
        relperm(const ECLGraph&       G,
                const ECLRestartData& rstrt,
                const ECLPhaseIndex   p) const;

    private:
        /// Implementation backend.
        class Impl;

        /// Pointer to actual backend/engine object.
        std::unique_ptr<Impl> pImpl_;
    };

} // namespace Opm

#endif // OPM_ECLSATURATIONFUNC_HEADER_INCLUDED
