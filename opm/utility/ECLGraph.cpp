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

#include <opm/utility/ECLGraph.hpp>
#include <opm/utility/ECLUtilities.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <exception>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <boost/filesystem.hpp>

#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_grid.h>
#include <ert/ecl/ecl_kw_magic.h>
#include <ert/ecl/ecl_nnc_export.h>
#include <ert/util/ert_unique_ptr.hpp>

/// \file
///
/// Implementation of \c ECLGraph interface.


/// Implementation of ECLGraph interface.
class Opm::ECLGraph::Impl
{
public:
    /// Constructor
    ///
    /// \param[in] grid Name or prefix of ECL grid (i.e., .GRID or
    ///                 .EGRID) file.
    ///
    /// \param[in] init Name of ECL INIT file corresponding to \p grid
    ///                 input.  Assumed to provide at least a complete set
    ///                 of pore-volume values (i.e., for all global cells
    ///                 defined in the \p grid).
    ///
    ///                 If available in the INIT file, the constructor will
    ///                 also leverage the transmissibility data when
    ///                 constructing the active cell neighbourship table.
    Impl(const Path& grid, const Path& init);

    /// Assign source object for phase flux calculation.
    ///
    /// \param[in] src Name of ECL restart file, possibly unified, from
    ///                which next set of phase fluxes should be retrieved.
    void assignDataSource(const Path& src);

    /// Retrieve number of active cells in graph.
    std::size_t numCells() const;

    /// Retrive number of connections in graph.
    std::size_t numConnections() const;

    /// Retrive neighbourship relations between active cells.
    ///
    /// The \c i-th connection is between active cells \code
    /// neighbours()[2*i + 0] \endcode and \code neighbours()[2*i + 1]
    /// \endcode.
    const std::vector<int>& neighbours() const;

    /// Retrive static pore-volume values on active cells only.
    ///
    /// Corresponds to the \c PORV vector in the INIT file, possibly
    /// restricted to those active cells for which the pore-volume is
    /// strictly positive.
    const std::vector<double>& activePoreVolume() const;

    /// Retrive phase flux on all connections defined by \code neighbours()
    /// \endcode.
    ///
    /// Non-"const" because this potentially loads new data from the backing
    /// store into internal cache data structures.
    ///
    /// \param[in] phase Canonical phase for which to retrive flux.
    ///
    /// \param[in] occurrence Selected temporal vector.  Essentially the
    ///                       report step number.
    ///
    /// \return Flux values corresponding to selected phase and report step.
    /// Empty if unavailable in the result set (e.g., by querying the gas
    /// flux in an oil/water system or if the specified \p occurrence is not
    /// reported due to report frequencies or no flux values are output at
    /// all).
    std::vector<double>
    flux(const BlackoilPhases::PhaseIndex phase,
         const int                        occurrence);

private:
    /// Collection of global (cell) IDs.
    using GlobalIDColl = std::vector<std::size_t>;

    /// Local convenience alias.
    using Direction = CellCollection::Direction;

    /// Collection of (global) cell IDs corresponding to the flow source of
    /// each connection.
    using OutCell = std::map<Direction, GlobalIDColl>;

    /// Handle to ECL result set.
    using CDataPtr = std::unique_ptr<ECL::GlobalCellData>;

    /// Flattened neighbourship relation (array of size \code
    /// 2*numConnections() \endcode).
    std::vector<int> neigh_;

    /// Static pore-volumes of all active cells.
    std::vector<double> activePVol_;

    /// Source cells for each Cartesian connection.
    OutCell outCell_;

    /// Subset of non-neighbouring connections that affect active cells in
    /// main grid.
    GlobalIDColl nncID_;

    /// Hook into ECL result set.  Pointer because we need deferred
    /// initialisation.
    CDataPtr data_;

    /// Derive connections in particular Cartesian direction on all active
    /// cells.
    ///
    /// Writes to \c neigh_ and \c outCell_.  Will possibly change the \code
    /// *data_ \endcode object.
    ///
    /// \param[in] d Cartesian direction in which to look for neighbouring
    ///              cells.
    ///
    /// \param[in] coll Backing data for neighbourship extraction.
    void deriveNeighbours(const CellCollection::Direction d,
                          const CellCollection&           coll);

    /// Extract explicit non-neighbouring connections from ECL output.
    ///
    /// Writes to \c neigh_ and \c nncID_.
    ///
    /// \param[in] G ERT Grid representation.
    ///
    /// \param[in] init ERT representation of INIT source.
    ///
    /// \param[in] coll Backing data for neighbourship extraction.
    void nnc(const ecl_grid_type*  G,
             const ecl_file_type*  init,
             const CellCollection& coll);

    /// Predicate that establishes availability of a particular phase flux
    /// in the current result set.
    ///
    /// \param[in] phase Canonical phase for which to check availability.
    ///
    /// \return Whether or not the current result set (\code *data_
    /// \endcode) contains phase flux data for the input phase.
    bool fluxAvailable(const BlackoilPhases::PhaseIndex phase) const;

    /// Compute ECL vector basename for particular phase flux.
    ///
    /// \param[in] phase Canonical phase for which to derive ECL vector
    /// basename.
    ///
    /// \return Basename for ECl vector corresponding to particular phase
    /// flux.
    std::string
    flowVector(const BlackoilPhases::PhaseIndex phase) const;

    /// Extract phase flux values for all Cartesian connections in a
    /// particular direction.
    ///
    /// Potentially modifies \code *data_ \endcode.
    ///
    /// \param[in] d Cartesian direction for which to extract phase flux
    ///              values.
    ///
    /// \param[in] vector Basename for ECL vector corresponding to
    ///              particular phase flux.  Intentially taken as mutable
    ///              copy which is modified in function body.
    ///
    /// \param[in] occurrence Selected temporal vector.  Essentially the
    ///              report step number.
    ///
    /// \param[in,out] flux Phase flux values.  Values corresponding to
    /// direction \c d are appended to the \p flux vector.
    void phaseFluxCartesian(const CellCollection::Direction d,
                            std::string                     vector,
                            const int                       occurrence,
                            std::vector<double>&            flux);

    /// Extract phase flux values for all non-neighbouring connections
    /// between active cells in main grid.
    ///
    /// Potentially modifies \code *data_ \endcode.
    ///
    /// \param[in] vector Basename for ECL vector corresponding to
    ///              particular phase flux.  Intentially taken as mutable
    ///              copy (modified in function body).
    ///
    /// \param[in] occurrence Selected temporal vector.  Essentially the
    ///              report step number.
    ///
    /// \param[in,out] flux Phase flux values.  Values corresponding to
    /// direction \c d are appended to the \p flux vector.
    void phaseFluxNNC(const std::string&   vector,
                      const int            occurrence,
                      std::vector<double>& flux);
};

Opm::ECLGraph::Impl::Impl(const Path& grid, const Path& init)
{
    const auto G = ECL::loadCase(grid);
    auto       I = ECL::loadFile(init);

    const auto pvol = ECL::getPVolVector(G.get(), I.get());

    this->data_.reset(new ECL::GlobalCellData(G.get(), pvol));
    this->data_->assignDataSource(std::move(I));

    const auto coll =
        CellCollection{ G.get(), this->data_->activeGlobalCells() };

    this->activePVol_.reserve(coll.activeGlobalCells().size());
    for (const auto& globCell : coll.activeGlobalCells()) {
        this->activePVol_.push_back(pvol[globCell]);
    }

    // Too large, but this is a quick estimate.
    this->neigh_.reserve(3 * (2 * this->numCells()));

    for (const auto d : { CellCollection::Direction::I ,
                          CellCollection::Direction::J ,
                          CellCollection::Direction::K })
    {
        this->deriveNeighbours(d, coll);
    }

    this->nnc(G.get(), this->data_->getDataSource(), coll);
}

void
Opm::ECLGraph::Impl::assignDataSource(const Path& src)
{
    this->data_->assignDataSource(src);
}

std::size_t
Opm::ECLGraph::Impl::numCells() const
{
    return this->activePVol_.size();
}

std::size_t
Opm::ECLGraph::Impl::numConnections() const
{
    assert (this->neighbours().size() % 2 == 0);

    return this->neighbours().size() / 2;
}

const std::vector<int>&
Opm::ECLGraph::Impl::neighbours() const
{
    return this->neigh_;
}

const std::vector<double>&
Opm::ECLGraph::Impl::activePoreVolume() const
{
    return this->activePVol_;
}

std::vector<double>
Opm::ECLGraph::Impl::
flux(const BlackoilPhases::PhaseIndex phase,
     const int                        occurrence)
{
    if (! this->fluxAvailable(phase)) {
        return {};
    }

    const auto vector = this->flowVector(phase);

    auto v = std::vector<double>{};

    for (const auto d : { CellCollection::Direction::I ,
                          CellCollection::Direction::J ,
                          CellCollection::Direction::K })
    {
        this->phaseFluxCartesian(d, vector, occurrence, v);
    }

    if (! this->nncID_.empty()) {
        this->phaseFluxNNC(vector, occurrence, v);
    }

    return v;
}

void
Opm::ECLGraph::Impl::
deriveNeighbours(const CellCollection::Direction d,
                 const CellCollection&           coll)
{
    auto tran = std::string{"TRAN"};

    switch (d) {
    case CellCollection::Direction::I:
        tran += 'X';
        break;

    case CellCollection::Direction::J:
        tran += 'Y';
        break;

    case CellCollection::Direction::K:
        tran += 'Z';
        break;

    default:
        throw std::invalid_argument("Input direction must be (I,J,K)");
    }

    const auto& T = this->data_->haveVector(tran)
        ? this->data_->getVector(tran)
        : std::vector<double>(coll.numGlobalCells(), 1.0);

    auto& ocell = this->outCell_[d];
    ocell.reserve(this->data_->activeGlobalCells().size());

    for (const auto& globID : coll.activeGlobalCells()) {
        if (T[globID] > 0.0) {
            const auto other = coll.cartesianNeighbour(globID, d);

            if (other >= 0) {
                const auto c1 = coll.activeCell(globID);

                assert (c1 != other);

                this->neigh_.push_back(c1);
                this->neigh_.push_back(other);

                ocell.push_back(globID);
            }
        }
    }
}

void
Opm::ECLGraph::Impl::nnc(const ecl_grid_type*  G,
                                  const ecl_file_type*  init,
                                  const CellCollection& coll)
{
    const auto nncData = ECL::loadNNC(G, init);
    const auto numNNC  = nncData.size();

    this->neigh_.reserve(this->neigh_.size() + (2 * numNNC));

    this->nncID_.clear();
    this->nncID_.reserve(numNNC);

    auto nncID = static_cast<std::size_t>(0);

    for (const auto& conn : nncData) {
        // Omit LGR.
        if ((conn.grid_nr1 == 0) && (conn.grid_nr2 == 0) &&
            (conn.trans > 0.0))
        {
            const auto c1 = coll.activeCell(conn.global_index1);
            const auto c2 = coll.activeCell(conn.global_index2);

            if ((c1 >= 0) && (c2 >= 0)) {
                assert (c1 != c2);

                this->neigh_.push_back(c1);
                this->neigh_.push_back(c2);

                this->nncID_.push_back(nncID);
            }
        }

        nncID += 1;
    }
}

bool
Opm::ECLGraph::Impl::
fluxAvailable(const BlackoilPhases::PhaseIndex phase) const
{
    const auto vector = this->flowVector(phase);

    auto haveVector = true;

    for (const auto* flowdir : { "I+", "J+", "K+" }) {
        haveVector = haveVector &&
            this->data_->haveVector(vector + flowdir);
    }

    if (haveVector && ! this->nncID_.empty()) {
        haveVector = this->data_->haveVector(vector + "N+");
    }

    return haveVector;
}

std::string
Opm::ECLGraph::Impl::
flowVector(const BlackoilPhases::PhaseIndex phase) const
{
    const auto vector = std::string("FLR"); // Flow-rate, reservoir

    if (phase == BlackoilPhases::PhaseIndex::Aqua) {
        return vector + "WAT";
    }

    if (phase == BlackoilPhases::PhaseIndex::Liquid) {
        return vector + "OIL";
    }

    if (phase == BlackoilPhases::PhaseIndex::Vapour) {
        return vector + "GAS";
    }

    {
        std::ostringstream os;

        os << "Invalid phase index '" << phase << '\'';

        throw std::invalid_argument(os.str());
    }
}

void
Opm::ECLGraph::Impl::
phaseFluxCartesian(const CellCollection::Direction d,
                   std::string                     vector,
                   const int                       occurrence,
                   std::vector<double>&            flux)
{
    switch (d) {
    case CellCollection::Direction::I:
        vector += 'I';
        break;

    case CellCollection::Direction::J:
        vector += 'J';
        break;

    case CellCollection::Direction::K:
        vector += 'K';
        break;
    }

    vector += '+';

    const auto& v = this->data_->getVector(vector, occurrence);

    flux.reserve(flux.size() + this->outCell_[d].size());

    for (const auto& globCell : this->outCell_[d]) {
        flux.push_back(v[globCell]);
    }
}

void
Opm::ECLGraph::Impl::
phaseFluxNNC(const std::string&   vector,
             const int            occurrence,
             std::vector<double>& flux)
{
    const auto& v =
        this->data_->getVector(vector + "N+", occurrence);

    flux.reserve(flux.size() + this->nncID_.size());

    for (const auto& nnc : this->nncID_) {
        flux.push_back(v[nnc]);
    }
}

// ======================================================================

Opm::ECLGraph::ECLGraph(ImplPtr pImpl)
    : pImpl_(std::move(pImpl))
{
}

Opm::ECLGraph::ECLGraph(ECLGraph&& rhs)
    : pImpl_(std::move(rhs.pImpl_))
{}

Opm::ECLGraph::~ECLGraph()
{}

Opm::ECLGraph&
Opm::ECLGraph::operator=(ECLGraph&& rhs)
{
    this->pImpl_ = std::move(rhs.pImpl_);

    return *this;
}

Opm::ECLGraph
Opm::ECLGraph::load(const Path& grid, const Path& init)
{
    auto pImpl = ImplPtr{new Impl(grid, init)};

    return { std::move(pImpl) };
}

void
Opm::ECLGraph::assignFluxDataSource(const Path& src)
{
    this->pImpl_->assignDataSource(src);
}

std::size_t
Opm::ECLGraph::numCells() const
{
    return this->pImpl_->numCells();
}

std::size_t
Opm::ECLGraph::numConnections() const
{
    return this->pImpl_->numConnections();
}

const std::vector<int>&
Opm::ECLGraph::neighbours() const
{
    return this->pImpl_->neighbours();
}

const std::vector<double>&
Opm::ECLGraph::poreVolume() const
{
    return this->pImpl_->activePoreVolume();
}

std::vector<double>
Opm::ECLGraph::
flux(const BlackoilPhases::PhaseIndex phase,
     const int                        occurrence)
{
    return this->pImpl_->flux(phase, occurrence);
}
