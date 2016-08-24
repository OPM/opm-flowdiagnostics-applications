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

namespace {
    namespace ECL {
        using GridPtr = ::ERT::ert_unique_ptr<ecl_grid_type, ecl_grid_free>;
        using FilePtr = ::ERT::ert_unique_ptr<ecl_file_type, ecl_file_close>;

        /// Retrieve global pore-volume vector from INIT source.
        ///
        /// Specialised tool needed to determine the active cells.
        ///
        /// \param[in] G ERT Grid representation.
        ///
        /// \param[in] init ERT representation of INIT source.
        ///
        /// \return Vector of pore-volumes for all global cells of \p G.
        std::vector<double>
        getPVolVector(const ecl_grid_type* G,
                      const ecl_file_type* init);

        /// Internalise on-disk representation of ECLIPSE grid.
        ///
        /// \param[in] grid Name or prefix of on-disk representation of
        ///                 ECLIPSE grid.  If using a prefix, the loader
        ///                 will consider both .EGRID and .GRID versions of
        ///                 the input.
        ///
        /// \return Internalised ERT Grid representation.
        GridPtr loadCase(const boost::filesystem::path& grid);

        /// Internalise on-disk representation of ECL file.
        ///
        /// \param[in] file Name of ECLIPSE output file.
        ///
        /// \return Internalised ERT file contents.
        FilePtr loadFile(const boost::filesystem::path& file);

        /// Extract non-neighbouring connections from ECLIPSE model
        ///
        /// \param[in] G ERT Grid representation.
        ///
        /// \param[in] init ERT representation of INIT source.
        ///
        /// \return Model's non-neighbouring connections, including those
        ///         between main and local grids.
        std::vector<ecl_nnc_type>
        loadNNC(const ecl_grid_type* G,
                const ecl_file_type* init);

        /// Facility for mapping a result set vector defined either on all
        /// global cells or on the explicitly active cells (ACTNUM != 0) to
        /// the global cells.
        class ScatterMap
        {
        public:
            /// Constructor
            ///
            /// \param[in] G ERT Grid representation.
            ///
            /// \param[in] pvol Vector of pore-volumes on all global cells
            ///                 of \p G.  Typically obtained through
            ///                 function getPVolVector().
            ScatterMap(const ecl_grid_type*       G,
                       const std::vector<double>& pvol);

            /// Retrive global cell indices of all active cells.
            std::vector<std::size_t> activeGlobalCells() const;

            /// Map input vector to all global cells.
            ///
            /// \param[in] x Input vector, defined on the explicitly active
            ///              cells, all global cells or some other subset
            ///              (e.g., all non-neighbouring connections).
            ///
            /// \return Input vector mapped to global cells or unchanged if
            /// input is defined on some other subset.
            template <typename T>
            std::vector<T>
            scatterToGlobal(const std::vector<T>& x) const;

        private:
            /// Explicit mapping between ACTNUM!=0 cells and global cells.
            struct ID {
                std::size_t act;
                std::size_t glob;
            };

            /// Number of explicitly active cells (SUM(ACTNUM != 0)).
            std::size_t nact_;

            /// Number of global cells.
            std::size_t nglob_;

            /// Active subset of global cells.
            std::vector<ID> subset_;
        };

        /// Facility for cached querying of ECLIPSE result set.
        class GlobalCellData
        {
        public:
            /// Constructor
            ///
            /// \param[in] G ERT Grid representation.
            ///
            /// \param[in] pvol Vector of pore-volumes on all global cells
            ///                 of \p G.  Typically obtained through
            ///                 function getPVolVector().
            GlobalCellData(const ecl_grid_type*       grid,
                           const std::vector<double>& pvol);

            /// Assign result set backing object.
            ///
            /// \param[in] src Name of result set (restart) file, possibly
            /// unified.
            void assignDataSource(const boost::filesystem::path& src);

            /// Assign result set backing object.
            ///
            /// Specialised convenience interface.
            ///
            /// \param[in] src Internalised result set.  Typically
            /// corresponds to the INIT object.
            void assignDataSource(FilePtr src);

            using DataVector = std::vector<double>;

            /// Retrive ERT representation of current result set.
            ///
            /// Convenience method that typically accesses the result set
            /// that was assigned through assignDataSource() with a \c
            /// FilePtr argument.
            const ecl_file_type* getDataSource() const;

            /// Retrive global cell indices of all active cells.
            std::vector<std::size_t> activeGlobalCells() const;

            /// Predicate for whether or not a particular vector exists in
            /// the current result set.
            bool haveVector(const std::string& vector) const;

            /// Retrieve particular vector from current result set.
            ///
            /// \param[in] vector Name of result vector.  Clients should
            ///                   ensure existence of data through
            ///                   haveVector() before calling getVector().
            ///
            /// \param[in] occurrence Selected temporal vector.  Essentially
            ///                       the report step number.
            ///
            /// \return Result vector corresponding to named quantity at
            /// given time.
            const DataVector&
            getVector(const std::string& vector,
                      const int          occurrence = 0);

        private:
            using DataHandle     = std::unique_ptr<DataVector>;
            using TemporalData   = std::map<int, DataHandle>;
            using DataCollection = std::map<std::string, TemporalData>;

            /// Map results from active to global cells.
            ScatterMap map_;

            /// Current result set.
            ECL::FilePtr src_;

            /// Result vector cache.
            DataCollection coll_;

            /// Load a result vector from backing store.
            ///
            /// \param[in] vector Name of result vector.  Clients should
            ///                   ensure existence of data through
            ///                   haveVector() before calling getVector().
            ///
            /// \param[in] occurrence Selected temporal vector.  Essentially
            ///                       the report step number.
            ///
            /// \return Raw result vector as represented in the backing
            /// store.  Must be mapped to global cells prior to access.
            DataVector
            load(const std::string& vector,
                 const int          occurrence);
        };
    } // namespace ECL

    /// Derive neighbourship relations between active cells from Cartesian
    /// relations.
    class CellCollection
    {
    public:
        /// Canonical directions of Cartesian neighbours.
        enum class Direction { I, J, K };

        /// Constructor.
        ///
        /// \param[in] G ERT Grid representation.
        ///
        /// \param[in] active_glob_cells Global indices of model's active
        ///       cells (\code ACTNUM != 0 && pore_volume > 0 \endcode).
        CellCollection(const ecl_grid_type*       G,
                       std::vector<std::size_t>&& active_glob_cells);

        /// Retrieve global indices of model's active cells.
        ///
        /// Convenience method that returns the second constructor argument.
        const std::vector<std::size_t>& activeGlobalCells() const;

        /// Retrieve total number of cells in model, including inactive ones.
        ///
        /// Needed to allocate result vectors on global cells.
        std::size_t numGlobalCells() const;

        /// Retrieve active cell ID of particular global cell.
        ///
        /// \param[in] globalCell Index of particular global cell.
        ///
        /// \return Active cell ID of \p globalCell.  Returns negative one
        /// (\code -1 \endcode) if \code globalCell >= numGlobalCells
        /// \endcode or if the global cell is inactive.
        int activeCell(const std::size_t globalCell) const;

        /// Retrieve active cell ID of particular global cell's neighbour in
        /// given Cartesian direction.
        ///
        /// \param[in] globalCell Index of particular global cell.
        ///
        /// \param[in] d Cartesian direction in which to look for a
        /// neighbouring cell.
        ///
        /// \return Active cell ID of \p globalCell's neighbour in direction
        /// \d.  Returns negative one (\code -1 \endcode) if \code
        /// globalCell >= numGlobalCells \endcode or if the global cell is
        /// inactive, or if there is no neighbour in direction \p d (e.g.,
        /// if purported neighbour would be outside model).
        int cartesianNeighbour(const std::size_t globalCell,
                               const Direction   d) const;

    private:
        using IJKTuple = std::array<std::size_t,3>;

        /// Number of global cells in X (I) direction.
        std::size_t nx_;

        /// Number of global cells in Y (J) direction.
        std::size_t ny_;

        /// Number of global cells in Z (K) direction.
        std::size_t nz_;

        /// Global indices of model's active cells.
        std::vector<std::size_t> glob_cell_;

        /// Active cell ID of model's global cells.  Negative one (\code -1
        /// \endcode) if inactive.
        std::vector<int> active_id_;

        /// Retrieve number of active cells.  Equivalent to \code
        /// glob_cell_.size() \endcode, but with a named expression.
        std::size_t numActive() const;

        /// Compute linear index of global cell from explicit (I,J,K) tuple.
        ///
        /// \param[in] ijk Explicit (I,J,K) tuple of global cell.
        ///
        /// \return Linear index (natural ordering) of global cell (I,J,K).
        std::size_t globIdx(const IJKTuple& ijk) const;

        /// Decompose global (linear) cell index into its (I,J,K) index
        /// tuple.
        ///
        /// \param[in] globalCell Index of particular global cell.  Must be
        ///    in the range \code [0 .. numGlobalCells()) \endcode.
        ///
        /// \return Index triplet of \p globalCell's location within model.
        IJKTuple ind2sub(const std::size_t globalCell) const;
    };
} // Anonymous namespace

// ======================================================================

std::vector<double>
ECL::getPVolVector(const ecl_grid_type* G,
                   const ecl_file_type* init)
{
    auto make_szt = [](const int i)
    {
        return static_cast<std::vector<double>::size_type>(i);
    };

    const auto nx = make_szt(ecl_grid_get_nx(G));
    const auto ny = make_szt(ecl_grid_get_ny(G));
    const auto nz = make_szt(ecl_grid_get_nz(G));

    const auto nglob = nx * ny * nz;

    auto pvol = std::vector<double>(nglob, 1.0);

    if (ecl_file_has_kw(init, "PORV")) {
        auto porv =
            ecl_file_iget_named_kw(init, "PORV", 0);

        assert ((make_szt(ecl_kw_get_size(porv)) == nglob)
                && "Pore-volume must be provided for all global cells");

        ecl_kw_get_data_as_double(porv, pvol.data());
    }

    return pvol;
}

ECL::GridPtr
ECL::loadCase(const boost::filesystem::path& grid)
{
    auto G = GridPtr{
        ecl_grid_load_case(grid.generic_string().c_str())
    };

    if (! G) {
        std::ostringstream os;

        os << "Failed to load ECL Grid from '"
           << grid.generic_string() << '\'';

        throw std::invalid_argument(os.str());
    }

    return G;
}

ECL::FilePtr
ECL::loadFile(const boost::filesystem::path& file)
{
    // Read-only, keep open between requests
    const auto open_flags = 0;

    auto F = FilePtr{
        ecl_file_open(file.generic_string().c_str(), open_flags)
    };

    if (! F) {
        std::ostringstream os;

        os << "Failed to load ECL File object from '"
           << file.generic_string() << '\'';

        throw std::invalid_argument(os.str());
    }

    return F;
}

std::vector<ecl_nnc_type>
ECL::loadNNC(const ecl_grid_type* G,
             const ecl_file_type* init)
{
    auto make_szt = [](const int n)
    {
        return static_cast<std::vector<ecl_nnc_type>::size_type>(n);
    };

    auto nncData = std::vector<ecl_nnc_type>{};

    const auto numNNC = make_szt(ecl_nnc_export_get_size(G));

    if (numNNC > 0) {
        nncData.resize(numNNC);

        ecl_nnc_export(G, init, nncData.data());
    }

    std::sort(nncData.begin(), nncData.end(),
              [](const ecl_nnc_type& nd1, const ecl_nnc_type& nd2) {
                  return nd1.input_index < nd2.input_index;
              });

    return nncData;
}

// ======================================================================

ECL::ScatterMap::ScatterMap(const ecl_grid_type*       G,
                            const std::vector<double>& pvol)
{
    auto make_szt = [](const int i)
    {
        return static_cast<std::size_t>(i);
    };

    this->nact_ = make_szt(ecl_grid_get_nactive(G));

    {
        const auto nx = make_szt(ecl_grid_get_nx(G));
        const auto ny = make_szt(ecl_grid_get_ny(G));
        const auto nz = make_szt(ecl_grid_get_nz(G));

        this->nglob_ = nx * ny * nz;
    }

    subset_.clear();
    subset_.reserve(this->nact_);

    if (pvol.empty()) {
        for (decltype(ecl_grid_get_nactive(G))
                 act = 0, nact = ecl_grid_get_nactive(G);
             act < nact; ++act)
        {
            const auto glob =
                make_szt(ecl_grid_get_global_index1A(G, act));

            this->subset_.push_back(ID{ make_szt(act), glob });
        }
    }
    else {
        for (decltype(ecl_grid_get_nactive(G))
                 act = 0, nact = ecl_grid_get_nactive(G);
             act < nact; ++act)
        {
            const auto glob =
                make_szt(ecl_grid_get_global_index1A(G, act));

            if (pvol[glob] > 0.0) {
                this->subset_.push_back(ID{ make_szt(act), glob });
            }
        }
    }
}

std::vector<std::size_t>
ECL::ScatterMap::activeGlobalCells() const
{
    auto active = std::vector<std::size_t>{};
    active.reserve(this->subset_.size());

    for (const auto& id : this->subset_) {
        active.push_back(id.glob);
    }

    return active;
}

template <typename T>
std::vector<T>
ECL::ScatterMap::scatterToGlobal(const std::vector<T>& x) const
{
    // Assume that input vector 'x' is either defined on explicit
    // notion of active cells (ACTNUM != 0) or on all global cells or
    // some other contiguous index set (e.g., the NNCs).

    if (x.size() != static_cast<decltype(x.size())>(this->nact_)) {
        // Input not defined on explictly active cells.  Let caller
        // deal with it.  This typically corresponds to the set of
        // global cells or the list of NNCs.
        return x;
    }

    auto y = std::vector<T>(this->nglob_);

    for (const auto& i : this->subset_) {
        y[i.glob] = x[i.act];
    }

    return y;
}

// ======================================================================

ECL::GlobalCellData::GlobalCellData(const ecl_grid_type*       grid,
                                    const std::vector<double>& pvol)
    : map_(grid, pvol)
{
}

void
ECL::GlobalCellData::assignDataSource(const boost::filesystem::path& src)
{
    this->src_ = ECL::loadFile(src);
}

void
ECL::GlobalCellData::assignDataSource(ECL::FilePtr src)
{
    this->src_ = std::move(src);
}

const ecl_file_type*
ECL::GlobalCellData::getDataSource() const
{
    return this->src_.get();
}

std::vector<std::size_t>
ECL::GlobalCellData::activeGlobalCells() const
{
    return this->map_.activeGlobalCells();
}

bool
ECL::GlobalCellData::haveVector(const std::string& vector) const
{
    return ecl_file_has_kw(this->src_.get(), vector.c_str());
}

const ECL::GlobalCellData::DataVector&
ECL::GlobalCellData::getVector(const std::string& vector,
                               const int          occurrence)
{
    auto& p = this->coll_[vector][occurrence];

    if (! p) {
        auto x = this->map_.scatterToGlobal(load(vector, occurrence));

        p = DataHandle{ new DataVector{ std::move(x) } };
    }

    return *p;
}

ECL::GlobalCellData::DataVector
ECL::GlobalCellData::load(const std::string& vector,
                          const int          occurrence)
{
    auto x = ecl_file_iget_named_kw(this->src_.get(),
                                    vector.c_str(),
                                    occurrence);

    auto y = DataVector(ecl_kw_get_size(x));

    ecl_kw_get_data_as_double(x, y.data());

    return y;
}

// ======================================================================

CellCollection::CellCollection(const ecl_grid_type*       G,
                               std::vector<std::size_t>&& active_glob_cells)
    : nx_       (static_cast<std::size_t>(ecl_grid_get_nx(G)))
    , ny_       (static_cast<std::size_t>(ecl_grid_get_ny(G)))
    , nz_       (static_cast<std::size_t>(ecl_grid_get_nz(G)))
    , glob_cell_(std::move(active_glob_cells))
    , active_id_(nx_ * ny_ * nz_, -1)
{
    auto act = 0;

    for (const auto& id : glob_cell_) {
        this->active_id_[static_cast<std::size_t>(id)] = act++;
    }
}

const std::vector<std::size_t>&
CellCollection::activeGlobalCells() const
{
    return this->glob_cell_;
}

std::size_t
CellCollection::numActive() const
{
    return this->glob_cell_.size();
}

std::size_t
CellCollection::numGlobalCells() const
{
    return this->active_id_.size();
}

int
CellCollection::activeCell(const std::size_t globalCell) const
{
    if (globalCell >= numGlobalCells()) { return -1; }

    return this->active_id_[globalCell];
}

int
CellCollection::cartesianNeighbour(const std::size_t globalCell,
                                   const Direction   d) const
{
    if (globalCell >= numGlobalCells()) { return -1; }

    auto ijk = ind2sub(globalCell);

    if      (d == Direction::I) { ijk[0] += 1; }
    else if (d == Direction::J) { ijk[1] += 1; }
    else if (d == Direction::K) { ijk[2] += 1; }
    else {
        return -1;
    }

    const auto globNeigh = globIdx(ijk);

    if (globNeigh >= numGlobalCells()) { return -1; }

    return this->active_id_[globNeigh];
}

std::size_t
CellCollection::globIdx(const IJKTuple& ijk) const
{
    if (ijk[0] >= this->nx_) { return -1; }
    if (ijk[1] >= this->ny_) { return -1; }
    if (ijk[2] >= this->nz_) { return -1; }

    return ijk[0] + nx_*(ijk[1] + ny_*ijk[2]);
}

CellCollection::IJKTuple
CellCollection::ind2sub(const std::size_t globalCell) const
{
    assert (globalCell < numGlobalCells());

    auto ijk = IJKTuple{};
    auto g   = globalCell;

    ijk[0] = g % this->nx_;  g /= this->nx_;
    ijk[1] = g % this->ny_;
    ijk[2] = g / this->ny_;  assert (ijk[2] < this->nz_);

    assert (globIdx(ijk) == globalCell);

    return ijk;
}

// =====================================================================

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
