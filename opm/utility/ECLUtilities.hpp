/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil ASA.

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

#ifndef OPM_ECLUTILITIES_HEADER_INCLUDED
#define OPM_ECLUTILITIES_HEADER_INCLUDED

#include <ert/util/ert_unique_ptr.hpp>
#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_grid.h>
#include <ert/ecl/ecl_nnc_export.h>
#include <boost/filesystem.hpp>
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
#include <vector>


namespace Opm {
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
} // namespace Opm


#endif // OPM_ECLUTILITIES_HEADER_INCLUDED
