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

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <opm/utility/ECLWellSolution.hpp>
#include <ert/ecl/ecl_kw_magic.h>
#include <stdexcept>
#include <sstream>

namespace Opm
{

    namespace {

        struct EclFileSelectReportBlock
        {
            EclFileSelectReportBlock(ecl_file_type* file, const int number)
                : file_(file)
            {
                ecl_file_push_block(file_);
                ecl_file_subselect_block(file_, SEQNUM_KW, number);
            }
            ~EclFileSelectReportBlock()
            {
                ecl_file_pop_block(file_);
            }
            ecl_file_type* file_;
        };

        ERT::ert_unique_ptr<ecl_file_type, ecl_file_close>
        load(const boost::filesystem::path& filename)
        {
            // Read-only, keep open between requests
            const auto open_flags = 0;
            using FilePtr = ERT::ert_unique_ptr<ecl_file_type, ecl_file_close>;
            FilePtr file(ecl_file_open(filename.generic_string().c_str(), open_flags));
            if (!file) {
                std::ostringstream os;
                os << "Failed to load ECL File object from '"
                   << filename.generic_string() << '\'';
                throw std::invalid_argument(os.str());
            }
            return file;
        }

    }




    ECLWellSolution::ECLWellSolution(const Path& restart_filename)
        : restart_(load(restart_filename))
    {
    }





    std::vector<ECLWellSolution::WellData>
    ECLWellSolution::solution(const int occurrence)
    {
        EclFileSelectReportBlock(restart_.get(), occurrence);
        {
            auto intehead = loadIntField(INTEHEAD_KW);
            const int unit = intehead[INTEHEAD_UNIT_INDEX];
            return {};
        }
    }




    std::vector<double>
    ECLWellSolution::loadDoubleField(const std::string& fieldname)
    {
        std::vector<double> field_data;
        const int local_occurrence = 0; // TODO: with LGRs this might need reconsideration.
        ecl_kw_type* keyword = ecl_file_iget_named_kw(restart_.get(), fieldname.c_str(), local_occurrence);
        field_data.resize(ecl_kw_get_size(keyword));
        ecl_kw_get_data_as_double(keyword, field_data.data());
        return field_data;
    }




    std::vector<int>
    ECLWellSolution::loadIntField(const std::string& fieldname)
    {
        std::vector<int> field_data;
        const int local_occurrence = 0; // TODO: with LGRs this might need reconsideration.
        ecl_kw_type* keyword = ecl_file_iget_named_kw(restart_.get(), fieldname.c_str(), local_occurrence);
        field_data.resize(ecl_kw_get_size(keyword));
        ecl_kw_get_memcpy_int_data(keyword, field_data.data());
        return field_data;
    }

} // namespace Opm
