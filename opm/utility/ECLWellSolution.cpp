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
#include <opm/utility/ECLUtilities.hpp>
#include <ert/ecl/ecl_kw_magic.h>

namespace Opm
{

    ECLWellSolution::ECLWellSolution(const Path& restart)
        : restart_path_(restart)
    {
    }





    std::vector<ECLWellSolution::WellData>
    ECLWellSolution::solution(const int occurrence)
    {
        auto restart = ECL::loadFile(restart_path_);
        auto intehead = loadIntField(restart.get(), "INTEHEAD", occurrence);
        return {};
    }




    std::vector<double>
    ECLWellSolution::loadDoubleField(ecl_file_type* restart,
                                     const std::string& fieldname,
                                     const int occurrence)
    {
        std::vector<double> field_data;
        ecl_file_push_block(restart);
        {
            // Select the block containing the requested occurrence.
            ecl_file_subselect_block(restart, SEQNUM_KW, occurrence);
            {
                const int local_occurrence = 0; // TODO: with LGRs this might need reconsideration.
                ecl_kw_type* keyword = ecl_file_iget_named_kw(restart, fieldname.c_str(), local_occurrence);
                field_data.resize(ecl_kw_get_size(keyword));
                ecl_kw_get_data_as_double(keyword, field_data.data());
            }
        }
        ecl_file_pop_block(restart);
        return field_data;
    }




    std::vector<int>
    ECLWellSolution::loadIntField(ecl_file_type* restart,
                                  const std::string& fieldname,
                                  const int occurrence)
    {
        std::vector<int> field_data;
        ecl_file_push_block(restart);
        {
            // Select the block containing the requested occurrence.
            ecl_file_subselect_block(restart, SEQNUM_KW, occurrence);
            {
                const int local_occurrence = 0; // TODO: with LGRs this might need reconsideration.
                ecl_kw_type* keyword = ecl_file_iget_named_kw(restart, fieldname.c_str(), local_occurrence);
                field_data.resize(ecl_kw_get_size(keyword));
                ecl_kw_get_memcpy_int_data(keyword, field_data.data());
            }
        }
        ecl_file_pop_block(restart);
        return field_data;
    }

} // namespace Opm
