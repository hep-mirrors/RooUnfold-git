# BEGIN ROOUNFOLD COPYRIGHT
# RooUnfold — Unfolding library for particle-physics inverse problems
#
# Copyright © 2021–2025 CERN and the authors’ respective research institutions
# Authors (by git history of this file):
#   - Archit Agrawal (2021)
#   - Vincent Croft (2021)
#   - Roel Aaij (2025)
#
# Note: Authorship is inferred from Git history. Copyright is held by CERN and by the
# respective research institutions employing the authors at the time of contribution.
#
# License: BSD-3-Clause
# SPDX-License-Identifier: BSD-3-Clause
#
# This file header was generated automatically from repository history.
# END ROOUNFOLD COPYRIGHT

from test_utils import perform_test, get_combination

if __name__ == "__main__":
    parms = {
        "method": ["1", "2", "4"],
        "ftrainx": ["7"],
        "doerror": ["2"],
        "dosys": ["1"],
        "verbose": ["3"],
        "seed": ["42"],
    }

    parms_name = list(parms.keys())
    parms_name.sort()
    combined_parm = get_combination(parms, parms_name)

    parms = {
        "method": ["2"],
        "ftrainx": ["7"],
        "doerror": ["3"],
        "dosys": ["1"],
        "ntoys": ["50", "500"],
        "verbose": ["3"],
        "seed": ["42"],
    }

    parms_name = list(parms.keys())
    parms_name.sort()
    combined_parm.extend(get_combination(parms, parms_name))

    ref_file_name = "../ref/test_uncertainty.ref"
    test_name = "test_uncertainty"
    field_to_compare = ["uncertainty"]
    perform_test(combined_parm, ref_file_name, test_name, field_to_compare, is_combined=True)
