# BEGIN ROOUNFOLD COPYRIGHT
# RooUnfold — Unfolding library for particle-physics inverse problems
#
# Copyright © 2021–2025 CERN and the authors’ respective research institutions
# Authors (by git history of this file):
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

from test_utils import perform_test
import os

if __name__ == "__main__":
    parms = {"method": ["1", "2", "3", "4", "5", "6"], "verbose": ["3"]}
    ref_file_name = "../ref/test_methods.ref"
    test_name = "test_methods"
    field_to_compare = ["unfold"]
    perform_test(parms, ref_file_name, test_name, field_to_compare)

    ## For plotting
    command_str = "../build/RooUnfoldTest ploterrors=2"
    os.system(command_str)
    command_str = "../build/RooUnfoldTest ploterrors=1"
    os.system(command_str)
    command_str = "../build/RooUnfoldTest plotparms=2"
    os.system(command_str)
