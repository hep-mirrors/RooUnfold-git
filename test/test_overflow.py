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

if __name__ == "__main__":
    parms = {"overflow": ["0", "1", "2"], "verbose": ["3"]}
    ref_file_name = "../ref/overflow.ref"
    test_name = "overflow"
    field_to_compare = ["unfoldoverflow"]
    perform_test(parms, ref_file_name, test_name, field_to_compare)
