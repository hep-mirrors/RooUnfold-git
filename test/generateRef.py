# BEGIN ROOUNFOLD COPYRIGHT
# RooUnfold — Unfolding library for particle-physics inverse problems
#
# Copyright © 2021–2025 CERN and the authors’ respective research institutions
# Authors (by git history of this file):
#   - Vincent Croft (2021)
#   - Archit Agrawal (2022)
#   - Carsten Burgard (2022)
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

import os
import json
from test_utils import get_combination, delete_files
from test_utils import get_field, write_field

basefolder = ".."


def get_parms_name(parms):
    parms_name = list(parms.keys())
    parms_name.sort()
    return parms_name


def genrate_ref(combined_parm, field_to_compare, ref_file_name):
    all_output = {}
    delete_files()
    for single_parm in combined_parm:
        command_str = basefolder + "/build/RooUnfoldTest " + single_parm
        if os.WEXITSTATUS(os.system(command_str)) != 0:
            print(
                "[ERROR] Cannot run RooUnfoldTest, make sure you compiled with '-D RooUnfoldTests=ON' and are executing this script from 'build'"
            )
            exit(1)
        u = get_field("RooUnfoldTest.root", field_to_compare)
        all_output[single_parm] = u
        delete_files()

    write_field(all_output, ref_file_name)


def geneate_ref_methods():
    parms = {"method": ["1", "2", "3", "4", "5", "6"], "verbose": ["3"]}
    ref_file_name = basefolder + "/ref/test_methods.ref"
    field_to_compare = ["unfold"]
    combined_parm = get_combination(parms, get_parms_name(parms))
    genrate_ref(combined_parm, field_to_compare, ref_file_name)


def generate_ref_uncertainity():
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
        "seed": ["42"],
        "verbose": ["3"],
    }

    parms_name = list(parms.keys())
    parms_name.sort()
    combined_parm.extend(get_combination(parms, parms_name))

    ref_file_name = basefolder + "/ref/test_uncertainty.ref"
    test_name = "test_uncertainty"
    field_to_compare = ["uncertainty"]
    genrate_ref(combined_parm, field_to_compare, ref_file_name)


def generate_ref_fakes():
    parms = {"method": ["1", "2", "3", "4", "5", "6"], "addfakes": ["1"], "verbose": ["3"]}

    ref_file_name = basefolder + "/ref/test_fakes.ref"
    field_to_compare = ["unfold"]
    combined_parm = get_combination(parms, get_parms_name(parms))
    genrate_ref(combined_parm, field_to_compare, ref_file_name)


def generate_ref_bin_correlation():
    parms = {"bincorr": ["1"], "seed": ["42"], "verbose": ["3"]}
    ref_file_name = basefolder + "/ref/test_correlation.ref"
    test_name = "test_correlation"
    field_to_compare = ["uncertainty"]
    combined_parm = get_combination(parms, get_parms_name(parms))
    genrate_ref(combined_parm, field_to_compare, ref_file_name)


def generate_ref_overflow():
    parms = {"overflow": ["0", "1", "2"], "verbose": ["3"]}
    ref_file_name = basefolder + "/ref/overflow.ref"
    field_to_compare = ["unfoldoverflow"]
    combined_parm = get_combination(parms, get_parms_name(parms))
    genrate_ref(combined_parm, field_to_compare, ref_file_name)


def generate_ref_2D():
    ref_file_name = basefolder + "/ref/test_2D.ref"
    field_to_compare = ["unfold2D"]
    all_output = {}
    command_str = basefolder + "/build/RooUnfoldTest2D  verbose=3"
    os.system(command_str)
    u = get_field("RooUnfoldTest2D.root", field_to_compare)
    all_output["default"] = u
    os.system("rm -f RooUnfoldTest2D.root")
    os.system("rm -f RooUnfoldTest2D.ps")
    write_field(all_output, ref_file_name)


def generate_ref_3D():
    ref_file_name = basefolder + "/ref/test_3D.ref"
    field_to_compare = ["unfold3D"]
    all_output = {}
    command_str = basefolder + "/build/RooUnfoldTest3D  verbose=3"
    os.system(command_str)
    u = get_field("RooUnfoldTest3D.root", field_to_compare)
    all_output["default"] = u
    os.system("rm -f RooUnfoldTest3D.root")
    os.system("rm -f RooUnfoldTest3D.ps")
    write_field(all_output, ref_file_name)


if __name__ == "__main__":
    geneate_ref_methods()
    generate_ref_uncertainity()
    generate_ref_fakes()
    generate_ref_bin_correlation()
    generate_ref_overflow()
    generate_ref_2D()
    generate_ref_3D()
