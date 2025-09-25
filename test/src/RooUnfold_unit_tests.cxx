/*===========================================================================*/
/*
 * BEGIN ROOUNFOLD COPYRIGHT
 * RooUnfold — Unfolding library for particle-physics inverse problems
 *
 * Copyright © 2021–2025 CERN and the authors’ respective research institutions
 * Authors (by git history of this file):
 *   - Vincent Croft (2021)
 *   - Tim Adye (2022)
 *   - Roel Aaij (2025)
 *
 * Note: Authorship is inferred from Git history. Copyright is held by CERN and by the
 * respective research institutions employing the authors at the time of contribution.
 *
 * License: BSD-3-Clause
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * This file header was generated automatically from repository history.
 * END ROOUNFOLD COPYRIGHT
 */
/*===========================================================================*/

#include <iostream>
#include <fstream>
#include <string>
#include "dbg.h"
#include "minunit.h"
using std::cout;
using std::endl;

#include "unittests.h"

int TestBayes(int mode = 0)
{
   auto u = BuildRooUnfoldBayes(mode);
   int n = u.GetNrows();

   std::vector<float> r;
   std::string ref_name;
   if (mode == 0)
      ref_name = "../ref/bayes.ref";
   else if (mode == 1)
      ref_name = "../ref/bayes_var.ref";
   std::ifstream ref(ref_name);
   if (ref.is_open()) {
      std::string line;
      while (getline(ref, line)) {
         r.push_back(std::stof(line));
      }
      ref.close();
   }
   check(n == Int_t(r.size()), "error in bayes: %d bins is not %zu", n, r.size());

   for (int i = 0; i < n; i++) {
      check(EQ(u[i], r[i], 1), "error in bayes bin: %d as %f is not %f", i, u[i], r[i]);
   }

   return 1;
error:
   return -1;
}

const char *run_bayes()
{
   mu_assert(TestBayes(0) == 1, "binwise test of Bayes failed");
   mu_assert(TestBayes(1) == 1, "binwise test of variable binning in Bayes failed");
   return NULL;
}

const char *test_bayes()
{
   RooUnfoldGenerate();
   RooUnfoldGenerateVariable();
   mu_suite_start();
   mu_run_test(run_bayes);
   return NULL;
}
