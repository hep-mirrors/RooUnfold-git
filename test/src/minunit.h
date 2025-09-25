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

#undef NDEBUG
#ifndef _minunit_h
#define _minunit_h

#include <stdio.h>
#include "dbg.h"
#include <stdlib.h>

#define mu_suite_start() const char *message = NULL

#define mu_assert(test, message) \
   if (!(test)) {                \
      log_err(message);          \
      return message;            \
   }
#define mu_run_test(test)         \
   debug("\n-----%s", " " #test); \
   message = test();              \
   if (message)                   \
      return message;

#define RUN_TESTS(name)                       \
   int main(int, char *argv[])                \
   {                                          \
      debug("----- RUNNING: %s", argv[0]);    \
      printf("----\nRUNNING: %s\n", argv[0]); \
      const char *result = name();            \
      if (result != 0) {                      \
         printf("FAILED: %s\n", result);      \
      } else {                                \
         printf("ALL TESTS PASSED\n");        \
      }                                       \
      exit(result != 0);                      \
   }

#endif
