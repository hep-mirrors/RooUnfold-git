/*===========================================================================*/
/*
 * BEGIN ROOUNFOLD COPYRIGHT
 * RooUnfold — Unfolding library for particle-physics inverse problems
 *
 * Copyright © 2007–2025 CERN and the authors’ respective research institutions
 * Authors (by git history of this file):
 *   - Tim Adye (2007–2011)
 *   - Vincent Croft (2021)
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

//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      2D test of RooUnfold package using toy MC generated according to PDFs
//      defined in RooUnfoldTestPdf.cxx or RooUnfoldTestPdfRooFit.cxx.
//      This is the main program. The actual tests are performed using the
//      RooUnfoldTestHarness2D class.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#include "RooUnfoldTestHarness2D.h"

RooUnfoldTestHarness2D *test2d = 0;

//==============================================================================
// Routine to run with parameters specified as a string
//==============================================================================

void RooUnfoldTest2D(const char *args = "")
{
   // If run interactively, remove canvas and all histograms that might have been
   // created with a previous invocation.
   delete test2d;
   test2d = 0;
   gDirectory->Clear();

   test2d = new RooUnfoldTestHarness2D("RooUnfoldTest2D", args);
   test2d->Run();
}

#ifndef __CINT__

//==============================================================================
// Main program when run stand-alone
//==============================================================================

int main(int argc, char **argv)
{
   RooUnfoldTestHarness2D maintest2d("RooUnfoldTest2D", argc, argv);
   return maintest2d.Run();
}

#endif
