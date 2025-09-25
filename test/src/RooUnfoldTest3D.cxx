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
//      3D test of RooUnfold package using toy MC generated according to PDFs
//      defined in RooUnfoldTestPdf.cxx or RooUnfoldTestPdfRooFit.cxx.
//      This is the main program. The actual tests are performed using the
//      RooUnfoldTestHarness3D class.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#include "RooUnfoldTestHarness3D.h"

RooUnfoldTestHarness3D *test3d = 0;

//==============================================================================
// Routine to run with parameters specified as a string
//==============================================================================

void RooUnfoldTest3D(const char *args = "")
{
   // If run interactively, remove canvas and all histograms that might have been
   // created with a previous invocation.
   delete test3d;
   test3d = 0;
   gDirectory->Clear();

   test3d = new RooUnfoldTestHarness3D("RooUnfoldTest3D", args);
   test3d->Run();
}

#ifndef __CINT__

//==============================================================================
// Main program when run stand-alone
//==============================================================================

int main(int argc, char **argv)
{
   RooUnfoldTestHarness3D maintest3d("RooUnfoldTest3D", argc, argv);
   return maintest3d.Run();
}

#endif
