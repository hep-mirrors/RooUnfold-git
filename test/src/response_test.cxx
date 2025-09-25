/*===========================================================================*/
/*
 * BEGIN ROOUNFOLD COPYRIGHT
 * RooUnfold — Unfolding library for particle-physics inverse problems
 *
 * Copyright © 2021–2025 CERN and the authors’ respective research institutions
 * Authors (by git history of this file):
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
//      Unit tests for generating the RooUnfoldResponse
//
// Authors: Vincent Croft <vincent.croft@cern.ch>
//
//==============================================================================

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "RooUnfoldResponse.h"
#include "unittests.h"
#include <string>

RooUnfoldResponse BuildRooUnfoldResponse(std::string filename = "response.root")
{
   TFile *f = new TFile(filename.c_str(), "OPEN");
   TH2D *h_response = (TH2D *)f->Get("res");
   TH1D *h_gen = (TH1D *)f->Get("gen");
   TH1D *h_sim = (TH1D *)f->Get("sim");
   RooUnfoldResponse response(h_sim, h_gen, h_response);
   f->Close();
   return response;
}