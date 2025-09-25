/*===========================================================================*/
/*
 * BEGIN ROOUNFOLD COPYRIGHT
 * RooUnfold — Unfolding library for particle-physics inverse problems
 *
 * Copyright © 2010–2025 CERN and the authors’ respective research institutions
 * Authors (by git history of this file):
 *   - Richard Claridge (2010)
 *   - Tim Adye (2010–2011)
 *   - Carsten Burgard (2019)
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
//      Graph Drawing Class for use with RooUnfold.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Richard Claridge <richard.claridge@stfc.ac.uk>
//
//==============================================================================

#ifndef ROOUNFOLDERRORS_H_
#define ROOUNFOLDERRORS_H_

#include "TNamed.h"

class TH1;
class TNtuple;

#include "RooUnfold.h"

class RooUnfoldErrors : public TNamed {

public:
   int toys;          // Number of toys
   RooUnfold *unfold; // Input unfolding object
   const TH1 *hTrue;
   RooUnfoldErrors(int NToys, RooUnfold *unfold, const TH1 *Truth = 0);
   virtual ~RooUnfoldErrors();
   TNtuple *Chi2();

   TH1 *RMSResiduals();
   TH1 *UnfoldingError();

private:
   void CreatePlots();
   void CreatePlotsWithChi2();
   TH1 *h_err;             // Output plot
   TH1 *h_err_res;         // Output plot
   TNtuple *hchi2;         // Output plot
   void GraphParameters(); //
   double xlo;             // Minimum x-axis value
   double xhi;             // Maximum x-axis value
   int ntx;                // Number of bins in true distribution

public:
   ClassDef(RooUnfoldErrors, 0) // Show unfolding errors
};

#endif
