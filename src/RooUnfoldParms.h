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
//      Optimisation of regularisation parameter class
//
// Author: Richard Claridge <richard.claridge@stfc.ac.uk>
//
//==============================================================================

#ifndef ROOUNFOLDPARMS_H_
#define ROOUNFOLDPARMS_H_

#include "TNamed.h"
#include "RooUnfold.h"

class TH1;
class TProfile;

class RooUnfoldParms : public TNamed {
public:
   RooUnfoldParms(const RooUnfold *unfold_in = 0, RooUnfolding::ErrorTreatment err = RooUnfolding::kCovariance,
                  const TH1 *truth = 0);
   virtual ~RooUnfoldParms();
   TProfile *GetChi2();
   TProfile *GetRMSError();
   TProfile *GetMeanResiduals();
   TH1 *GetRMSResiduals();
   const RooUnfold *unfold;              // Input object from RooUnfold
   RooUnfolding::ErrorTreatment doerror; // Set error calculation method
   const TH1 *hTrue;                     // Truth Distribution
   void SetMinParm(double min);
   void SetMaxParm(double max);
   void SetStepSizeParm(double size);

private:
   bool _done_math;
   TH1 *hrms;      // Output plot
   TProfile *hch2; // Output plot
   TProfile *herr; // Output plot
   TProfile *hres; // Output plot
   void DoMath();
   void Init();
   Double_t _maxparm;      // Maximum parameter
   Double_t _minparm;      // Minimum parameter
   Double_t _stepsizeparm; // Step size
public:
   ClassDef(RooUnfoldParms, 0) // Optimisation of unfolding regularisation parameter
};
#endif /*ROOUNFOLDPARMS_H_*/
