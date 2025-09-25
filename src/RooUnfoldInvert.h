/*===========================================================================*/
/*
 * BEGIN ROOUNFOLD COPYRIGHT
 * RooUnfold — Unfolding library for particle-physics inverse problems
 *
 * Copyright © 2010–2025 CERN and the authors’ respective research institutions
 * Authors (by git history of this file):
 *   - Richard Claridge (2010)
 *   - Tim Adye (2010–2011, 2022)
 *   - Carsten Burgard (2019, 2022, 2024)
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
//! \class RooUnfoldInvertT
//! \brief Unfolding class using inversion of the response matrix. This does not produce
//!      good results and is designed to illustrate the need for more sophisticated
//!      unfolding techniques
//! \author Richard Claridge <richard.claridge@stfc.ac.uk> & Tim Adye <T.J.Adye@rl.ac.uk>
//==============================================================================

#ifndef ROOUNFOLDINVERT_H_
#define ROOUNFOLDINVERT_H_

#include "RooUnfold.h"
#include "RooUnfoldResponse.h"

class TDecompSVD;

template <class Hist, class Hist2D>
class RooUnfoldInvertT : public RooUnfoldT<Hist, Hist2D> {

public:
   RooUnfoldInvertT();                                                     // default constructor
   RooUnfoldInvertT(const char *name, const char *title);                  // named constructor
   RooUnfoldInvertT(const TString &name, const TString &title);            // named constructor
   RooUnfoldInvertT(const RooUnfoldInvertT<Hist, Hist2D> &rhs);            // copy constructor
   virtual ~RooUnfoldInvertT();                                            // destructor
   RooUnfoldInvertT &operator=(const RooUnfoldInvertT<Hist, Hist2D> &rhs); // assignment operator
   RooUnfoldInvertT(const RooUnfoldResponseT<Hist, Hist2D> *res, const Hist *meas, const char *name = 0,
                    const char *title = 0);

   virtual RooUnfolding::Algorithm GetAlgorithm() const override;
   virtual void Reset() override;
   TDecompSVD *Impl();
   const TMatrixD &InverseResponse() const;

protected:
   virtual void Unfold() const override;
   virtual void GetCov() const override;

private:
   void Init();
   Bool_t InvertResponse() const;

protected:
   // instance variables
   mutable TDecompSVD *_svd;
   mutable TMatrixD *_resinv;

public:
   ClassDefOverride(RooUnfoldInvertT, 1) // Unregularised unfolding
};

//! \class RooUnfoldInvert
//! \brief specialization of RooUnfoldInvertT for TH1/TH2 objects
typedef RooUnfoldInvertT<TH1, TH2> RooUnfoldInvert;
#ifndef NOROOFIT
//! \class RooFitUnfoldInvert
//! \brief specialization of RooUnfoldInvertT for RooAbsReal objects
typedef RooUnfoldInvertT<RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> RooFitUnfoldInvert;
#endif

#endif /*ROOUNFOLDINVERT_H_*/
