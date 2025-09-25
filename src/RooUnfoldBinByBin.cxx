/*===========================================================================*/
/*
 * BEGIN ROOUNFOLD COPYRIGHT
 * RooUnfold — Unfolding library for particle-physics inverse problems
 *
 * Copyright © 2007–2025 CERN and the authors’ respective research institutions
 * Authors (by git history of this file):
 *   - Tim Adye (2007, 2009–2011, 2022)
 *   - Richard Claridge (2010)
 *   - Carsten Burgard (2019, 2022)
 *   - Pim Verschuuren (2019–2020)
 *   - David Hutchcroft (2021)
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

/*! \class RooUnfoldBinByBinT
Uses the correction factor method to unfold the distribution by looking at each bin individually.
This method cannot account for bin migration and as such cannot unfold reliably if a bias/smearing effects are applied.
Can only handle 1 dimensional distributions
True and measured distributions must have the same binning
*/

#include "RooUnfoldBinByBin.h"
#include "RooUnfoldHelpers.h"

#include <iostream>
#include "TH1.h"
#include "TH2.h"

#include "RooUnfoldResponse.h"

template <class Hist, class Hist2D>
RooUnfoldBinByBinT<Hist, Hist2D>::RooUnfoldBinByBinT(const RooUnfoldBinByBinT<Hist, Hist2D> &rhs)
   : RooUnfoldT<Hist, Hist2D>(rhs)
{
   // Copy constructor.
   this->GetSettings();
}

template <class Hist, class Hist2D>
RooUnfolding::Algorithm RooUnfoldBinByBinT<Hist, Hist2D>::GetAlgorithm() const
{
   //! return the unfolding algorithm used
   return RooUnfolding::kBinByBin;
}

template <class Hist, class Hist2D>
RooUnfoldBinByBinT<Hist, Hist2D>::RooUnfoldBinByBinT(const RooUnfoldResponseT<Hist, Hist2D> *res, const Hist *meas,
                                                     const char *name, const char *title)
   : RooUnfoldT<Hist, Hist2D>(res, meas, name, title)
{
   // Constructor with response matrix object and measured unfolding input histogram.
   this->GetSettings();
}

template <class Hist, class Hist2D>
RooUnfoldBinByBinT<Hist, Hist2D>::~RooUnfoldBinByBinT()
{
}

template <class Hist, class Hist2D>
TVectorD *RooUnfoldBinByBinT<Hist, Hist2D>::Impl()
{
   return &_specialcache._factors;
}

#include <unistd.h>

template <class Hist, class Hist2D>
void RooUnfoldBinByBinT<Hist, Hist2D>::Unfold() const
{

   const TVectorD &vmeas(this->Vmeasured());
   const TVectorD &vtrain(this->_res->Vmeasured());
   const TVectorD &vtruth(this->_res->Vtruth());
   const TVectorD &fakes(this->_res->Vfakes());

   int nm = vtrain.GetNrows();
   int nt = vtruth.GetNrows();

   Double_t fac = 0.0;
   if (this->_res->HasFakes()) {
      fac = vtrain.Sum();
      if (fac != 0.0)
         fac = vmeas.Sum() / fac;
      if (this->_verbose >= 1)
         std::cout << "Subtract " << fac * fakes.Sum() << " fakes from measured distribution" << std::endl;
   }

   this->_cache._rec.ResizeTo(nt);
   this->_specialcache._factors.ResizeTo(nt);
   Int_t nb = std::min(nm, nt);
   for (int i = 0; i < nb; i++) {
      Double_t train = vtrain[i] - fakes[i];
      if (train == 0.0)
         continue;
      Double_t c = vtruth[i] / train;
      this->_specialcache._factors[i] = c;
      this->_cache._rec[i] = c * (vmeas[i] - fac * fakes[i]);
   }

   this->_cache._unfolded = true;
}

template <class Hist, class Hist2D>
void RooUnfoldBinByBinT<Hist, Hist2D>::GetCov() const
{
   //! Get covariance matrix
   int nt = this->response()->Vtruth().GetNrows();
   int nm = this->response()->Vmeasured().GetNrows();
   const TMatrixD &covmeas(this->GetMeasuredCov());
   this->_cache._cov.ResizeTo(nt, nt);
   Int_t nb = std::min(nm, nt);
   for (int i = 0; i < nb; i++)
      for (int j = 0; j < nb; j++)
         this->_cache._cov(i, j) = pow(this->_specialcache._factors[i], 2) * covmeas(i, j);
   this->_cache._haveCov = true;
}

template <class Hist, class Hist2D>
RooUnfoldBinByBinT<Hist, Hist2D>::RooUnfoldBinByBinT() : RooUnfoldT<Hist, Hist2D>()
{
   // Default constructor. Use Setup() to prepare for unfolding.
   this->GetSettings();
}

template <class Hist, class Hist2D>
RooUnfoldBinByBinT<Hist, Hist2D>::RooUnfoldBinByBinT(const char *name, const char *title)
   : RooUnfoldT<Hist, Hist2D>(name, title)
{
   // Basic named constructor. Use Setup() to prepare for unfolding.
   this->GetSettings();
}

template <class Hist, class Hist2D>
RooUnfoldBinByBinT<Hist, Hist2D>::RooUnfoldBinByBinT(const TString &name, const TString &title)
   : RooUnfoldT<Hist, Hist2D>(name, title)
{
   // Basic named constructor. Use Setup() to prepare for unfolding.
   this->GetSettings();
}

template <class Hist, class Hist2D>
RooUnfoldBinByBinT<Hist, Hist2D> &
RooUnfoldBinByBinT<Hist, Hist2D>::operator=(const RooUnfoldBinByBinT<Hist, Hist2D> &rhs)
{
   // Assignment operator for copying RooUnfoldBinByBinT settings.
   this->Assign(rhs);
   return *this;
}

template class RooUnfoldBinByBinT<TH1, TH2>;
ClassImp(RooUnfoldBinByBin)

#ifndef NOROOFIT
   typedef RooUnfoldBinByBinT<RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> RooFitUnfoldBinByBin;
template class RooUnfoldBinByBinT<RooUnfolding::RooFitHist, RooUnfolding::RooFitHist>;
ClassImp(RooFitUnfoldBinByBin)
#endif
