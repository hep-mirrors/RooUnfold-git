/*===========================================================================*/
/*
 * BEGIN ROOUNFOLD COPYRIGHT
 * RooUnfold — Unfolding library for particle-physics inverse problems
 *
 * Copyright © 2007–2025 CERN and the authors’ respective research institutions
 * Authors (by git history of this file):
 *   - Tim Adye (2007, 2010–2011, 2013, 2017–2018)
 *   - Richard Claridge (2010)
 *   - Carsten Burgard (2019)
 *   - Pim Verschuuren (2019–2020)
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

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class RooUnfold - ;
#pragma link C++ class RooUnfoldBayes + ;
#pragma link C++ class RooUnfoldSvd + ;
#pragma link C++ class RooUnfoldSvd::SVDUnfold + ;
#pragma link C++ class RooUnfoldBinByBin + ;
#pragma link C++ class RooUnfoldIds + ;
#pragma link C++ class RooUnfoldGP + ;
#pragma link C++ class RooUnfoldPoisson + ;
#pragma link C++ class TUnfoldV17 + ;
#pragma link C++ class TUnfoldSysV17 + ;
#pragma link C++ class TUnfoldBinningV17 + ;
#pragma link C++ class TUnfoldDensityV17 + ;
#pragma link C++ class TUnfoldBinningXMLV17 + ;
#pragma link C++ class TUnfoldIterativeEMV17 + ;
#ifndef NOTUNFOLD
#pragma link C++ class RooUnfoldTUnfold + ;
#endif
#pragma link C++ class RooUnfoldResponseT < TH1, TH2> - ;
#pragma link C++ class RooUnfoldResponse + ;
#ifndef NOROOFIT
#pragma link C++ class RooUnfolding::RooFitHist + ;
#pragma link C++ class RooUnfoldResponseT < RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> - ;
#pragma link C++ class RooUnfoldT < RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> - ;
#pragma link C++ class RooUnfoldFunc + ;
#pragma link C++ class RooUnfoldSpec + ;
#pragma link C++ class RooUnfoldInvertT < RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> + ;
#pragma link C++ class RooUnfoldGPT < RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> + ;
#pragma link C++ class RooUnfoldPoissonT < RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> + ;
#pragma link C++ class RooUnfoldBayesT < RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> + ;
#pragma link C++ class RooUnfoldSvdT < RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> + ;
#pragma link C++ class RooUnfoldSvdT < RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> ::SVDUnfold + ;
#pragma link C++ class RooUnfoldBinByBinT < RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> + ;
#pragma link C++ class RooUnfoldIdsT < RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> + ;
#pragma link C++ class RooUnfoldTUnfoldT < RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> + ;
#pragma link C++ class RooFitUnfoldResponse + ;
#endif
#pragma link C++ class RooUnfoldErrors + ;
#pragma link C++ class RooUnfoldParms + ;
#pragma link C++ class RooUnfoldInvert + ;

#endif
