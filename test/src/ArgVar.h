/*===========================================================================*/
/*
 * BEGIN ROOUNFOLD COPYRIGHT
 * RooUnfold — Unfolding library for particle-physics inverse problems
 *
 * Copyright © 2010–2025 CERN and the authors’ respective research institutions
 * Authors (by git history of this file):
 *   - Tim Adye (2010–2011)
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
//      Parse argument list for parameter settings
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef ARGVAR_H
#define ARGVAR_H

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TNamed.h"
#include "TString.h"
#endif

class ArgVar : public TNamed {
public:
   const char *defhelp; //!
   Int_t *ivar;         //!
   Int_t idef;
   Double_t *fvar; //!
   Double_t fdef;
   TString *svar; //!
   TString sdef;
   bool setdef;
   ArgVar() : defhelp(0), ivar(0), idef(0), fvar(0), fdef(0), svar(0), setdef(false) {}
   ArgVar(const char *n, Int_t *v)
      : TNamed(n, ""), defhelp(0), ivar(v), idef(0), fvar(0), fdef(0), svar(0), setdef(false)
   {
   }
   ArgVar(const char *n, Double_t *v)
      : TNamed(n, ""), defhelp(0), ivar(0), idef(0), fvar(v), fdef(0), svar(0), setdef(false)
   {
   }
   ArgVar(const char *n, TString *v)
      : TNamed(n, ""), defhelp(0), ivar(0), idef(0), fvar(0), fdef(0), svar(v), setdef(false)
   {
   }
   ArgVar(const char *n, Int_t *v, Int_t d, const char *h = 0, const char *dh = 0)
      : TNamed(n, h), defhelp(dh), ivar(v), idef(d), fvar(0), fdef(0), svar(0), setdef(true)
   {
   }
   ArgVar(const char *n, Double_t *v, Double_t d, const char *h = 0, const char *dh = 0)
      : TNamed(n, h), defhelp(dh), ivar(0), idef(0), fvar(v), fdef(d), svar(0), setdef(true)
   {
   }
   ArgVar(const char *n, TString *v, const TString &d, const char *h = 0, const char *dh = 0)
      : TNamed(n, h), defhelp(dh), ivar(0), idef(0), fvar(0), fdef(0), svar(v), sdef(d), setdef(true)
   {
   }
};

#endif
