/*===========================================================================*/
/*
 * BEGIN ROOUNFOLD COPYRIGHT
 * RooUnfold — Unfolding library for particle-physics inverse problems
 *
 * Copyright © 2010–2025 CERN and the authors’ respective research institutions
 * Authors (by git history of this file):
 *   - Tim Adye (2010–2011, 2018)
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

#ifndef ARGVARS_H
#define ARGVARS_H

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include "TObject.h"
#include "TNamed.h"
#include "TList.h"
#include "TString.h"
#endif

#include "ArgVar.h"

class ArgVars : public TObject {
private:
   TList lst;
   static bool CmpOpt(const char *p, const char *opt, const char *s);
   ArgVar *Find(const char *name) const;
   ArgVars &Add(ArgVar *arg);

public:
   ArgVars() { lst.SetOwner(); }
   virtual ~ArgVars();
   ArgVars &Add(const char *name, Int_t *var) { return Add(new ArgVar(name, var)); }
   ArgVars &Add(const char *name, Double_t *var) { return Add(new ArgVar(name, var)); }
   ArgVars &Add(const char *name, Int_t *var, Int_t def, const char *help = 0, const char *defhelp = 0)
   {
      return Add(new ArgVar(name, var, def, help, defhelp));
   }
   ArgVars &Add(const char *name, Double_t *var, Double_t def, const char *help = 0, const char *defhelp = 0)
   {
      return Add(new ArgVar(name, var, def, help, defhelp));
   }
   ArgVars &Add(const char *name, TString *var, const TString &def, const char *help = 0, const char *defhelp = 0)
   {
      return Add(new ArgVar(name, var, def, help, defhelp));
   }
   ArgVars &Add(const ArgVars &args);
   ArgVars &SetDefault(const char *name, Int_t def);
   ArgVars &SetDefault(const char *name, Double_t def);
   ArgVars &SetDefault(const char *name, const TString &def);
   Int_t SetArgs(int argc, const char *const *argv, bool split = false) const;
   void SetDefaults() const;
   virtual void Print(std::ostream &o, const char *sep = " ") const;
   virtual void Print(const char *sep = " ") const { Print(std::cout, sep); }
   void Usage(const char *prog) const;
   void ArgHelp(std::ostream &o) const;
};

#ifndef NOINLINE
#include "ArgVars.cxx"
#endif

#endif
