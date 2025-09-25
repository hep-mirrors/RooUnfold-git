/*===========================================================================*/
/*
 * BEGIN ROOUNFOLD COPYRIGHT
 * RooUnfold — Unfolding library for particle-physics inverse problems
 *
 * Copyright © 2021–2022 CERN and the authors’ respective research institutions
 * Authors (by git history of this file):
 *   - Vincent Croft (2021)
 *   - Tim Adye (2022)
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

#ifndef __UNITTESTS_H__
#define __UNITTESTS_H__

#include "RooUnfoldResponse.h"
#include "TVector.h"
#include <string>

void RooUnfoldGenerate();
void RooUnfoldGenerateVariable();
RooUnfoldResponse BuildRooUnfoldResponse(std::string);
TVector BuildRooUnfoldBayes(int);
void WriteRooUnfoldBayes(int);
int TestBayes(int);
const char *test_bayes();
#endif
