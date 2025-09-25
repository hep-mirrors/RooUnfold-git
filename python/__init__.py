# BEGIN ROOUNFOLD COPYRIGHT
# RooUnfold — Unfolding library for particle-physics inverse problems
#
# Copyright © 2021–2025 CERN and the authors’ respective research institutions
# Authors (by git history of this file):
#   - Carsten Burgard (2021)
#   - Roel Aaij (2025)
#
# Note: Authorship is inferred from Git history. Copyright is held by CERN and by the
# respective research institutions employing the authors at the time of contribution.
#
# License: BSD-3-Clause
# SPDX-License-Identifier: BSD-3-Clause
#
# This file header was generated automatically from repository history.
# END ROOUNFOLD COPYRIGHT

import ROOT

ROOT.gSystem.Load("libRooUnfold.so")
RooUnfoldResponse = ROOT.RooUnfoldResponse
RooUnfoldBayes = ROOT.RooUnfoldBayes
RooUnfoldBinByBin = ROOT.RooUnfoldBayes
RooUnfoldInvert = ROOT.RooUnfoldInvert
RooUnfoldSvd = ROOT.RooUnfoldSvd
RooUnfoldTUnfold = ROOT.RooUnfoldTUnfold
RooUnfoldSpec = ROOT.RooUnfoldSpec
RooUnfoldFunc = ROOT.RooUnfoldFunc
RooFitUnfoldResponse = ROOT.RooFitUnfoldResponse
RooFitUnfoldBayes = ROOT.RooFitUnfoldBayes
RooFitUnfoldSvd = ROOT.RooFitUnfoldSvd
RooFitUnfoldBinByBin = ROOT.RooFitUnfoldBinByBin
RooFitHist = ROOT.RooUnfolding.RooFitHist
RooUnfolding = ROOT.RooUnfolding
