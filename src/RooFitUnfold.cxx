/*===========================================================================*/
/*
 * BEGIN ROOUNFOLD COPYRIGHT
 * RooUnfold — Unfolding library for particle-physics inverse problems
 *
 * Copyright © 2019–2025 CERN and the authors’ respective research institutions
 * Authors (by git history of this file):
 *   - Carsten Burgard (2019, 2021, 2024)
 *   - Lydia Brenner (2019)
 *   - Pim Verschuuren (2019–2021)
 *   - Tim Adye (2022, 2024)
 *   - Callum McCracken (2024)
 *   - Mars Lyukova (2024)
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

#ifndef NOROOFIT
#include <sstream>
#include "RooUnfold.h"
#include "RooFitUnfold.h"
#include "RooUnfoldHelpers.h"
#include "RooUnfoldTH1Helpers.h"
#include "RooUnfoldFitHelpers.h"
#include "RooHistFunc.h"
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 12, 0)
#include "RooRealSumFunc.h"
#else
#pragma message("RooFitUnfold requires ROOT version 6.12 or later (missing RooRealSumFunc.h)")
#endif
#include "RooPrintable.h"
#include "RooHistPdf.h"
#include "RooExtendPdf.h"
#include "RooProdPdf.h"
#include "RooGaussian.h"
#include "RooPoisson.h"
#include "RooStats/HistFactory/PiecewiseInterpolation.h"
#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include "RooStats/HistFactory/ParamHistFunc.h"
#include "RooProduct.h"
#ifndef NO_WRAPPERPDF
#include "RooWrapperPdf.h"
#endif
#include "RooBinning.h"

#include "THStack.h"
#include <RooFit/Detail/JSONInterface.h>

using namespace RooUnfolding;

RooUnfoldFunc::RooUnfoldFunc(const char *name, const char *title,
                             const RooUnfoldT<RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> *unf, bool useDensity)
   : RooAbsReal(name, title), _useDensity(useDensity)
{
   //! constructor
   auto unfolding = RooUnfoldT<RooUnfolding::RooFitHist, RooUnfolding::RooFitHist>::New(
      unf->GetAlgorithm(), unf->response(), unf->Hmeasured(), unf->GetRegParm(), unf->GetName(), unf->GetTitle());
   if (unf->Hbkg())
      unfolding->SetBkg(unf->Hbkg());
   unfolding->SetTruth(unf->Htruth());
   unfolding->SetMeasuredCov(unf->GetMeasuredCov());
   this->_unfolding = dynamic_cast<RooUnfoldT<RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> *>(unfolding);
   this->_unfolding->SetVerbose(0);
   const RooUnfoldResponseT<RooFitHist, RooFitHist> *res = this->_unfolding->response();
   if (res) {
      const RooFitHist *htruth = res->Htruth();
      if (htruth) {
         this->addServer(*(htruth->func()));
         for (int i = 0; i < dim(htruth); ++i) {
            this->addServer(*htruth->obs(i));
         }
      }
      const RooFitHist *hfakes = res->Hfakes();
      if (hfakes) {
         this->addServer(*(hfakes->func()));
         for (int i = 0; i < dim(hfakes); ++i) {
            this->addServer(*hfakes->obs(i));
         }
      }
      const RooFitHist *hresponse = res->Hresponse();
      if (hresponse) {
         this->addServer(*(hresponse->func()));
         for (int i = 0; i < dim(hresponse); ++i) {
            this->addServer(*hresponse->obs(i));
         }
      }
      const RooFitHist *hmeasured = res->Hmeasured();
      if (hmeasured) {
         this->addServer(*(hmeasured->func()));
         for (int i = 0; i < dim(hmeasured); ++i) {
            this->addServer(*hmeasured->obs(i));
         }
      }
   }

   const RooFitHist *hmeasured = this->_unfolding->Hmeasured();

   if (hmeasured) {
      this->addServer(*(hmeasured->func()));
      for (int i = 0; i < dim(hmeasured); ++i) {
         this->addServer(*hmeasured->obs(i));
      }
   }
}
RooUnfoldFunc::RooUnfoldFunc() : _unfolding(NULL)
{
   //! constructor
}
RooUnfoldFunc::~RooUnfoldFunc()
{
   //! destructor
   delete _unfolding;
}

Bool_t RooUnfoldFunc::redirectServersHook(const RooAbsCollection &newServerList, Bool_t mustReplaceAll,
                                          Bool_t nameChange, Bool_t isRecursive)
{
   //! redirect servers
   RooUnfoldResponseT<RooFitHist, RooFitHist> *res = this->_unfolding->response();
   if (res) {
      RooFitHist *htruth = res->Htruth();
      if (htruth) {
         htruth->replace(newServerList);
      }
      RooFitHist *hfakes = res->Hfakes();
      if (hfakes) {
         hfakes->replace(newServerList);
      }
      RooFitHist *hresponse = res->Hresponse();
      if (hresponse) {
         hresponse->replace(newServerList);
      }
      RooFitHist *hmeasured = res->Hmeasured();
      if (hmeasured) {
         hmeasured->replace(newServerList);
      }
   }
   RooFitHist *hmeasured = this->_unfolding->Hmeasured();
   if (hmeasured) {
      hmeasured->replace(newServerList);
   }
   return RooAbsReal::redirectServersHook(newServerList, mustReplaceAll, nameChange, isRecursive);
}

const RooUnfoldT<RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> *RooUnfoldFunc::unfolding() const
{
   //! retrieve the unfolding object
   return this->_unfolding;
}

std::list<Double_t> *RooUnfoldFunc::binBoundaries(RooAbsRealLValue &obs, Double_t xlo, Double_t xhi) const
{
   //! retrieve the list of bin boundaries
   return this->_unfolding->response()->Htruth()->func()->binBoundaries(obs, xlo, xhi);
}

std::list<Double_t> *RooUnfoldFunc::plotSamplingHint(RooAbsRealLValue &obs, Double_t xlo, Double_t xhi) const
{
   //! retrieve the sampling hint
   return this->_unfolding->response()->Htruth()->func()->plotSamplingHint(obs, xlo, xhi);
}

Double_t RooUnfoldFunc::getValV(const RooArgSet *set) const
{
   //! return the value
   this->_curNormSet = set;
   return RooAbsReal::getValV(set);
}

RooArgList *RooUnfoldFunc::makeParameterList() const
{
   //! return a list of all parameters in this function
   RooArgSet obs;
   for (size_t d = 0; d < this->unfolding()->response()->Hresponse()->dim(); ++d) {
      obs.add(*(this->unfolding()->response()->Hresponse()->obs(d)));
   }
   RooArgSet *pset = this->getParameters(&obs);

   RooArgList *list = new RooArgList(*pset);
   delete pset;
   return list;
}

bool RooUnfoldFunc::isDensity() const
{
   //! return true if the return value is density-corrected, false otherwise
   return this->_useDensity;
}
void RooUnfoldFunc::setDensity(bool d)
{
   //! set if the return value should be density-corrected
   this->_useDensity = d;
}

Double_t RooUnfoldFunc::evaluate() const
{
   //! call getVal on the internal function
   std::map<std::string, double> snapshot;
   this->_unfolding->response()->Hresponse()->saveSnapshot(snapshot);
   int bin = this->_unfolding->response()->Htruth()->bin();
   this->_unfolding->ForceRecalculation();
   this->_unfolding->response()->Htruth()->checkValidity();
   double v = this->_unfolding->Vunfold()[bin];
   if (this->_useDensity) {
      v /= binVolume(this->_unfolding->response()->Htruth(), bin, false);
   }
   this->_unfolding->response()->Hresponse()->loadSnapshot(snapshot);
   return v;
}

Bool_t RooUnfoldFunc::isBinnedDistribution(const RooArgSet &obs) const
{
   //! check if this PDF is a binned distribution in the given observable
   return this->_unfolding->response()->Hresponse()->func()->isBinnedDistribution(obs);
}

Bool_t RooUnfoldFunc::checkObservables(const RooArgSet *nset) const
{
   //! call checkOvservables on the response
   return this->_unfolding->response()->Hresponse()->func()->checkObservables(nset);
}

Bool_t RooUnfoldFunc::forceAnalyticalInt(const RooAbsArg &arg) const
{
   //! force the analytical integral
   return this->_unfolding->response()->Htruth()->func()->forceAnalyticalInt(arg);
}

Int_t RooUnfoldFunc::getAnalyticalIntegralWN(RooArgSet &allVars, RooArgSet &numVars, const RooArgSet *normSet,
                                             const char *rangeName) const
{
   //! retrieve the analytical integral status
   return this->_unfolding->response()->Htruth()->func()->getAnalyticalIntegralWN(allVars, numVars, normSet, rangeName);
}

Double_t
RooUnfoldFunc::analyticalIntegralWN(Int_t /*code*/, const RooArgSet * /*normSet*/, const char * /*rangeName*/) const
{
   //! retrieve the analytical integral status
   double val = 0;
   auto vec = this->_unfolding->Vunfold();
   for (int i = 0; i < vec.GetNrows(); ++i) {
      // assuming that density correction has been applied already
      val += vec[i];
   }
   return val;
}

void RooUnfoldFunc::printMetaArgs(std::ostream &os) const
{
   //! printing helper function
   return this->_unfolding->response()->Htruth()->func()->printMetaArgs(os);
}

RooAbsArg::CacheMode RooUnfoldFunc::canNodeBeCached() const
{
   return this->_unfolding->response()->Htruth()->func()->canNodeBeCached();
}

void RooUnfoldFunc::setCacheAndTrackHints(RooArgSet &arg)
{
   this->_unfolding->response()->Htruth()->func()->setCacheAndTrackHints(arg);
}
TObject *RooUnfoldFunc::clone(const char *newname) const
{
   //! produce a clone (deep copy) of this object
   return new RooUnfoldFunc(newname ? newname : this->GetName(), this->GetTitle(), this->_unfolding, this->_useDensity);
}

namespace {
bool readToken(TString &instr, std::vector<TString> &tokens)
{
   int pos = instr.First(",");
   if (pos == -1) {
      tokens.push_back(instr);
      return false;
   } else {
      tokens.push_back(instr(0, pos));
      instr.Remove(0, pos + 1);
      return true;
   }
}
} // namespace

#include <RooAbsCollection.h>

#if ROOT_VERSION_CODE < ROOT_VERSION(6, 18, 0)
namespace {

struct IteratorHelper {
   RooFIter itr;
   RooAbsArg *nextItem;
   IteratorHelper(const RooAbsCollection &c);
   IteratorHelper();
   RooAbsArg *operator++();
   bool operator!=(const IteratorHelper &other);
   bool operator!=(const RooAbsArg *other);
   RooAbsArg *operator*();
};

IteratorHelper::IteratorHelper(const RooAbsCollection &c) : itr(c.fwdIterator()), nextItem(itr.next()) {}
IteratorHelper::IteratorHelper() : itr(), nextItem(NULL) {}
RooAbsArg *IteratorHelper::operator++()
{
   nextItem = itr.next();
   return nextItem;
}
bool IteratorHelper::operator!=(const IteratorHelper &other)
{
   return this->nextItem != other.nextItem;
}
bool IteratorHelper::operator!=(const RooAbsArg *other)
{
   return this->nextItem != other;
}

RooAbsArg *IteratorHelper::operator*()
{
   return nextItem;
}
} // namespace

::IteratorHelper begin(const RooAbsCollection &c)
{
   return ::IteratorHelper(c);
}

::IteratorHelper end(const RooAbsCollection &)
{
   return ::IteratorHelper();
}
#endif

RooUnfoldSpec::RooUnfoldSpec(const char *name, const char *title, const TH1 *truth, const char *obs_truth,
                             const TH1 *reco, const char *obs_reco, const TH2 *response, const TH1 *data,
                             bool includeUnderflowOverflow, double errorThreshold, bool useDensity)
   : RooUnfoldSpec(name, title, truth, obs_truth, reco, obs_reco, response, 0, data, includeUnderflowOverflow,
                   errorThreshold, useDensity)
{
   //! constructor forwarding
}

namespace {
template <class Hist>
RooRealVar *makeObs(const char *name, RooUnfolding::Variable<Hist> v)
{
   RooRealVar *var = new RooRealVar(name, name, v._min, v._min, v._max);
   var->setConstant(true);
   if (v.irregular()) {
      var->setBinning(RooBinning(v._nBins, &(v._bounds[0])));
   } else {
      var->setBins(v._nBins);
   }
   return var;
}
} // namespace

RooUnfoldSpec::RooUnfoldSpec(const char *name, const char *title, const TH1 *truth, const char *obs_truth,
                             const TH1 *reco, const char *obs_reco, const TH2 *response, const TH1 *bkg,
                             const TH1 *data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity)
   : TNamed(name, title)
{
   //! constructor
   int d = dim(truth);
   if (d != dim(reco)) {
      throw std::runtime_error("inconsistent dimensionality between truth and reco histograms!");
   }

   TString obs_truth_s(obs_truth);
   TString obs_reco_s(obs_reco);
   std::vector<TString> obs_truth_v;
   std::vector<TString> obs_reco_v;
   bool more_truth = false;
   bool more_reco = false;
   for (int i = 0; i < d; ++i) {
      more_truth = ::readToken(obs_truth_s, obs_truth_v);
      more_reco = ::readToken(obs_reco_s, obs_reco_v);
      if (!more_truth || !more_reco)
         break;
   }
   if (more_truth)
      throw std::runtime_error(
         TString::Format("encountered additional characters on truth observable list: '%s'", obs_truth_s.Data())
            .Data());
   if (more_reco)
      throw std::runtime_error(
         TString::Format("encountered additional characters on reco observable list: '%s'", obs_reco_s.Data()).Data());
   if ((int)obs_truth_v.size() != d)
      throw std::runtime_error(
         TString::Format("truth observable list is too short for %d dimensions: '%s'", d, obs_truth).Data());
   if ((int)obs_reco_v.size() != d)
      throw std::runtime_error(
         TString::Format("reco observable list is too short for %d dimensions: '%s'", d, obs_truth).Data());

   RooArgList truth_vars;
   for (int i = 0; i < d; ++i) {
      auto v = var(truth, (RooUnfolding::Dimension)i);
      if (i > 0 && v._nBins == 1)
         continue;
      RooRealVar *obs = makeObs(obs_truth_v[i], v);
      truth_vars.add(*obs);
   }

   RooArgList reco_vars;
   for (int i = 0; i < d; ++i) {
      auto v = var(reco, (RooUnfolding::Dimension)i);
      if (i > 0 && v._nBins == 1)
         continue;
      RooRealVar *obs = makeObs(obs_reco_v[i], v);
      reco_vars.add(*obs);
   }
   this->setup(truth, truth_vars, reco, reco_vars, response, bkg, data, includeUnderflowOverflow, errorThreshold,
               useDensity);
}

RooUnfoldSpec::~RooUnfoldSpec()
{
   //! destructor
}

#ifndef NOROOFIT

RooUnfoldSpec::RooUnfoldSpec(const char *name, const char *title, const TH1 *truth_th1, const RooArgList &obs_truth,
                             const TH1 *reco_th1, const RooArgList &obs_reco, const TH2 *response_th1, const TH1 *bkg,
                             const TH1 *data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity)
   : TNamed(name, title)
{
   //! constructor
   this->setup(truth_th1, obs_truth, reco_th1, obs_reco, response_th1, bkg, data, includeUnderflowOverflow,
               errorThreshold, useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char *name, const char *title, const TH1 *truth_th1, const RooArgList &obs_truth,
                             const TH1 *reco_th1, const RooArgList &obs_reco, const TH2 *response_th1, RooAbsReal *bkg,
                             RooDataHist *data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity)
   : TNamed(name, title)
{
   //! constructor
   this->_bkg.setNominal(bkg, obs_reco);
   this->_data.setNominal(data, obs_reco);
   this->setup(truth_th1, obs_truth, reco_th1, obs_reco, response_th1, NULL, NULL, includeUnderflowOverflow,
               errorThreshold, useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char *name, const char *title, const TH1 *truth_th1, const RooArgList &obs_truth,
                             RooAbsReal *reco, const RooArgList &obs_reco, const TH2 *response_th1, RooAbsReal *bkg,
                             RooDataHist *data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity)
   : TNamed(name, title)
{
   //! constructor
   this->_reco.setNominal(reco, obs_reco);
   this->_bkg.setNominal(bkg, obs_reco);
   this->_data.setNominal(data, obs_reco);
   this->setup(truth_th1, obs_truth, NULL, obs_reco, response_th1, NULL, NULL, includeUnderflowOverflow, errorThreshold,
               useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char *name, const char *title, const TH1 *truth_th1, RooAbsArg *obs_truth,
                             RooAbsReal *reco, RooAbsArg *obs_reco, const TH2 *response_th1, RooAbsReal *bkg,
                             RooDataHist *data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity)
   : RooUnfoldSpec(name, title, truth_th1, RooArgList(*obs_truth), reco, RooArgList(*obs_reco), response_th1, bkg, data,
                   includeUnderflowOverflow, errorThreshold, useDensity)
{
   //! constructor
}

namespace {
RooAbsReal *makeParamHistFunc(const char *name, const char *title, const RooArgList &observables,
                              const RooArgList &parameters, bool multiplyDensity)
{
   ParamHistFunc *phf = new ParamHistFunc(name, title, observables, parameters);
   if (!multiplyDensity)
      return phf;
   RooArgList inverseBinWidths;
   RooRealVar *obs = (RooRealVar *)(observables.first());
   for (int i = 0; i < obs->getBinning().numBins(); ++i) {
      double width = obs->getBinWidth(i);
      RooRealVar *bw =
         new RooRealVar(TString::Format("%s_bin_%d_widthCorr", obs->GetName(), i),
                        TString::Format("density correction of %s in bin %d", obs->GetTitle(), i), 1. / width);
      bw->setConstant(true);
      inverseBinWidths.add(*bw);
   }
   ParamHistFunc *bwcorr =
      new ParamHistFunc(TString::Format("%s_widthCorr", obs->GetName()),
                        TString::Format("density corrections for %s", obs->GetTitle()), observables, inverseBinWidths);
   RooArgList elems;
   elems.add(*phf);
   elems.add(*bwcorr);
   RooProduct *prod = new RooProduct(TString::Format("%s_x_%s", phf->GetName(), bwcorr->GetName()), title, elems);
   return prod;
}
} // namespace

RooUnfoldSpec::RooUnfoldSpec(const char *name, const char *title, const TH1 *truth_th1, RooAbsArg *obs_truth,
                             const RooArgList &reco_bins, RooAbsArg *obs_reco, const TH2 *response_th1,
                             const RooArgList &bkg_bins, RooDataHist *data, bool includeUnderflowOverflow,
                             double errorThreshold, bool useDensity)
   : TNamed(name, title)
{
   //! constructor
   RooArgList obs_reco_list(*obs_reco);
   RooArgList obs_truth_list(*obs_truth);
   this->_reco.setNominal(
      ::makeParamHistFunc(TString::Format("signal_reco_%s_differential", obs_reco->GetName()).Data(),
                          obs_reco->GetTitle(), obs_reco_list, reco_bins, useDensity),
      obs_reco_list);
   this->_bkg.setNominal(::makeParamHistFunc(TString::Format("bkg_reco_%s_differential", obs_reco->GetName()).Data(),
                                             obs_reco->GetTitle(), obs_reco_list, bkg_bins, useDensity),
                         obs_reco_list);
   this->_data.setNominal(data, obs_reco_list);
   this->setup(truth_th1, obs_truth_list, NULL, obs_reco_list, response_th1, NULL, NULL, includeUnderflowOverflow,
               errorThreshold, useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char *name, const char *title, const TH1 *truth_th1, RooAbsArg *obs_truth,
                             const TH1 *reco_th1, RooAbsArg *obs_reco, const TH2 *response_th1, RooAbsReal *bkg,
                             RooDataHist *data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity)
   : RooUnfoldSpec(name, title, truth_th1, RooArgList(*obs_truth), reco_th1, RooArgList(*obs_reco), response_th1, bkg,
                   data, includeUnderflowOverflow, errorThreshold, useDensity)
{
}

RooUnfoldSpec::RooUnfoldSpec(const char *name, const char *title, const TH1 *truth_th1, RooAbsArg *obs_truth,
                             const TH1 *reco_th1, RooAbsArg *obs_reco, const TH2 *response_th1,
                             const RooArgList &bkg_bins, RooDataHist *data, bool includeUnderflowOverflow,
                             double errorThreshold, bool useDensity)
   : TNamed(name, title)
{
   //! constructor
   RooArgList obs_reco_list(*obs_reco);
   RooArgList obs_truth_list(*obs_truth);
   this->_bkg.setNominal(::makeParamHistFunc(TString::Format("bkg_reco_%s_differential", obs_reco->GetName()),
                                             obs_reco->GetTitle(), obs_reco_list, bkg_bins, useDensity),
                         obs_reco_list);
   this->_data.setNominal(data, obs_reco_list);
   this->setup(truth_th1, obs_truth_list, reco_th1, obs_reco_list, response_th1, NULL, NULL, includeUnderflowOverflow,
               errorThreshold, useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char *name, const char *title, const TH1 *truth_th1, RooAbsArg *obs_truth,
                             const TH1 *reco_th1, RooAbsArg *obs_reco, const TH2 *response_th1, RooAbsReal *measured,
                             bool includeUnderflowOverflow, double errorThreshold, bool useDensity)
   : TNamed(name, title)
{
   //! constructor
   RooArgList obs_reco_list(*obs_reco);
   RooArgList obs_truth_list(*obs_truth);
   this->_data.setNominal(measured, obs_reco_list);
   this->setup(truth_th1, obs_truth_list, reco_th1, obs_reco_list, response_th1, NULL, NULL, includeUnderflowOverflow,
               errorThreshold, useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char *name, const char *title, const TH1 *truth_th1, RooAbsArg *obs_truth,
                             const TH1 *reco_th1, RooAbsArg *obs_reco, const TH2 *response_th1,
                             const RooArgList &measured_bins, bool includeUnderflowOverflow, double errorThreshold,
                             bool useDensity)
   : TNamed(name, title)
{
   //! constructor
   RooArgList obs_reco_list(*obs_reco);
   RooArgList obs_truth_list(*obs_truth);
   this->_data.setNominal(::makeParamHistFunc(TString::Format("measured_reco_%s_differential", obs_reco->GetName()),
                                              obs_reco->GetTitle(), obs_reco_list, measured_bins, useDensity),
                          obs_reco_list);
   this->setup(truth_th1, obs_truth_list, reco_th1, obs_reco_list, response_th1, NULL, NULL, includeUnderflowOverflow,
               errorThreshold, useDensity);
}

namespace {
RooRealSumFunc *
makeRooRealSumFunc(const char *name, const char *title, const RooArgSet &contributions, double densityCorrection)
{
   RooRealVar *binWidth =
      new RooRealVar(TString::Format("%s_binWidthCorrection", name), "bin width correction", densityCorrection);
   binWidth->setConstant(true);
   RooArgList functions;
   RooArgList coefs;
   for (RooAbsArg *obj : contributions) {
      functions.add(*obj);
      coefs.add(*binWidth);
   }
   return new RooRealSumFunc(name, title, functions, coefs);
}
} // namespace

RooUnfoldSpec::RooUnfoldSpec(const char *name, const char *title, const TH1 *truth_th1, RooAbsArg *obs_truth,
                             RooAbsReal *reco, RooAbsArg *obs_reco, const TH2 *response_th1,
                             const RooArgSet &bkg_contributions, RooDataHist *data, bool includeUnderflowOverflow,
                             double errorThreshold, bool useDensity)
   : TNamed(name, title)
{
   //! constructor
   RooArgList obs_reco_list(*obs_reco);
   RooArgList obs_truth_list(*obs_truth);
   this->_reco.setNominal(reco, obs_reco_list);
   this->_data.setNominal(data, obs_reco_list);
   double densityCorr = 1;
   if (useDensity && obs_reco->InheritsFrom(RooRealVar::Class())) {
      densityCorr = 1. / ((RooRealVar *)(obs_reco))->getBinning().averageBinWidth();
   }
   this->_bkg.setNominal(::makeRooRealSumFunc(TString::Format("bkg_reco_%s_differential", obs_reco->GetName()),
                                              obs_reco->GetTitle(), bkg_contributions, densityCorr),
                         obs_reco_list);
   this->setup(truth_th1, obs_truth_list, NULL, obs_reco_list, response_th1, NULL, NULL, includeUnderflowOverflow,
               errorThreshold, useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char *name, const char *title, RooAbsReal *truth, RooAbsArg *obs_truth,
                             RooAbsReal *reco, RooAbsArg *obs_reco, const TH2 *response_th1,
                             const RooArgSet &bkg_contributions, RooDataHist *data, bool includeUnderflowOverflow,
                             double errorThreshold, bool useDensity)
   : TNamed(name, title)
{
   //! constructor
   RooArgList obs_reco_list(*obs_reco);
   RooArgList obs_truth_list(*obs_truth);
   this->_truth.setNominal(truth, obs_truth_list);
   this->_reco.setNominal(reco, obs_reco_list);
   this->_data.setNominal(data, obs_reco_list);
   double densityCorr = 1;
   if (useDensity && obs_reco->InheritsFrom(RooRealVar::Class())) {
      densityCorr = 1. / ((RooRealVar *)(obs_reco))->getBinning().averageBinWidth();
   }
   this->_bkg.setNominal(::makeRooRealSumFunc(TString::Format("bkg_reco_%s_differential", obs_reco->GetName()),
                                              obs_reco->GetTitle(), bkg_contributions, densityCorr),
                         obs_reco_list);
   this->setup(NULL, obs_truth_list, NULL, obs_reco_list, response_th1, NULL, NULL, includeUnderflowOverflow,
               errorThreshold, useDensity);
}

void RooUnfoldSpec::setup(const TH1 *truth_th1, const RooArgList &obs_truth, const TH1 *reco_th1,
                          const RooArgList &obs_reco, const TH2 *response_th1, const TH1 *bkg_th1, const TH1 *data_th1,
                          bool includeUnderflowOverflow, double errorThreshold, bool useDensity)
{
   //! setup helper function
   this->_includeUnderflowOverflow = includeUnderflowOverflow;
   this->_useDensity = useDensity;
   this->_errorThreshold = errorThreshold;
   if (truth_th1)
      this->_truth.setNominal(truth_th1, obs_truth, errorThreshold, includeUnderflowOverflow, !this->_useDensity);
   if (reco_th1)
      this->_reco.setNominal(reco_th1, obs_reco, errorThreshold, includeUnderflowOverflow, !this->_useDensity);
   this->_obs_reco.add(obs_reco);
   this->_obs_all.add(obs_reco);
   this->_obs_truth.add(obs_truth);
   this->_obs_all.add(obs_truth);
   if (response_th1)
      this->_res.setNominal(response_th1, this->_obs_all, errorThreshold, includeUnderflowOverflow, !this->_useDensity);
   if (bkg_th1)
      this->_bkg.setNominal(bkg_th1, obs_reco, errorThreshold, includeUnderflowOverflow, !this->_useDensity);
   if (data_th1)
      this->_data.setNominal(data_th1, obs_reco, 0., includeUnderflowOverflow, !this->_useDensity);
}

RooUnfolding::RooFitHist *RooUnfoldSpec::makeHistogram(const TH1 *hist)
{
   // convert a TH1 into a RooFitHist
   return new RooUnfolding::RooFitHist(hist, this->_obs_truth, this->_includeUnderflowOverflow, this->_errorThreshold,
                                       this->_useDensity);
}

#endif

namespace {
int countBins(const RooArgList &set)
{
   int bins = 1;
   for (int i = 0; i < set.getSize(); ++i) {
      RooRealVar *obs = dynamic_cast<RooRealVar *>(set.at(i));
      if (obs) {
         bins *= obs->numBins();
      } else {
         throw std::runtime_error("Observable is not of type RooRealVar");
      }
   }
   return bins;
}
} // namespace

TMatrixD RooUnfoldSpec::makeCovarianceMatrix() const
{
   const int bins = countBins(_obs_reco);

   // Initialize the covariance matrix with zeros
   TMatrixD covarianceMatrix(bins, bins);

   // Process shape and norm uncertainties for _reco and _bkg
   addToCovarianceMatrix(_reco, covarianceMatrix);
   addToCovarianceMatrix(_bkg, covarianceMatrix);

   return covarianceMatrix;
}

TH2 *RooUnfoldSpec::makeCovarianceHistogram() const
{
   // convert a matrix into a 2d histogram
   auto matrix = makeCovarianceMatrix();
   RooUnfolding::Variable<TH2> obs = RooUnfolding::Variable<TH2>(static_cast<RooRealVar *>(&_obs_reco[0]));
   TH2 *hist = createHist<TH2>(matrix, "covariances", "covariances", obs, obs);
   return hist;
}

namespace {
TH1 *stathist(std::vector<TH1 *> &hists, const std::vector<Variable<TH1>> &vars,
              const RooUnfoldSpec::HistContainer &cont, const char *name, const char *title, int color, bool useDensity,
              bool includeUnderflowOverflow)
{
   // here, we assume poisson statistics, so the variance is equal to the mean
   const auto vec = h2v(cont._nom, false, useDensity);
   TH1 *variance = createHist<TH1>(vec, name, title, vars, includeUnderflowOverflow);
   variance->SetFillColor(color);
   variance->SetDirectory(0);
   variance->GetYaxis()->SetTitle("Variance");
   hists.push_back(variance);
   return variance;
};

void addshapes(std::vector<TH1 *> &hists, const std::vector<Variable<TH1>> &vars,
               const RooUnfoldSpec::HistContainer &cont, const char *name, const char *title, int color,
               bool useDensity, bool includeUnderflowOverflow)
{
   if (!cont._shapes.empty()) {
      int isys = 1;
      for (const auto &var : cont._shapes) {
         auto nominal = h2v(cont._nom, false, useDensity);
         TVectorD variances(nominal.GetNrows());
         auto shapeUp = h2v(var.second[0], false, useDensity);
         auto shapeDown = h2v(var.second[1], false, useDensity);
         for (int i = 0; i < shapeUp.GetNrows(); ++i) {
            variances[i] = std::pow(shapeUp[i] - nominal[i], 2) + std::pow(shapeDown[i] - nominal[i], 2);
         }
         TH1 *varhist = createHist<TH1>(variances, (std::string(name) + " " + var.first).c_str(),
                                        (std::string(title) + " " + var.first).c_str(), vars, includeUnderflowOverflow);
         varhist->SetLineColor(0);
         varhist->SetFillColor(color + isys);
         varhist->SetDirectory(0);
         hists.push_back(varhist);
         ++isys;
      }
   }
};

void addnorms(std::vector<TH1 *> &hists, const std::vector<Variable<TH1>> &vars,
              const RooUnfoldSpec::HistContainer &cont, const char *name, const char *title, int color, bool useDensity,
              bool includeUnderflowOverflow)
{
   if (!cont._norms.empty()) {
      int isys = 1;
      for (const auto &var : cont._norms) {
         // Calculate sum of squares of norm uncertainties
         auto nominal = h2v(cont._nom, false, useDensity);
         for (int i = 0; i < nominal.GetNrows(); ++i) {
            nominal[i] = pow(nominal[i] * (var.second[0] - 1), 2) + pow(nominal[i] * (var.second[1] - 1), 2);
         }
         TH1 *varhist = createHist<TH1>(nominal, (std::string(name) + " " + var.first).c_str(),
                                        (std::string(title) + " " + var.first).c_str(), vars, includeUnderflowOverflow);
         varhist->SetLineColor(0);
         varhist->SetFillColor(color - isys);
         varhist->SetDirectory(0);
         hists.push_back(varhist);
         ++isys;
      }
   }
};
} // namespace

THStack *RooUnfoldSpec::makeMeasuredBreakdownHistogram() const
{
   std::vector<TH1 *> hists;
   std::vector<Variable<TH1>> vars;
   for (auto &obs : _obs_reco) {
      vars.push_back(Variable<TH1>(static_cast<RooRealVar *>(obs)));
   }

   auto *sig_nom = stathist(hists, vars, _reco, "sig_stat", "signal statistics", kRed, this->_useDensity,
                            this->_includeUnderflowOverflow);
   auto *bkg_nom = stathist(hists, vars, _bkg, "bkg_stat", "background statistics", kBlue, this->_useDensity,
                            this->_includeUnderflowOverflow);

   addshapes(hists, vars, _reco, "sig_shape", "signal shape systematic", kRed, this->_useDensity,
             this->_includeUnderflowOverflow);
   addshapes(hists, vars, _bkg, "bkg_shape", "background shape systematic", kBlue, this->_useDensity,
             this->_includeUnderflowOverflow);

   addnorms(hists, vars, _reco, "sig_norm", "signal normalization systematic", kRed, this->_useDensity,
            this->_includeUnderflowOverflow);
   addnorms(hists, vars, _bkg, "bkg_norm", "background normalization systematic", kBlue, this->_useDensity,
            this->_includeUnderflowOverflow);

   std::sort(hists.begin(), hists.end(), [](const TH1 *ha, const TH1 *hb) { return ha->Integral() < hb->Integral(); });

   THStack *hstack = new THStack("breakdown", "breakdown");
   for (auto *h : hists) {
      hstack->Add(h);
   }
   hstack->SetTitle(TString(
      "Breakdown of systematic uncertainties on expected measurement;Observable;Variance / Number of Events Squared"));
   return hstack;
}

THStack *RooUnfoldSpec::makeTruthBreakdownHistogram() const
{
   std::vector<TH1 *> hists;
   std::vector<Variable<TH1>> vars;
   for (auto &obs : _obs_truth) {
      vars.push_back(Variable<TH1>(static_cast<RooRealVar *>(obs)));
   }

   auto *sig_nom = stathist(hists, vars, _truth, "sig_stat", "signal statistics", kRed, this->_useDensity,
                            this->_includeUnderflowOverflow);

   addshapes(hists, vars, _truth, "sig_shape", "signal shape systematic", kRed, this->_useDensity,
             this->_includeUnderflowOverflow);

   addnorms(hists, vars, _truth, "sig_norm", "signal normalization systematic", kRed, this->_useDensity,
            this->_includeUnderflowOverflow);

   std::sort(hists.begin(), hists.end(), [](const TH1 *ha, const TH1 *hb) { return ha->Integral() < hb->Integral(); });

   THStack *hstack = new THStack("breakdown", "breakdown");
   for (auto *h : hists) {
      hstack->Add(h);
   }
   hstack->SetTitle(
      TString("Breakdown of systematic uncertainties on truth model;Observable;Variance / Number of Events Squared"));
   return hstack;
}

void RooUnfoldSpec::addStatToCovarianceMatrix(const HistContainer &histContainer, TMatrixD &covarianceMatrix) const
{
   const int bins = covarianceMatrix.GetNrows();
   auto counts = h2v(histContainer._nom, false, this->_useDensity);
   for (size_t i = 0; i < bins; ++i) {
      // asusming sqrt(n) uncertainties, so the variance is the central value
      covarianceMatrix(i, i) += counts[i];
   }
}

void RooUnfoldSpec::addShapeToCovarianceMatrix(const HistContainer &cont, const ShapeSys &var,
                                               TMatrixD &covarianceMatrix) const
{
   const int bins = covarianceMatrix.GetNrows();
   if (var.size() != 2) {
      throw std::runtime_error("unable to process systematic with size " + std::to_string(var.size()) + " != 2");
   }

   // Calculate sum of squares of shape uncertainties
   auto nominal = h2v(cont._nom, false, this->_useDensity);
   auto shapeUp = h2v(var[0], false, this->_useDensity);
   auto shapeDown = h2v(var[1], false, this->_useDensity);
   for (int i = 0; i < bins; ++i) {
      double shapeUncertainty = std::pow(shapeUp[i] - nominal[i], 2) + std::pow(shapeDown[i] - nominal[i], 2);
      covarianceMatrix(i, i) += shapeUncertainty;
   }
}

void RooUnfoldSpec::addNormToCovarianceMatrix(const HistContainer &cont, const NormSys &var,
                                              TMatrixD &covarianceMatrix) const
{
   const int bins = covarianceMatrix.GetNrows();
   if (var.size() != 2) {
      throw std::runtime_error("unable to process systematic with size " + std::to_string(var.size()) + " != 2");
   }

   // Calculate sum of squares of norm uncertainties
   auto nominal = h2v(cont._nom, false, this->_useDensity);
   for (int i = 0; i < bins; ++i) {
      for (int i = 0; i < nominal.GetNrows(); ++i) {
         covarianceMatrix(i, i) += pow(nominal[i] * (var[0] - 1), 2) + pow(nominal[i] * (var[1] - 1), 2);
      }
   }
}

void RooUnfoldSpec::addToCovarianceMatrix(const HistContainer &histContainer, TMatrixD &covarianceMatrix) const
{
   this->addStatToCovarianceMatrix(histContainer, covarianceMatrix);

   // Process shape uncertainties
   if (!histContainer._shapes.empty()) {
      for (const auto &var : histContainer._shapes) {
         RooUnfoldSpec::addShapeToCovarianceMatrix(histContainer, var.second, covarianceMatrix);
      }
   }
   // Process norm uncertainties
   if (!histContainer._norms.empty()) {
      for (const auto &var : histContainer._norms) {
         RooUnfoldSpec::addNormToCovarianceMatrix(histContainer, var.second, covarianceMatrix);
      }
   }
}

namespace {
template <typename THist>
THist *
createHistogram(RooAbsReal *hist, bool useDensity, const std::string &name, const std::vector<Variable<THist>> &vars);

template <>
TH1 *createHistogram<TH1>(RooAbsReal *hist, bool useDensity, const std::string &name,
                          const std::vector<Variable<TH1>> &vars)
{
   // Create histogram for TH1
   return createHist<TH1>(h2v(hist, false, useDensity), name.c_str(), name.c_str(), vars, true);
}

template <>
TH2 *createHistogram<TH2>(RooAbsReal *hist, bool useDensity, const std::string &name,
                          const std::vector<Variable<TH2>> &vars)
{
   // Create histogram for TH2
   return createHist<TH2>(h2m(hist, false, useDensity), name.c_str(), name.c_str(), vars, true);
}

template <typename THist>
void addSampleToDictionary(const std::string &histName, const RooUnfoldSpec::HistContainer &cont,
                           std::map<std::string, TH1 *> &histograms, const std::vector<RooRealVar *> &observables,
                           bool useDensity)
{
   // Create variables for the histogram
   std::vector<Variable<THist>> vars;
   for (auto &obs : observables) {
      vars.push_back(Variable<THist>(obs));
   }

   histograms[histName] = createHistogram<THist>(cont._nom, useDensity, histName, vars);
   for (const auto &var : cont._shapes) {
      std::string shape_up = histName + "_shape_up_" + var.first;
      std::string shape_down = histName + "_shape_down_" + var.first;
      histograms[shape_up] = createHistogram<THist>(var.second[0], useDensity, shape_up, vars);
      histograms[shape_down] = createHistogram<THist>(var.second[1], useDensity, shape_down, vars);
   }
}
} // namespace

std::map<std::string, TH1 *> RooUnfoldSpec::createHistogramDictionary() const
{
   std::vector<RooRealVar *> observables;
   for (auto &obs : _obs_reco) {
      observables.push_back(static_cast<RooRealVar *>(obs));
   }

   std::map<std::string, TH1 *> histograms;
   addSampleToDictionary<TH1>("reco_sig", _reco, histograms, observables, _useDensity);
   addSampleToDictionary<TH1>("bkg", _bkg, histograms, observables, _useDensity);
   addSampleToDictionary<TH1>("truth_sig", _truth, histograms, observables, _useDensity);
   addSampleToDictionary<TH2>("response", _res, histograms, observables, _useDensity);
   return histograms;
}

std::string RooUnfoldSpec::createLikelihoodConfig() const
{
   auto tree = RooFit::Detail::JSONTree::create();
   auto &root = tree->rootnode();
   root.set_map();

   // Settings
   auto &settings = root["settings"].set_map();
   settings["include_systematics"] << false;
   settings["prune_systematics_threshold"] << 0;
   settings["prune_migration_threshold"] << 0;

   // Channels
   auto &channels = root["channels"].set_map();
   auto &channel = channels["reco_SR"].set_map();

   // Variables
   auto &variables = channel["variables"].set_seq();
   variables.append_child() << "reco_x";

   // Samples
   auto &samples = channel["samples"].set_map();

   auto addSampleToJSON = [&](const std::string &histString, const HistContainer &cont,
                              RooFit::Detail::JSONNode &target) {
      target["data"] << histString;

      auto &modifiers = target["modifiers"].set_seq();

      for (const auto &var : cont._shapes) {
         std::string shape_up = histString + "_shape_up_" + var.first;
         std::string shape_down = histString + "_shape_down_" + var.first;

         auto &mod = modifiers.append_child().set_map();
         mod["name"] << histString + "_shape_var_" + var.first;
         mod["type"] << "histosys";
         auto &data = mod["data"].set_map();
         data["hi"] << shape_up;
         data["lo"] << shape_down;
      }

      for (const auto &var : cont._norms) {
         std::string norm_up = histString + "_norm_up_" + var.first;
         std::string norm_down = histString + "_norm_down_" + var.first;

         auto &mod = modifiers.append_child().set_map();
         mod["name"] << histString + "_norm_var_" + var.first;
         mod["type"] << "normsys";
         auto &data = mod["data"].set_map();
         data["hi"] << norm_up;
         data["lo"] << norm_down;
      }
   };

   auto &reco_sig_sample = samples["signal"].set_map();
   addSampleToJSON("reco_sig", _reco, reco_sig_sample);
   auto &reco_bkg_sample = samples["background"].set_map();
   addSampleToJSON("bkg", _bkg, reco_bkg_sample);

   // Unfolding
   auto &unfolding = reco_sig_sample["unfolding"].set_map();
   auto &regularize = unfolding["regularize"].set_map();
   regularize["type"] << "tikhonov";
   regularize["strength"] << 0.1;
   regularize["curvature"] << "ss";
   unfolding["poi"] << "xs";
   unfolding["poitype"] << "cs";

   auto &truth = unfolding["truth"].set_map();
   addSampleToJSON("truth_sig", _truth, truth);

   auto &migration = unfolding["migration"].set_map();
   addSampleToJSON("response", _res, migration);

   // Convert the JSON configuration to a string
   std::stringstream ss;
   ss << root;
   return ss.str();
}

RooUnfolding::RooFitHist *RooUnfoldSpec::makeHistogram(const HistContainer &source, double /*errorThreshold*/)
{
   // build a new RooFitHist based on the source HistContainer. relative bin errors above errorThreshold will be
   // modelled as gamma parameters.
   RooAbsReal *hf = source._nom;
   RooAbsReal *func = hf;
   if (!hf)
      return 0;
   std::vector<RooRealVar *> obs;
   RooArgList obslist;
   if (this->_obs_all.getSize() < 1) {
      throw std::runtime_error("in RooUnfoldSpec::makeHistogram: no observables known!");
   }
   for (RooAbsArg *arg : this->_obs_all) {
      if (!arg)
         continue;
      if (!hf->dependsOn(*arg))
         continue;
      obs.push_back(dynamic_cast<RooRealVar *>(arg));
      obslist.add(*arg);
   }
   if (obslist.getSize() < 1) {
      std::stringstream ss;
      ss << "in RooUnfoldSpec::makeHistogram: function '";
      hf->printStream(ss, 0, RooPrintable::kStandard, "");
      ss << "' does not depend on any of the known observables '";
      this->_obs_all.printStream(ss, 0, RooPrintable::kStandard, "");
      ss << "'!";
      throw std::runtime_error(ss.str());
   }
   if (source._shapes.size() > 0) {
      RooArgList up, dn;
      RooArgList params;
      for (auto var : source._shapes) {
         TString sysname(var.first);
         if (var.second.size() != 2) {
            throw std::runtime_error(TString::Format("unable to process systematics '%s' with size %d != 2",
                                                     var.first.c_str(), (int)(var.second.size()))
                                        .Data());
         }
         up.add(*var.second[0]);
         dn.add(*var.second[1]);
         TString name = TString::Format("alpha_%s", var.first.c_str());
         RooRealVar *p = (RooRealVar *)(this->_alphas.find(name));
         if (!p) {
            p = new RooRealVar(name, name, 0, -5, 5);
            p->setError(1);
            this->addGaussNP(p);
         }
         params.add(*p);
      }
      TString name = TString::Format("%s_HistoSystematics", hf->GetName());
      hf->SetName(TString::Format("%s_Nominal", hf->GetName()));
      func = new PiecewiseInterpolation(name.Data(), name.Data(), *hf, up, dn, params);
   }
   RooArgList components;
   if (source._norms.size() > 0) {
      std::vector<double> up, dn;
      RooArgList params;
      for (auto var : source._norms) {
         TString sysname(var.first);
         up.push_back(var.second[0]);
         dn.push_back(var.second[1]);
         TString name = TString::Format("alpha_%s", var.first.c_str());
         RooRealVar *p = (RooRealVar *)(this->_alphas.find(name));
         if (!p) {
            p = new RooRealVar(name, name, 0, -5, 5);
            p->setError(1);
            this->addGaussNP(p);
         }
         params.add(*p);
      }
      TString name = TString::Format("%s_OverallSystematics", hf->GetName());
      components.add(*(new RooStats::HistFactory::FlexibleInterpVar(name.Data(), name.Data(), params, 1., up, dn)));
   }
   for (auto g : source._gammas) {
      this->addPoissonNP(g);
   }
   if (source._staterror)
      components.add(*source._staterror);
   if (components.getSize() > 0) {
      components.add(*func);
      TString name(hf->GetName());
      hf->SetName(TString::Format("%s_hist", hf->GetName()));
      func = new RooProduct(name.Data(), hf->GetTitle(), components);
   }
   return new RooUnfolding::RooFitHist(func, obs, source._gammas);
}

RooHistFunc *RooUnfoldSpec::makeHistFuncTruth(const TH1 *hist)
{
   //! create a new truth hist func
   if (!hist)
      return NULL;
   return RooUnfolding::makeHistFunc(hist, this->_obs_truth, this->_includeUnderflowOverflow, !this->_useDensity);
}

RooHistFunc *RooUnfoldSpec::makeHistFuncMeasured(const TH1 *hist)
{
   //! create a new measured hist func
   if (!hist)
      return NULL;
   return RooUnfolding::makeHistFunc(hist, this->_obs_reco, this->_includeUnderflowOverflow, !this->_useDensity);
}

RooProdPdf *RooUnfoldSpec::makeConstraints()
{
   //! create all the constraint terms
   RooArgList constraints;
   for (auto a : this->_alphas) {
      RooRealVar *p = (RooRealVar *)(a);
      RooRealVar *mean = new RooRealVar(TString::Format("%s_nom", p->GetName()),
                                        TString::Format("Mean value for %s", p->GetName()), p->getVal());
      RooRealVar *sigma = new RooRealVar(TString::Format("%s_sigma", p->GetName()),
                                         TString::Format("Sigma for %s", p->GetName()), p->getError());
      mean->setConstant(true);
      sigma->setConstant(true);
      RooGaussian *gaus =
         new RooGaussian(TString::Format("%s_constraint", p->GetName()),
                         TString::Format("Gaussian constraint term for %s", p->GetName()), *p, *mean, *sigma);
      constraints.add(*gaus);
   }
   for (auto g : this->_gammas) {
      RooRealVar *p = (RooRealVar *)(g);
      double n = pow(p->getError(), -2);
      RooRealVar *tau = new RooRealVar(TString::Format("%s_tau", p->GetName()),
                                       TString::Format("Tau parameter value for %s", p->GetName()), n);
      tau->setConstant(true);
      RooArgList params(*p, *tau);
      RooProduct *prod = new RooProduct(TString::Format("%s_nEvents", p->GetName()),
                                        TString::Format("Number of events for %s", p->GetName()), params);
      RooPoisson *pois = new RooPoisson(TString::Format("%s_constraint", p->GetName()),
                                        TString::Format("Poisson constraint term for %s", p->GetName()), *prod, *tau);
      constraints.add(*pois);
   }
   return new RooProdPdf(TString::Format("%s_constraints", this->GetName()), "Unfolding constraint terms", constraints);
}

void RooUnfoldSpec::addGaussNP(RooRealVar *v)
{
   //! add a new gaussian NP
   if (v)
      this->_alphas.add(*v);
}
void RooUnfoldSpec::addPoissonNP(RooRealVar *v)
{
   //! add a new poisson NP
   if (v)
      this->_gammas.add(*v);
}

void RooUnfoldSpec::makeBackground()
{
   //! create the background
   this->_locked = true;
   if (!this->_cache._bkg) {
      this->_cache._bkg = this->makeHistogram(this->_bkg, this->_errorThreshold);
   }
}
void RooUnfoldSpec::makeData()
{
   //! create the data
   this->_locked = true;
   if (!this->_cache._data) {
      this->_cache._data = this->makeHistogram(this->_data, 0);
   }
}
void RooUnfoldSpec::makeResponse()
{
   //! create the response
   this->_locked = true;
   if (!this->_cache._res) {
      this->makeReco();
      this->makeTruth();
      if (!this->_res._nom)
         throw std::runtime_error("no response input given!");
      this->_cache._res = this->makeHistogram(this->_res, this->_errorThreshold);
      this->_cache._response = new RooFitUnfoldResponse(this->GetName(), this->GetTitle(), this->_cache._res,
                                                        this->_cache._truth, this->_cache._reco, true);
   }
}
void RooUnfoldSpec::makeTruth()
{
   //! create the truth
   this->_locked = true;
   if (!this->_cache._truth) {
      if (!this->_truth._nom)
         throw std::runtime_error("no truth input given!");
      this->_cache._truth = this->makeHistogram(this->_truth, this->_errorThreshold);
   }
}
void RooUnfoldSpec::makeReco()
{
   //! create the reconstructed
   this->_locked = true;
   if (!this->_cache._reco) {
      if (!this->_reco._nom)
         throw std::runtime_error("no measured input given!");
      this->_cache._reco = this->makeHistogram(this->_reco, this->_errorThreshold);
   }
}
void RooUnfoldSpec::makeDataMinusBackground()
{
   //! create the data minus background
   this->_locked = true;
   this->makeData();
   if (!this->_cache._data_minus_bkg) {
      if (this->_bkg._nom) {
         this->makeResponse();
         this->makeBackground();
         this->_cache._data_minus_bkg =
            this->_cache._response->makeHistSum(this->_cache._data->func(), this->_cache._bkg->func(), 1., -1.);
      } else {
         this->_cache._data_minus_bkg = this->_cache._data;
      }
      TString name(TString::Format("%s_data_minus_bkg", this->GetName()));
      this->_cache._data_minus_bkg->func()->SetName(name);
      this->_cache._data_minus_bkg->func()->SetTitle(name);
   }
}

RooAbsReal *RooUnfoldSpec::getBackground()
{
   //! retrieve the background
   this->makeBackground();
   if (this->_cache._bkg)
      return this->_cache._bkg->func();
   else
      return NULL;
}
RooAbsReal *RooUnfoldSpec::getData()
{
   //! retrieve the background
   this->makeData();
   return this->_cache._data->func();
}
RooAbsReal *RooUnfoldSpec::getResponse()
{
   //! retrieve the response
   this->makeResponse();
   return this->_cache._res->func();
}
RooAbsReal *RooUnfoldSpec::getTruth()
{
   //! retrieve the truth
   this->makeTruth();
   return this->_cache._truth->func();
}
RooAbsReal *RooUnfoldSpec::getReco()
{
   //! retrieve the reconstructed
   this->makeReco();
   return this->_cache._reco->func();
}
RooAbsReal *RooUnfoldSpec::getDataMinusBackground()
{
   //! retrieve the data minus background
   this->makeDataMinusBackground();
   return this->_cache._data_minus_bkg->func();
}

RooUnfoldT<RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> *RooUnfoldSpec::unfold(Algorithm alg, Double_t regparam)
{
   //! create the unfolding object
   this->makeResponse();
   this->makeDataMinusBackground();

   RooUnfoldT<RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> *unfolding =
      RooUnfoldT<RooUnfolding::RooFitHist, RooUnfolding::RooFitHist>::New(alg, this->_cache._response,
                                                                          this->_cache._data_minus_bkg, regparam);
   if (this->_cache._bkg) {
      unfolding->SetBkg(this->_cache._bkg);
   }
   unfolding->SetMeasuredCov(this->makeCovarianceMatrix());
   unfolding->SetTruth(this->_cache._truth);
   unfolding->SetOverflow(this->_includeUnderflowOverflow);

   return unfolding;
}

void RooUnfoldSpec::HistContainer::setNominal(RooAbsReal *nom, const RooArgList &obslist)
{
   this->_nom = nom;
   this->_obs.add(obslist);
}

void RooUnfoldSpec::HistContainer::setNominal(RooDataHist *data, const RooArgList &obslist)
{
   this->_nom = RooUnfolding::makeHistFunc(data, obslist);
   this->_gammas = RooUnfolding::createGammas(data, obslist, 0.);
   this->_obs.add(obslist);
   if (_gammas.size() > 0) {
      this->_staterror = RooUnfolding::makeParamHistFunc(TString::Format("%s_staterrors", _nom->GetName()),
                                                         _nom->GetTitle(), obslist, _gammas);
   }
}

void RooUnfoldSpec::HistContainer::setNominal(const TH1 *nom, const RooArgList &obslist, double errorThreshold,
                                              bool includeUnderflowOverflow, bool correctDensity)
{
   // if useDensity is true, the inputs are already in density space - then we don't need to correct anymore
   this->_nom = RooUnfolding::makeHistFunc(nom, obslist, includeUnderflowOverflow, correctDensity);
   this->_obs.add(obslist);
   if (errorThreshold >= 0) {
      this->_gammas = RooUnfolding::createGammas(nom, includeUnderflowOverflow, errorThreshold);
      if (_gammas.size() > 0) {
         this->_staterror = RooUnfolding::makeParamHistFunc(TString::Format("%s_staterrors", nom->GetName()),
                                                            nom->GetTitle(), obslist, _gammas);
      }
   }
}

void RooUnfoldSpec::HistContainer::addShape(const char *name, RooAbsReal *up, RooAbsReal *dn)
{
   this->_shapes[name] = {up, dn};
}

void RooUnfoldSpec::HistContainer::addNorm(const char *name, double up, double dn)
{
   this->_norms[name] = {up, dn};
}

RooUnfoldSpec::HistContainer::~HistContainer() {}

void RooUnfoldSpec::lockCheck()
{
   if (this->_locked) {
      throw std::runtime_error("this instance of RooUnfoldSpec is locked - it has already been used to produce results "
                               "and can no longer be modified. please create a new instance for modifications!");
   }
}

void RooUnfoldSpec::checkConsistency(const HistContainer &cont, const TH1 *hist)
{
   int bins = countBins(cont._obs);
   if (bins != (hist->GetNbinsX() * hist->GetNbinsY() * hist->GetNbinsZ())) {
      throw std::runtime_error(std::string("unable to register systematic histogram '") + hist->GetName() +
                               ", number of bins does not match (" + std::to_string(bins) + " vs." +
                               std::to_string(hist->GetNbinsX()) + ")");
   }
}

void RooUnfoldSpec::registerSystematic(Contribution c, const char *name, const TH1 *up, const TH1 *down)
{
   //! register a new shape systematic
   this->lockCheck();
   // if useDensity is true, the inputs are already in density space - then we don't need to correct anymore
   switch (c) {
   case kTruth:
      this->checkConsistency(this->_truth, up);
      this->checkConsistency(this->_truth, down);
      this->_truth.addShape(
         name,
         RooUnfolding::makeHistFunc(TString::Format("truth_%s_%s_up", up->GetName(), name), up, this->_obs_truth,
                                    this->_includeUnderflowOverflow, !this->_useDensity),
         RooUnfolding::makeHistFunc(TString::Format("truth_%s_%s_dn", down->GetName(), name), down, this->_obs_truth,
                                    this->_includeUnderflowOverflow, !this->_useDensity));
      break;
   case kMeasured:
      this->checkConsistency(this->_reco, up);
      this->checkConsistency(this->_reco, down);
      this->_reco.addShape(
         name,
         RooUnfolding::makeHistFunc(TString::Format("meas_%s_%s_up", up->GetName(), name), up, this->_obs_reco,
                                    this->_includeUnderflowOverflow, !this->_useDensity),
         RooUnfolding::makeHistFunc(TString::Format("meas_%s_%s_dn", down->GetName(), name), down, this->_obs_reco,
                                    this->_includeUnderflowOverflow, !this->_useDensity));
      break;
   case kData:
      this->checkConsistency(this->_data, up);
      this->checkConsistency(this->_data, down);
      this->_data.addShape(
         name,
         RooUnfolding::makeHistFunc(TString::Format("data_%s_%s_up", up->GetName(), name), up, this->_obs_reco,
                                    this->_includeUnderflowOverflow, !this->_useDensity),
         RooUnfolding::makeHistFunc(TString::Format("data_%s_%s_dn", down->GetName(), name), down, this->_obs_reco,
                                    this->_includeUnderflowOverflow, !this->_useDensity));
      break;
   case kResponse:
      this->checkConsistency(this->_res, up);
      this->checkConsistency(this->_res, down);
      this->_res.addShape(
         name,
         RooUnfolding::makeHistFunc(TString::Format("resp_%s_%s_up", up->GetName(), name), up, this->_obs_all,
                                    this->_includeUnderflowOverflow, !this->_useDensity),
         RooUnfolding::makeHistFunc(TString::Format("resp_%s_%s_dn", down->GetName(), name), down, this->_obs_all,
                                    this->_includeUnderflowOverflow, !this->_useDensity));
      break;
   case kBackground:
      this->checkConsistency(this->_bkg, up);
      this->checkConsistency(this->_bkg, down);
      this->_bkg.addShape(
         name,
         RooUnfolding::makeHistFunc(TString::Format("bkg_%s_%s_up", up->GetName(), name), up, this->_obs_reco,
                                    this->_includeUnderflowOverflow, !this->_useDensity),
         RooUnfolding::makeHistFunc(TString::Format("bkg_%s_%s_dn", down->GetName(), name), down, this->_obs_reco,
                                    this->_includeUnderflowOverflow, !this->_useDensity));
      break;
   }
}

void RooUnfoldSpec::registerSystematic(Contribution c, const char *name, double up, double down)
{
   //! register a new normalization systematic
   this->lockCheck();
   switch (c) {
   case kTruth: this->_truth.addNorm(name, up, down); break;
   case kMeasured: this->_reco.addNorm(name, up, down); break;
   case kData: this->_data.addNorm(name, up, down); break;
   case kResponse: this->_res.addNorm(name, up, down); break;
   case kBackground: this->_bkg.addNorm(name, up, down); break;
   }
}

#ifdef NO_WRAPPERPDF
RooAbsPdf *RooUnfoldSpec::makePdf(Algorithm /*alg*/, Double_t /*regparam*/)
{
   throw std::runtime_error("need RooWrapperPdf to create unfolding Pdfs, upgrade ROOT version!");
   return NULL;
}
#else
RooAbsPdf *RooUnfoldSpec::makePdf(Algorithm alg, Double_t regparam)
{
   //! create an unfolding pdf
   RooUnfoldFunc *func = static_cast<RooUnfoldFunc *>(this->makeFunc(alg, regparam));
   func->setDensity(true);
   RooWrapperPdf *pdf =
      new RooWrapperPdf(TString::Format("%s_pdf", func->GetName()), TString::Format("%s Pdf", func->GetTitle()), *func);
   RooAbsReal *integral = func->createIntegral(this->_obs_truth);
   RooExtendPdf *extpdf = new RooExtendPdf(TString::Format("%s_extpdf", func->GetName()),
                                           TString::Format("%s Extended Pdf", func->GetTitle()), *pdf, *integral);

   RooProdPdf *constraints = this->makeConstraints();
   RooArgList comps(*extpdf, *constraints);
   RooProdPdf *prod = new RooProdPdf(TString::Format("%s_x_constraints", this->GetName()),
                                     "Unfolding pdf, including constraints", comps);
   prod->setStringAttribute("source", func->GetName());
   return prod;
}
#endif

RooAbsReal *RooUnfoldSpec::makeFunc(Algorithm alg, Double_t regparam)
{
   //! create an unfolding function
   auto *unfold = this->unfold(alg, regparam);
   RooAbsReal *func = new RooUnfoldFunc(this->GetName(), this->GetTitle(), unfold, false);
   func->setStringAttribute("source", func->GetName());
   return func;
}

ClassImp(RooUnfoldSpec)

#endif
