/*===========================================================================*/
/*
 * BEGIN ROOUNFOLD COPYRIGHT
 * RooUnfold — Unfolding library for particle-physics inverse problems
 *
 * Copyright © 2019–2025 CERN and the authors’ respective research institutions
 * Please refer to the CONTRIBUTORS file for details.
 *
 * License: BSD-3-Clause
 * SPDX-License-Identifier: BSD-3-Clause
 *
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

#include "TVectorD.h"
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

namespace {
std::vector<std::vector<double>> getBinCenters(const RooArgList &vars)
{
   size_t numVars = vars.getSize();
   std::vector<std::vector<double>> result;
   std::vector<size_t> counters(numVars, 0);

   while (true) {
      std::vector<double> binCenter(numVars);
      bool done = true;

      for (size_t varIdx = 0; varIdx < numVars; ++varIdx) {
         auto *var = dynamic_cast<RooRealVar *>(vars.at(varIdx));
         if (!var)
            throw std::runtime_error("Non-RooRealVar object in RooArgList.");

         size_t numBins = var->numBins();
         binCenter[varIdx] = var->getBinning().binCenter(counters[varIdx]);

         if (counters[varIdx] < numBins - 1)
            done = false;
      }

      result.push_back(std::move(binCenter));

      if (done)
         break;

      for (size_t varIdx = 0; varIdx < numVars; ++varIdx) {
         auto *var = dynamic_cast<RooRealVar *>(vars.at(varIdx));
         if (++counters[varIdx] < var->numBins()) {
            break;
         } else {
            counters[varIdx] = 0;
         }
      }
   }

   return result;
}
void setValues(const RooArgList &vars, const std::vector<double> &values)
{
   for (size_t i = 0; i < vars.size(); ++i) {
      static_cast<RooRealVar *>(&vars[i])->setVal(values[i]);
   }
}
std::vector<double> v2v(const TVectorD &v)
{
   return std::vector<double>(v.GetMatrixArray(), v.GetMatrixArray() + v.GetNrows());
}
} // namespace

namespace {
template <typename K, typename V>
void expand_key_set(std::unordered_set<std::string> &key_set, const std::map<const K, V> &new_map)
{
   for (const auto &[key, _] : new_map) {
      key_set.insert(key);
   }
}
bool isCloseTo(double val, double ref, double tolerance)
{
   return std::fabs(val - ref) <= tolerance;
}
bool isCloseToOne(double val, double tolerance)
{
   return isCloseTo(val, 1., tolerance);
}
void clearIfCloseToOne(std::vector<double> &vec, double tolerance)
{
   for (double val : vec) {
      if (!isCloseToOne(val, tolerance)) {
         return; // Found a value that deviates too much — do nothing
      }
   }
   vec.clear(); // All values are within tolerance — clear the vector
}
double sanitize_ratio(double ratio, double bineff)
{
   const double a = 0.001; // full suppression below this
   const double b = 0.1;   // no effect above this

   double w;
   if (bineff <= a) {
      w = 0.0;
   } else if (bineff >= b) {
      w = 1.0;
   } else {
      double t = (bineff - a) / (b - a);
      w = t * t * (3 - 2 * t); // smoothstep
   }

   return 1.0 + (ratio - 1.0) * w;
}

std::pair<double, std::vector<double>>
factorize_sys(const std::vector<double> &nominal, const std::vector<double> &variation, bool sanitize)
{
   if (nominal.size() != variation.size()) {
      throw std::invalid_argument("Input vectors must have the same length.");
   }

   double sum_nominal = std::accumulate(nominal.begin(), nominal.end(), 0.0);
   double sum_variation = std::accumulate(variation.begin(), variation.end(), 0.0);

   std::vector<double> component_ratios;
   component_ratios.reserve(nominal.size());
   if (sum_nominal == 0.0 || sum_variation == 0.0) {
      return {1., component_ratios};
   }

   double total_ratio = sum_variation / sum_nominal;
   for (size_t i = 0; i < nominal.size(); ++i) {
      const double v = variation[i] / sum_variation * sum_nominal;
      const double n = nominal[i];
      double ratio = v / n;
      const double bineff = fabs(n / sum_nominal * nominal.size());
      if (n == 0.) {
         ratio = 1.;
      }
      component_ratios.push_back(sanitize ? sanitize_ratio(ratio, bineff) : ratio);
   }
   clearIfCloseToOne(component_ratios, 1e-6);
   std::vector<double> components;
   return {total_ratio, component_ratios};
}
void fill_vector(RooFit::Detail::JSONNode &node, const std::vector<double> &v)
{
   node.set_seq();
   for (size_t i = 0; i < v.size(); ++i) {
      node.append_child() << v[i];
   }
}
void fill_vector_product(RooFit::Detail::JSONNode &node, const std::vector<double> &v, const std::vector<double> &w,
                         bool flip = false)
{
   node.set_seq();
   for (size_t i = 0; i < v.size(); ++i) {
      node.append_child() << (flip ? w[i] * (2. - v[i]) : w[i] * v[i]);
   }
}
} // namespace

std::string RooUnfoldSpec::createLikelihoodJSON(double tau, bool include_sys, bool xs_pois, bool sanitize_sys) const
{
   // Create the JSON tree and set up the root node
   using namespace RooFit::Detail;
   auto tree = JSONTree::create();
   auto &root = tree->rootnode();
   root.set_map();

   // Add metadata
   root["metadata"].set_map()["hs3_version"] << "0.2";

   // Initialize JSON structure elements
   auto &distributions = root["distributions"].set_seq();
   auto &functions = root["functions"].set_seq();
   auto &domains = root["domains"].set_seq();
   auto &parameter_points = root["parameter_points"].set_seq();
   auto &data = root["data"].set_seq();
   auto &likelihoods = root["likelihoods"].set_seq();
   auto &analyses = root["analyses"].set_seq();

   auto writeObservable = [&](RooRealVar *r, JSONNode &axes) {
      auto &obs = axes.append_child().set_map();
      obs["name"] << r->GetName();
      if (!r->getBinning().isUniform()) {
         auto &bounds = obs["edges"].set_seq();
         bounds.append_child() << r->getBinning().binLow(0);
         for (size_t i = 0; i < r->numBins(); ++i) {
            bounds.append_child() << r->getBinning().binHigh(i);
         }
      } else {
         obs["nbins"] << r->getBinning().numBins();
         obs["min"] << r->getMin();
         obs["max"] << r->getMax();
      }
   };
   auto writeHistogram = [&](auto &h, JSONNode &node) {
      node.set_map();
      fill_vector(node["contents"], v2v(h2v(h, false, _useDensity)));
   };
   auto writeAxes = [&](JSONNode &axes) {
      for (auto &obs : this->_obs_reco) {
         writeObservable(static_cast<RooRealVar *>(obs), axes);
      }
   };
   auto writeParameter = [&](JSONNode &paramlist, const std::string &name, double value, bool constant = false) {
      auto &p = paramlist.append_child().set_map();
      p["name"] << name;
      p["value"] << value;
      if (constant) {
         p["const"] << 1;
      }
   };
   auto writeDomain = [&](JSONNode &paramlist, const std::string &name, double min, double max, bool constant = false) {
      auto &p = paramlist.append_child().set_map();
      p["name"] << name;
      p["min"] << min;
      p["max"] << max;
      if (constant) {
         p["const"] << 1;
      }
   };
   auto binVolume = [&](const RooArgList &observables) {
      double v = 1;
      for (auto *x : observables) {
         RooRealVar *obs = static_cast<RooRealVar *>(x);
         const int binIdx = obs->getBinning().binNumber(obs->getVal());
         v *= obs->getBinning().binWidth(binIdx);
      }
      return v;
   };

   // Add likelihood and analyses
   auto &simultaneous_likelihood = likelihoods.append_child().set_map();
   simultaneous_likelihood["name"] << "simultaneous_likelihood";
   simultaneous_likelihood["data"].set_seq();
   simultaneous_likelihood["distributions"].set_seq();
   simultaneous_likelihood["aux_distributions"].set_seq();

   auto &unfolded_analysis = analyses.append_child().set_map();
   std::string model_name = "unfolded_model";
   unfolded_analysis["name"] << model_name;
   unfolded_analysis["likelihood"] << "simultaneous_likelihood";
   auto &analysis_domains = unfolded_analysis["domains"].set_seq();
   auto &pois = unfolded_analysis["parameters_of_interest"].set_seq();

   auto createDomain = [&](const std::string &name, bool attach) -> JSONNode & {
      auto &domain = domains.append_child().set_map();
      if (attach)
         analysis_domains.append_child() << name;
      domain["name"] << name;
      domain["type"] << "product_domain";
      return domain["axes"].set_seq();
   };

   // Add domains and parameter points
   auto &default_axes = createDomain("default_domain", false);
   auto &poi_axes = createDomain(model_name + "_parameters_of_interest", true);
   auto &obs_axes = createDomain(model_name + "_observables", true);
   auto &globs_axes = createDomain(model_name + "_global_observables", true);
   auto &np_axes = createDomain(model_name + "_nuisance_parameters", true);

   for (auto &v : this->_obs_reco) {
      RooRealVar *obs = static_cast<RooRealVar *>(v);
      writeDomain(obs_axes, obs->GetName(), obs->getMin(), obs->getMax());
   }

   auto &default_values = parameter_points.append_child().set_map();
   default_values["name"] << "default_values";
   auto &default_parameters = default_values["parameters"].set_seq();

   auto &background_values = parameter_points.append_child().set_map();
   background_values["name"] << "background_only";
   auto &background_parameters = background_values["parameters"].set_seq();

   // add distributions
   auto &signalregion = distributions.append_child().set_map();
   signalregion["type"] << "histfactory_dist";
   writeAxes(signalregion["axes"].set_seq());
   auto &pdf_name = "reco_SR";
   signalregion["name"] << pdf_name;
   simultaneous_likelihood["distributions"].append_child() << pdf_name;
   auto &samples = signalregion["samples"].set_seq();

   std::unordered_set<std::string> np_names;

   // add the background sample
   std::unordered_set<std::string> background_systematics;
   if (include_sys) {
      expand_key_set(background_systematics, _bkg._shapes);
   }
   if (_bkg._nom) {
      auto &background = samples.append_child().set_map();
      background["name"] << "background";
      writeHistogram(_bkg._nom, background["data"]);
      auto &modifiers = background["modifiers"].set_seq();
      auto &bkg_nf = modifiers.append_child().set_map();
      bkg_nf["type"] << "normfactor";
      bkg_nf["name"] << "bkg_norm";
      writeParameter(default_parameters, "bkg_norm", 1., true);
      writeDomain(np_axes, "bkg_norm", 0, 50.);
      auto nominal = v2v(h2v(_bkg._nom, false, _useDensity));
      double sum_nom = std::accumulate(nominal.begin(), nominal.end(), 0.0);
      for (const auto &k : background_systematics) {
         const auto &shape = *_bkg._shapes.find(k);
         if (!shape.second.size() == 2) {
            throw std::runtime_error("cannot deal with bkg systematics that are not up/down");
         }
         auto up = factorize_sys(v2v(h2v(shape.second[0], false, _useDensity)), nominal, sanitize);
         auto dn = factorize_sys(v2v(h2v(shape.second[1], false, _useDensity)), nominal, sanitize);
         std::string pname;
         if (!isCloseToOne(up.first, 1e-6) || !isCloseToOne(dn.first, 1e-6)) {
            auto &normsys = modifiers.append_child().set_map();
            normsys["name"] << shape.first;
            normsys["type"] << "normsys";
            pname = "alpha_" + shape.first;
            auto &normdata = normsys["data"].set_map();
            normdata["hi"] << up.first;
            normdata["lo"] << dn.first;
            normdata["parameter"] << pname;
            normsys["constraint_name"] << pname + "Constraint";
         }
         if (up.second.size() > 0 || dn.second.size() > 0) {
            auto &shapesys = modifiers.append_child().set_map();
            shapesys["name"] << shape.first;
            shapesys["type"] << "histosys";
            pname = "alpha_" + shape.first;
            auto &shapedata = shapesys["data"].set_map();
            if (up.second.size() > 0)
               fill_vector_product(shapedata["hi"].set_map()["contents"], up.second, nominal);
            else
               fill_vector_product(shapedata["hi"].set_map()["contents"], dn.second, nominal, true);
            if (dn.second.size() > 0)
               fill_vector_product(shapedata["lo"].set_map()["contents"], dn.second, nominal);
            else
               fill_vector_product(shapedata["hi"].set_map()["contents"], up.second, nominal, true);
            shapedata["parameter"] << pname;
            shapesys["constraint_name"] << pname + "Constraint";
         }
         if (!pname.empty())
            np_names.insert(pname);
      }
   }

   // this is the part where the actual unfolding happens
   const auto &truth_bins = getBinCenters(_obs_truth);
   const auto &reco_bins = getBinCenters(_obs_reco);
   struct Folding {
      std::vector<std::vector<double>> response;
      std::vector<double> fiducial;
   };

   auto calculate_folding = [&](RooAbsReal *truth, RooAbsReal *res, RooAbsReal *reco) {
      Folding result;
      for (const auto &truth_bin : truth_bins) {
         setValues(_obs_truth, truth_bin);
         double fiducial_yield = truth->getVal() * binVolume(_obs_truth);
         result.fiducial.push_back(fiducial_yield);
         result.response.push_back(std::vector<double>());
      }
      result.response.push_back(std::vector<double>()); // this is for fakes
      std::vector<double> fakes;
      double total_fakes = 0;
      for (const auto &reco_bin : reco_bins) {
         setValues(_obs_reco, reco_bin);
         double observed_yield = reco->getVal() * binVolume(_obs_reco);
         double sum_fiducial_yields = 0;
         int i_truth = 0;
         for (const auto &truth_bin : truth_bins) {
            setValues(_obs_truth, truth_bin);
            double fiducial_yield = truth->getVal() * binVolume(_obs_truth);
            double response_yield = res->getVal() * binVolume(_obs_truth) * binVolume(_obs_reco);
            double response_function = response_yield / fiducial_yield;
            result.response[i_truth].push_back(xs_pois ? response_function : response_function * fiducial_yield);
            sum_fiducial_yields += response_function * fiducial_yield;
            ++i_truth;
         }
         double fake_yield = std::max(0., observed_yield - sum_fiducial_yields);
         total_fakes += fake_yield;
         fakes.push_back(fake_yield);
      }
      result.fiducial.push_back(total_fakes);
      for (size_t i = 0; i < reco_bins.size(); ++i) {
         result.response[truth_bins.size()].push_back(xs_pois ? fakes[i] / total_fakes : fakes[i]);
      }
      return result;
   };

   auto nominal_folding = calculate_folding(_truth._nom, _res._nom, _reco._nom);

   std::unordered_set<std::string> signal_systematics;
   if (include_sys) {
      expand_key_set(signal_systematics, _truth._shapes);
      expand_key_set(signal_systematics, _res._shapes);
      expand_key_set(signal_systematics, _reco._shapes);
   }
   std::map<std::string, Folding> up_systematics;
   std::map<std::string, Folding> dn_systematics;
   for (const auto &k : signal_systematics) {
      auto *truth_up = _truth._nom;
      auto *reco_up = _reco._nom;
      auto *res_up = _res._nom;
      auto *truth_dn = _truth._nom;
      auto *reco_dn = _reco._nom;
      auto *res_dn = _res._nom;
      auto truth_sys = _truth._shapes.find(k);
      auto reco_sys = _reco._shapes.find(k);
      auto res_sys = _res._shapes.find(k);
      if (truth_sys != _truth._shapes.end()) {
         if (truth_sys->second.size() != 2)
            throw std::runtime_error("cannot deal with truth systematics that are not up/down");
         truth_up = truth_sys->second[0];
         truth_dn = truth_sys->second[1];
      }
      if (reco_sys != _reco._shapes.end()) {
         if (reco_sys->second.size() != 2)
            throw std::runtime_error("cannot deal with reco systematics that are not up/down");
         reco_up = reco_sys->second[0];
         reco_dn = reco_sys->second[1];
      }
      if (res_sys != _res._shapes.end()) {
         if (res_sys->second.size() != 2)
            throw std::runtime_error("cannot deal with resolution systematics that are not up/down");
         res_up = res_sys->second[0];
         res_dn = res_sys->second[1];
      }
      up_systematics[k] = calculate_folding(truth_up, res_up, reco_up);
      dn_systematics[k] = calculate_folding(truth_dn, res_dn, reco_dn);
   }

   // add everything to the json structure
   int i_truth = 0;
   std::vector<std::string> poi_names, poi_nomnames;
   for (size_t i_truth = 0; i_truth <= truth_bins.size(); ++i_truth) {
      auto &signal = samples.append_child().set_map();
      std::string truth_cat = "fakes";
      if (i_truth < truth_bins.size())
         truth_cat = TString::Format("bin_%d", i_truth + 1).Data();
      signal["name"] << "signal_" + truth_cat;

      auto &data = signal["data"].set_map();
      const auto &nominal = nominal_folding.response[i_truth];
      fill_vector(data["contents"], nominal);

      std::string poi = (xs_pois ? "xs_" : "mu_") + truth_cat;
      std::string poi_nom = "nom_" + poi;

      auto &modifiers = signal["modifiers"].set_seq();

      auto &xs = modifiers.append_child().set_map();
      pois.append_child() << poi;
      xs["name"] << poi;
      xs["type"] << "normfactor";
      for (const auto &k : signal_systematics) {
         auto up = factorize_sys(nominal, up_systematics[k].response[i_truth], sanitize_sys);
         auto dn = factorize_sys(nominal, dn_systematics[k].response[i_truth], sanitize_sys);
         std::string pname;
         if (!isCloseToOne(up.first, 1e-6) || !isCloseToOne(dn.first, 1e-6)) {
            auto &normsys = modifiers.append_child().set_map();
            normsys["name"] << k;
            normsys["type"] << "normsys";
            pname = "alpha_" + k;
            auto &normdata = normsys["data"].set_map();
            normdata["hi"] << up.first;
            normdata["lo"] << dn.first;
            normdata["parameter"] << pname;
            normsys["constraint_name"] << pname + "Constraint";
         }
         if (up.second.size() > 0 || dn.second.size() > 0) {
            auto &shapesys = modifiers.append_child().set_map();
            shapesys["name"] << k;
            shapesys["type"] << "histosys";
            pname = "alpha_" + k;
            auto &shapedata = shapesys["data"].set_map();
            if (up.second.size() > 0)
               fill_vector_product(shapedata["hi"].set_map()["contents"], up.second, nominal);
            else
               fill_vector_product(shapedata["hi"].set_map()["contents"], dn.second, nominal, true);
            if (dn.second.size() > 0)
               fill_vector_product(shapedata["lo"].set_map()["contents"], dn.second, nominal);
            else
               fill_vector_product(shapedata["hi"].set_map()["contents"], up.second, nominal, true);
            shapedata["parameter"] << pname;
            shapesys["constraint_name"] << pname + "Constraint";
         }
         if (!pname.empty())
            np_names.insert(pname);
      }

      double poival = xs_pois ? nominal_folding.fiducial[i_truth] : 1.;
      writeParameter(default_parameters, poi, poival);
      writeParameter(background_parameters, poi, 0);
      writeParameter(default_parameters, poi_nom, poival, true);
      writeDomain(poi_axes, poi, 0, fabs(10 * poival));
      writeDomain(globs_axes, poi_nom, 0, fabs(10 * poival));
      poi_names.push_back(poi);
      poi_nomnames.push_back(poi_nom);
   }

   for (const auto &np : np_names) {
      auto &constraint = distributions.append_child().set_map();
      std::string glob = "nom_" + np;
      constraint["type"] << "gaussian_dist";
      constraint["mean"] << glob;
      constraint["x"] << np;
      constraint["sigma"] << 1.;
      constraint["name"] << np + "Constraint";
      writeParameter(default_parameters, glob, 0., true);
      writeParameter(default_parameters, np, 0.);
      writeDomain(np_axes, np, -5., 5.);
      writeDomain(globs_axes, glob, -5., 5.);
   }

   // create regularization term
   if (tau > 0) {
      std::stringstream tikhonov_differential;
      for (size_t i = 1; i < poi_names.size() - 2; ++i) { // skip fakes
         if (i > 1)
            tikhonov_differential << "+";
         tikhonov_differential << "pow(" << poi_names[i - 1] << "/" << poi_nomnames[i - 1] << "-2*" << poi_names[i]
                               << "/" << poi_nomnames[i] << "+" << poi_names[i + 1] << "/" << poi_nomnames[i + 1]
                               << ",2)";
      }
      auto &regularization_dist = distributions.append_child().set_map();
      regularization_dist["type"] << "exponential_dist";
      regularization_dist["x"] << "tikhonov_differential";
      regularization_dist["c"] << "tau";
      regularization_dist["name"] << "regularization";
      auto &regularization_func = functions.append_child().set_map();
      regularization_func["type"] << "generic_function";
      regularization_func["name"] << "tikhonov_differential";
      regularization_func["expression"] << tikhonov_differential.str().c_str();
      simultaneous_likelihood["aux_distributions"].append_child() << "regularization";
      auto &tau_val = default_parameters.append_child().set_map();
      tau_val["name"] << "tau";
      tau_val["value"] << tau;
      ;
      tau_val["const"] << 1;
   }

   // Add data
   std::string data_name = "measured_data";
   auto &data_entry = data.append_child().set_map();
   writeHistogram(_data._nom, data_entry);
   writeAxes(data_entry["axes"].set_seq());
   data_entry["name"] << data_name;
   data_entry["type"] << "binned";
   simultaneous_likelihood["data"].append_child() << data_name;

   // Serialize the tree to a string
   std::stringstream ss;
   root.writeJSON(ss);
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
   case kSignalTruth:
      this->checkConsistency(this->_truth, up);
      this->checkConsistency(this->_truth, down);
      this->_truth.addShape(
         name,
         RooUnfolding::makeHistFunc(TString::Format("truth_%s_%s_up", up->GetName(), name), up, this->_obs_truth,
                                    this->_includeUnderflowOverflow, !this->_useDensity),
         RooUnfolding::makeHistFunc(TString::Format("truth_%s_%s_dn", down->GetName(), name), down, this->_obs_truth,
                                    this->_includeUnderflowOverflow, !this->_useDensity));
      break;
   case kSignalMeasured:
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
   case kSignalResponse:
      this->checkConsistency(this->_res, up);
      this->checkConsistency(this->_res, down);
      this->_res.addShape(
         name,
         RooUnfolding::makeHistFunc(TString::Format("resp_%s_%s_up", up->GetName(), name), up, this->_obs_all,
                                    this->_includeUnderflowOverflow, !this->_useDensity),
         RooUnfolding::makeHistFunc(TString::Format("resp_%s_%s_dn", down->GetName(), name), down, this->_obs_all,
                                    this->_includeUnderflowOverflow, !this->_useDensity));
      break;
   case kBackgroundMeasured:
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
   case kSignalTruth: this->_truth.addNorm(name, up, down); break;
   case kSignalMeasured: this->_reco.addNorm(name, up, down); break;
   case kData: this->_data.addNorm(name, up, down); break;
   case kSignalResponse: this->_res.addNorm(name, up, down); break;
   case kBackgroundMeasured: this->_bkg.addNorm(name, up, down); break;
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
