#ifndef ROOUNFOLDHELPERS_ROOABSREAL_HH
#define ROOUNFOLDHELPERS_ROOABSREAL_HH

#include "RooUnfoldHelpers.h"

#ifndef NOROOFIT

#include "RooRealVar.h"
#include "RooAbsData.h"
#include "TH1.h"

class RooHistFunc;
class RooHistPdf;

namespace RooUnfolding {
class RooFitHist : public TObject {
public:
   RooFitHist();
   RooFitHist(const RooFitHist *h);
   RooFitHist(RooAbsReal *f, const std::vector<RooRealVar *> &obs);
   RooFitHist(RooAbsReal *f, RooRealVar *obs);
   RooFitHist(RooAbsReal *f, RooRealVar *obs1, RooRealVar *obs2);
   RooFitHist(RooAbsReal *f, const std::vector<RooRealVar *> &obs, const std::vector<RooRealVar *> &nps);
   RooFitHist(RooAbsReal *f, RooRealVar *obs, const std::vector<RooRealVar *> &nps);
   RooFitHist(RooAbsReal *f, RooRealVar *obs1, RooRealVar *obs2, const std::vector<RooRealVar *> &nps);
   RooFitHist(RooHistFunc *f, const std::vector<RooRealVar *> &obs);
   RooFitHist(RooHistFunc *f, RooRealVar *obs);
   RooFitHist(RooHistFunc *f, RooRealVar *obs1, RooRealVar *obs2);
   RooFitHist(RooHistFunc *hist, const RooArgList &obslist, double uncThreshold = -1);
   RooFitHist(RooDataHist *hist, const std::vector<RooRealVar *> &obs);
   RooFitHist(RooDataHist *hist, RooRealVar *obs);
   RooFitHist(RooDataHist *hist, RooRealVar *obs1, RooRealVar *obs2);
   RooFitHist(RooDataHist *hist, const RooArgList &obslist, double uncThreshold = -1);
   RooFitHist(const TH1 *hist, const std::vector<RooRealVar *> &obs, bool includeUnderflowOverflow,
              double errorThreshold, bool correctDensity = false);
   RooFitHist(const TH1 *hist, const RooArgList &obs, bool includeUnderflowOverflow, double errorThreshold,
              bool correctDensity = false);
   virtual ~RooFitHist() = default;

   virtual const char *name() const;
   virtual const char *title() const;
   virtual RooRealVar *obs(size_t) const;
   virtual size_t dim() const;
   RooFitHist *clone() const;
   RooFitHist *asimovClone(bool correctDensity) const;
   RooFitHist *asimov1DClone(bool correctDensity, TVectorD &val, TVectorD &err) const;

   virtual double error() const;
   virtual double value() const;
   virtual bool weighted() const;
   virtual int bin() const;
   virtual const std::vector<RooRealVar *> &nps() const;

   virtual void printHist() const;
   virtual RooAbsReal *func() const;
   virtual bool checkValidity() const;
   virtual void replace(const RooAbsCollection &newlist);

   virtual void saveSnapshot(std::map<std::string, double> &snsh) const;
   virtual void loadSnapshot(const std::map<std::string, double> &snsh);

   virtual const char *GetName() const override;
   virtual const char *GetTitle() const override;
   virtual void Print(const char *opts = 0) const override;

protected:
   RooAbsReal *setupErrors(const RooHistFunc *hf, const RooDataHist *dh, double uncThreshold);
   RooAbsReal *_func;
   std::vector<RooRealVar *> _obs;
   std::vector<RooRealVar *> _gamma;
   ClassDefOverride(RooFitHist, 1)
};

template <>
struct Variable<RooFitHist> {
   RooRealVar *_var;
   Variable(RooRealVar *var);
   bool irregular() const;
   void print() const;
};

RooDataHist *convertTH1(const TH1 *histo, const std::vector<RooRealVar *> &vars, bool includeUnderflowOverflow,
                        bool correctDensity = false, double scale = 1.);
double getIntegral(const TH1 *histo, bool includeUnderflowOverflow, bool correctDensity);
RooDataHist *convertTH1(const TH1 *histo, const RooArgList &obs, bool includeUnderflowOverflow,
                        bool correctDensity = false, double scale = 1.);
std::vector<RooRealVar *> createGammas(const TH1 *histo, bool includeUnderflowOverflow, double uncThreshold);
std::vector<RooRealVar *> createGammas(const RooDataHist *dh, const RooArgList &obs, double uncThreshold);
RooAbsReal *makeParamHistFunc(const char *name, const char *title, const RooArgList &obslist,
                              const std::vector<RooRealVar *> &gamma);
const RooArgSet *getObservables(const RooHistFunc *f);
void setGammaUncertainties(RooWorkspace *ws);
RooHistFunc *
makeHistFunc(const TH1 *histo, const RooArgList &obs, bool includeUnderflowOverflow, bool correctDensity = false);
RooHistFunc *makeHistFunc(const char *name, const TH1 *histo, const RooArgList &obs, bool includeUnderflowOverflow,
                          bool correctDensity = false);
RooHistFunc *makeHistFunc(RooDataHist *dhist, const std::vector<RooRealVar *> &obs);
RooHistFunc *makeHistFunc(RooDataHist *dhist, const RooArgList &obs);
RooAbsPdf *
makeHistPdf(const TH1 *histo, const RooArgList &obs, bool includeUnderflowOverflow, bool correctDensity = false);
RooAbsPdf *makeHistPdf(const char *name, const TH1 *histo, const RooArgList &obs, bool includeUnderflowOverflow,
                       bool correctDensity = false);
RooAbsPdf *makeHistPdf(RooDataHist *dhist, const std::vector<RooRealVar *> &obs);
void importToWorkspace(RooWorkspace *ws, RooAbsReal *object);
void importToWorkspace(RooWorkspace *ws, RooAbsData *object);
RooArgSet allVars(RooWorkspace *ws, const char *pattern);
template <class ObjT, class ListT>
std::vector<ObjT *> matchingObjects(const ListT *c, const char *pattern);
std::vector<RooAbsReal *> matchingObjects(const RooAbsCollection *c, const char *pattern);
void printClients(const RooAbsArg *obj);
void printServers(const RooAbsArg *obj);
TH1 *convertTH1(const TVectorD &values, const TVectorD &errors, const RooUnfolding::RooFitHist *hist);
TH1 *convertTH1(const TVectorD &values, const RooUnfolding::RooFitHist *hist);
TH2 *convertTH2(const TMatrixD &values, const TMatrixD &errors, const RooUnfolding::RooFitHist *hist);

TH1 *get_unfolded_histogram(const RooArgList &parameters, const TH1 *hist_template,
                            const std::string &parameter_template);
TH1 *get_unfolded_histogram(RooWorkspace &workspace, const TH1 *hist_template, const std::string &parameter_template);
} // namespace RooUnfolding

#endif

#endif
