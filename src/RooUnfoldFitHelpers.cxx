#include "RooUnfoldHelpers.h"

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"

namespace {
  RooRealVar* makeVar(const RooUnfolding::Variable& x){
    return new RooRealVar(x._name,x._name,x._nBins,x._min,x._max);
  }
}


namespace RooUnfolding {
  template<> void reset<RooAbsReal>(RooAbsReal* r){
    // TODO
  }
  template<> int findBin<RooAbsReal>(const RooAbsReal* h, double x, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  template<> int findBin<RooAbsReal>(const RooAbsReal* h, Double_t x){
    // TODO
    return 0;
  }
  template<> int findBin<RooAbsReal>(const RooAbsReal* h, Double_t x, Double_t y){
    // TODO
    return 0;
  }
  template<> int findBin<RooAbsReal>(const RooAbsReal* h, Double_t x, Double_t y, Double_t z){
    // TODO
    return 0;
  }

  template<> double min<RooAbsReal>(const RooAbsReal* hist, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  template<> double max<RooAbsReal>(const RooAbsReal* hist, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  template<> int sumW2N<RooAbsReal>(const RooAbsReal* hist){
    return 0;
  }
  template<> void add<RooAbsReal>(RooAbsReal* hista, const RooAbsReal* histb){
    // TODO
  }  
  template<> void projectY<RooAbsReal>(RooAbsReal* _res, RooAbsReal* _tru, bool overflow){
    // TODO
  } 
  template<> void projectX<RooAbsReal>(RooAbsReal* _res, RooAbsReal* _mes, bool overflow){
    // TODO
  }  
  template<> void subtractProjectX<RooAbsReal>(RooAbsReal* _res, RooAbsReal* _mes, RooAbsReal* _fak, bool overflow){
    // TODO
  }
  template<> int fill<RooAbsReal>(RooAbsReal* hist, double x, double w){
    return 0;
  }
  template<> int fill<RooAbsReal>(RooAbsReal* hist, double x, double y, double w){
    return 0;
  }  
  template<> int fill<RooAbsReal>(RooAbsReal* hist, double x, double y, double z, double w){
    return 0;
  }  
  template<> RooAbsReal* copy<RooAbsReal>(const RooAbsReal* r, bool reset, const char* name, const char* title){
    RooAbsReal* retval = (RooAbsReal*)(r->clone());
    if(name) retval->SetName(name);
    if(title) retval->SetTitle(title);
    return retval;
  }
  template<> int entries<RooAbsReal>(const RooAbsReal* hist){
    // TODO
    return 0;
  }
  template<> int dim<RooAbsReal>(const RooAbsReal* hist){
    // TODO
    return 0;
  }
  template<> int nBins<RooAbsReal>(const RooAbsReal* hist, bool overflow){
    // TODO
    return 0;
  }
  template<> int nBins<RooAbsReal>(const RooAbsReal* hist, RooUnfolding::Dimension d, bool overflow){
    // TODO
    return 0;
  }
  template<> double binCenter<RooAbsReal>(const RooAbsReal*h, int i, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  template<> double binWidth<RooAbsReal>(const RooAbsReal*h, int i, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  template<> double binHighEdge<RooAbsReal>(const RooAbsReal*h, int i, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  template<> double binLowEdge<RooAbsReal>(const RooAbsReal*h, int i, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  template<> void binXYZ<RooAbsReal>(const RooAbsReal* tru, int i, int& jx, int& jy, int& jz){
    // TODO
  }
  template<> double binError<RooAbsReal>(const RooAbsReal* h, Int_t i, Bool_t overflow)
  {
    // Bin error   by vector index
    // TODO
    return 0;
  }  
  template<> double binContent<RooAbsReal> (const RooAbsReal* h, Int_t i, Bool_t overflow){
    // TODO
    return 0;
  }
  template<class Hist2D> Hist2D* createHist(const char* name, const char* title, const Variable& x, const Variable& y);
  template<> RooAbsReal* createHist<RooAbsReal>(const char* name, const char* title, const Variable& x, const Variable& y){
    RooRealVar* rx = ::makeVar(x);
    RooRealVar* ry = ::makeVar(y);
    RooArgSet vars(*rx,*ry);
    RooDataHist* hist = new RooDataHist (name,title,vars);
    return new RooHistFunc(name,title,vars,vars,*hist);
  }
  template<> RooAbsReal* createHist<RooAbsReal>(const char* name, const char* title, const std::vector<Variable>& x){
    return NULL;
  }
  template<> RooAbsReal* createHist<RooAbsReal>(const TMatrixD& m, const char* name, const char* title, const Variable& x, const Variable& y){  
    // Sets the bin content of the histogram as that element of the input vector
    // TODO
    return NULL;
  }
  template<> RooAbsReal* createHist<RooAbsReal>(const TMatrixD& m, const TMatrixD& me, const char* name, const char* title, const Variable& x, const Variable& y){  
    // Sets the bin content of the histogram as that element of the input vector
    // TODO
    return NULL;
  }

  template<> RooAbsReal* h2h1d<RooAbsReal>(const RooAbsReal* h, int nb){
    // TODO
    return 0;
  }
  template<> RooAbsReal* copyHistogram<RooAbsReal>(const RooAbsReal* h, bool includeOverflow){
    // TODO
    return 0;
  }
  template<> RooAbsReal* createHist<RooAbsReal>(const TVectorD& v, const char* name, const char* title, const std::vector<Variable>& x, bool overflow){
    // TODO
    return 0;
  }
  template<> RooAbsReal* createHist<RooAbsReal>(const TVectorD& v, const TVectorD& ve, const char* name, const char* title, const std::vector<Variable>& x, bool overflow){
    // TODO
    return 0;
  }
  template<> const char* varname<RooAbsReal>(const RooAbsReal* h, Dimension d){
    // TODO
    return "";
  }
  template<> RooAbsReal* histNoOverflow<RooAbsReal>(const RooAbsReal* h, bool overflow){
    // overflow not supported anyway, do nothing
    // TODO
    return NULL;
  }

  template<> TVectorD subtract<RooAbsReal,TVectorD>(const TVectorD& orig, const RooAbsReal* hist, bool overflow) {
    // TODO
    TVectorD res;
    return res;
  }
  template<> void printTable<RooAbsReal> (std::ostream& o, const RooAbsReal* hTrainTrue, const RooAbsReal* hTrain,
                   const RooAbsReal* hTrue, const RooAbsReal* hMeas, const RooAbsReal* hReco,
                   Int_t _nm, Int_t _nt, Bool_t _overflow,
                   RooUnfolding::ErrorTreatment withError, Double_t chi_squ){
    // TODO
  }
  template<> RooAbsReal* resize<RooAbsReal> (RooAbsReal* h, Int_t nx, Int_t ny, Int_t nz){
    // TOOD
    return h;
  }

  template<> void printHistogram<RooAbsReal>(const RooAbsReal* h){
    // TODO
  }

  template<> void subtract<RooAbsReal>(RooAbsReal* hist, const TVectorD& vec, double fac){
    // TODO
  }

  template<> void h2v<RooAbsReal>  (const RooAbsReal* h, TVectorD& v, bool overflow){}
  template<> void h2ve<RooAbsReal>  (const RooAbsReal* h, TVectorD& v, bool overflow){}    

  template<> TVectorD h2v<RooAbsReal>  (const RooAbsReal* h, bool overflow){
    TVectorD v;
    h2v(h,v);
    return v;    
  }
  template<> TVectorD h2ve<RooAbsReal>  (const RooAbsReal* h, bool overflow){
    TVectorD v;
    h2ve(h,v);
    return v;    
  }
  template<> void h2mNorm<RooAbsReal>  (const RooAbsReal* h, TMatrixD& m, const RooAbsReal* norm, bool overflow){
    // sets Matrix to values of bins in a 2D input histogram    
  }
  template<> void h2meNorm<RooAbsReal>  (const RooAbsReal* h, TMatrixD& m, const RooAbsReal* norm, bool overflow){
    // sets Matrix to errors of bins in a 2D input histogram    
  }
  template<> TMatrixD h2mNorm<RooAbsReal>  (const RooAbsReal* h, const RooAbsReal* norm, bool overflow){
    // Returns Matrix of values of bins in a 2D input histogram
    TMatrixD m;
    h2mNorm(h,m,norm, overflow);
    return m;
  }
  template<> TMatrixD h2meNorm<RooAbsReal>  (const RooAbsReal* h, const RooAbsReal* norm, bool overflow){
    // Returns Matrix of errors of bins in a 2D input histogram
    TMatrixD m;
    h2meNorm(h,m,norm, overflow);
    return m;
  }
  template<> void h2m  (const RooAbsReal* h, TMatrixD& m, bool overflow) { h2mNorm (h,m,(const RooAbsReal*)NULL,overflow); }
  template<> void h2me  (const RooAbsReal* h, TMatrixD& m, bool overflow){ h2meNorm(h,m,(const RooAbsReal*)NULL,overflow); };  
  template<> TMatrixD h2m<RooAbsReal>  (const RooAbsReal* h, bool overflow){
    // Returns Matrix of values of bins in a 2D input histogram
    TMatrixD m;
    h2m(h,m,overflow);
    return m;
  }
  template<> TMatrixD h2me<RooAbsReal>  (const RooAbsReal* h, bool overflow){
    // Returns Matrix of errors of bins in a 2D input histogram
    TMatrixD m;
    h2me(h,m,overflow);
    return m;
  }
}  

  


