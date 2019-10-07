/*! \class RooUnfoldSvdT
Links to TSVDUnfold class which unfolds using Singular Value Decomposition (SVD).
Regularisation parameter defines the level at which values are deemed to be due to statistical fluctuations and are cut out. (Default= number of bins/2)
Returns errors as a full matrix of covariances
Can only handle 1 dimensional distributions
Can account for both smearing and biasing
*/

#include "RooUnfoldSvd.h"

#include <iostream>
#include <iomanip>

#include "TClass.h"
#include "TNamed.h"
#include "TBuffer.h"
#include "TH1.h"
#include "TH2.h"
#include "TVectorD.h"
#include "TMatrixD.h"

#include "RooUnfoldHelpers.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldTH1Helpers.h"
#include "RooUnfoldFitHelpers.h"

using namespace RooUnfolding;

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::RooUnfoldSvdT (const RooUnfoldSvdT<Hist,Hist2D>& rhs)
  : RooUnfoldT<Hist,Hist2D> (rhs)
{
  //! Copy constructor.
  Init();
  CopyData (rhs);
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::RooUnfoldSvdT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t kreg,
                            const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D> (res, meas, name, title), _kreg(kreg ? kreg : res->GetNbinsTruth()/2)
{
  //! Constructor with response matrix object and measured unfolding input histogram.
  //! The regularisation parameter is kreg.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::RooUnfoldSvdT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t kreg, Int_t ntoyssvd,
                            const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D> (res, meas, name, title), _kreg(kreg ? kreg : res->GetNbinsTruth()/2)
{
  //! Constructor with old ntoyssvd argument. No longer required.
  Init();
  this->_NToys = ntoyssvd;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::Reset()
{
  // destroy and re-initialize this object
  Destroy();
  Init();
  RooUnfoldT<Hist,Hist2D>::Reset();
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::Destroy()
{
  //! delete all members of this object
  delete this->_svd;
  delete this->_meas1d;
  delete this->_train1d;
  delete this->_truth1d;
  delete this->_reshist;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::Init()
{
  //! initialize this object with zero values
  this->_svd= 0;
  this->_meas1d= this->_train1d= this->_truth1d= 0;
  this->_reshist= this->_meascov= 0;
  GetSettings();
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::Assign (const RooUnfoldSvdT<Hist,Hist2D>& rhs)
{
  //! assign data from another instance
  RooUnfoldT<Hist,Hist2D>::Assign (rhs);
  CopyData (rhs);
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::CopyData (const RooUnfoldSvdT<Hist,Hist2D>& rhs)
{
  //! copy data from another instance
  this->_kreg= rhs._kreg;
}

template<class Hist,class Hist2D> typename RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold*
RooUnfoldSvdT<Hist,Hist2D>::Impl()
{
  //! retrieve the SVDUnfold object
  return this->_svd;
}


template<class Hist,class Hist2D> void RooUnfoldSvdT<Hist,Hist2D>::PrepareHistograms() const{
  //! fill all the histogram members
  this->_meas1d = this->_meas;
  this->_train1d= this->_res->Hmeasured();
  this->_truth1d= this->_res->Htruth();
  this->_reshist= this->_res->Hresponse();
}

namespace{
  TH1* histNoOverflow(const TH1* hist, bool overflow){
    return createHist<TH1>(h2v(hist,overflow),h2ve(hist,overflow),name(hist),title(hist),vars(hist),overflow);
  }
  void subtract(TH1* hist, const TVectorD& vec, double fac){
    for (Int_t i= 1; i<=hist->GetNbinsX()+1; i++){
      hist->SetBinContent (i, hist->GetBinContent(i)-(fac*vec[i-1]));
    }
  }
}
template<> void RooUnfoldSvdT<TH1,TH2>::PrepareHistograms() const {
  //! fill all the histogram members
  TH1* meas1d = ::histNoOverflow (this->_meas,             this->_overflow);
  TH1* train1d= ::histNoOverflow (this->_res->Hmeasured(), this->_overflow);
  TH1* truth1d= ::histNoOverflow (this->_res->Htruth(),    this->_overflow);
  TH2* reshist= this->_res->HresponseNoOverflow();
  RooUnfolding::resize (meas1d,  this->_nb);
  RooUnfolding::resize (train1d, this->_nb);
  RooUnfolding::resize (truth1d, this->_nb);
  RooUnfolding::resize (reshist, this->_nb, this->_nb);

  // Subtract fakes from measured distribution
  if (this->_res->HasFakes()) {
    TVectorD fakes(this->_res->Vfakes());
    Double_t fac= this->_res->Vmeasured().Sum();
    if (fac!=0.0) fac=  this->Vmeasured().Sum() / fac;
    if (this->_verbose>=1) std::cout << "Subtract " << fac*fakes.Sum() << " fakes from measured distribution" << std::endl;
    ::subtract(meas1d,fakes,fac);
  }

  delete this->_meas1d  ;
  delete this->_train1d ;
  delete this->_truth1d ;
  delete this->_reshist ;
  
  this->_meas1d  = meas1d  ;
  this->_train1d = train1d;
  this->_truth1d = truth1d;
  this->_reshist = reshist ;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::Unfold() const
{
  //! perform the unfolding
  if (this->_res->GetDimensionTruth() != 1 || this->_res->GetDimensionMeasured() != 1) {
    std::cerr << "RooUnfoldSvdT may not work very well for multi-dimensional distributions" << std::endl;
  }
  if (this->_kreg < 0) {
    std::cerr << "RooUnfoldSvdT invalid kreg: " << this->_kreg << std::endl;
    return;
  }

  this->_nb= this->_nm > this->_nt ? this->_nm : this->_nt;
 
  if (this->_kreg > this->_nb) {
    std::cerr << "RooUnfoldSvdT invalid kreg=" << this->_kreg << " with " << this->_nb << " bins" << std::endl;
    return;
  }

  this->PrepareHistograms();
  
  const TMatrixD& cov(this->GetMeasuredCov());

  if (this->_verbose>=1) std::cout << "SVD init " << nBins(this->_reshist,X) << " x " << nBins(this->_reshist,Y) << " bins, kreg=" << this->_kreg << std::endl;
  if(!this->_meas1d) throw std::runtime_error("no meas1d given!");
  if(!this->_train1d) throw std::runtime_error("no train1d given!");
  if(!this->_truth1d) throw std::runtime_error("no truth1d given!");
  if(!this->_reshist) throw std::runtime_error("no reshist given!");

  this->_svd= new SVDUnfold (this->_meas1d, cov, this->_train1d, this->_truth1d, this->_reshist);

  this->_cache._rec.ResizeTo (this->_nt);
  this->_cache._rec = this->_svd->UnfoldV (this->_kreg);

  if (this->_verbose>=2) {
    printTable (std::cout, h2v(this->_truth1d), h2v(this->_train1d), h2v(this->_meas1d), this->_cache._rec);
    TMatrixD resmat(h2m(this->_reshist));
    printMatrix(resmat,"SVDUnfold response matrix");
  }

  this->_cache._unfolded= true;
  this->_cache._haveCov=  false;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::GetCov() const
{
  //! Get covariance matrix
  if (!this->_svd) return;
  this->_cache._cov.ResizeTo(this->_nb,this->_nb);
  //Get the covariance matrix for statistical uncertainties on the measured distribution
  if (this->_dosys!=2) this->_cache._cov = this->_svd->GetXtau();
  //Get the covariance matrix for statistical uncertainties on the response matrix
  //Uses Poisson or Gaussian-distributed toys, depending on response matrix histogram's Sumw2 setting.

  if (this->_dosys){
    add(this->_cache._cov,this->_svd->GetAdetCovMatrix (this->_NToys));
  }
  this->_cache._haveCov= true;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::GetWgt() const
{
  //! Get weight matrix
  if (this->_dosys) RooUnfoldT<Hist,Hist2D>::GetWgt();   // can't add sys errors to weight, so calculate weight from covariance
  if (!this->_svd) return;

  //Get the covariance matrix for statistical uncertainties on the measured distribution
  this->_cache._wgt = this->_svd->GetXinv();
  
  this->_cache._haveWgt= true;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::GetSettings() const {
    this->_cache._minparm=0;
    this->_cache._maxparm= this->_meas ? nBins(this->_meas,X) : 0;
    this->_cache._stepsizeparm=1;
    this->_cache._defaultparm=this->_cache._maxparm/2;
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::RooUnfoldSvdT()
  : RooUnfoldT<Hist,Hist2D>()
{
  //! Default constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::RooUnfoldSvdT (const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D>(name,title)
{
  //! Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::RooUnfoldSvdT (const TString& name, const TString& title)
  : RooUnfoldT<Hist,Hist2D>(name,title)
{
  //! Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>& RooUnfoldSvdT<Hist,Hist2D>::operator= (const RooUnfoldSvdT<Hist,Hist2D>& rhs)
{
  //! Assignment operator for copying RooUnfoldSvdT settings.
  Assign(rhs);
  return *this;
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::~RooUnfoldSvdT()
{
  Destroy();
}


template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::SetKterm (Int_t kreg)
{
  //! Set regularisation parameter
  this->_kreg= kreg;
}


template<class Hist,class Hist2D> Int_t
RooUnfoldSvdT<Hist,Hist2D>::GetKterm() const
{
  //! Return regularisation parameter
  return this->_kreg;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::SetRegParm (Double_t parm)
{
  //! Set regularisation parameter
  SetKterm(Int_t(parm+0.5));
}

template<class Hist,class Hist2D> Double_t
RooUnfoldSvdT<Hist,Hist2D>::GetRegParm() const
{
  //! Return regularisation parameter
  return GetKterm();
}

template class RooUnfoldSvdT<TH1,TH2>;
ClassImp (RooUnfoldSvd)

#ifndef NOROOFIT
typedef RooUnfoldSvdT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist> RooFitUnfoldSvd;
template class RooUnfoldSvdT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>;
ClassImp (RooFitUnfoldSvd)
#endif

