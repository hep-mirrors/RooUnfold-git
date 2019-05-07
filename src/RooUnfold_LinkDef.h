#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class RooUnfold-;
#pragma link C++ class RooFitUnfold-;
#pragma link C++ class RooUnfoldBayes+;
#pragma link C++ class RooUnfoldSvd-;
#pragma link C++ class RooFitUnfoldSvd-;
#pragma link C++ class RooUnfoldSvd::SVDUnfold-;
#pragma link C++ class RooFitUnfoldSvd::SVDUnfold-;
#pragma link C++ class RooUnfoldBinByBin+;
#pragma link C++ class RooFitUnfoldBinByBin+;
#pragma link C++ class RooUnfoldResponse-;
#pragma link C++ class RooFitUnfoldResponse-;
#pragma link C++ class RooUnfoldErrors+;
#pragma link C++ class RooUnfoldParms+;
#pragma link C++ class RooUnfoldInvert+;
#ifndef NOTUNFOLD
#pragma link C++ class RooUnfoldTUnfold+;
#endif
#ifdef HAVE_DAGOSTINI
#pragma link C++ class RooUnfoldDagostini+;
#endif
#pragma link C++ class RooUnfoldIds-;

#endif
