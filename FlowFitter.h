#ifndef FLOWFITTER_H
#define FLOWFITTER_H

#include <TNamed.h>
#include <TMath.h>

#include <TF1.h>
class TH1F;

Double_t mpi = 0.13957;

class FlowFitter : public TNamed {

 public:
  FlowFitter();
  FlowFitter(const TH1F* massHistoToFit, const TH1F* flowHistoToFit, Double_t minvalue, Double_t maxvalue, Bool_t isPtFit=kFALSE, Int_t typeofbkg=2, Int_t typeofv2bkg=0);
  FlowFitter(const FlowFitter &mfit);
  virtual ~FlowFitter();
  
  //setters
  void     SetHisto(const TH1F *massHistoToFit, const TH1F *flowHistoToFit);
  void     SetRangeFit(Double_t minvalue, Double_t maxvalue){fminMass=minvalue; fmaxMass=maxvalue;}
  void     SetTypeOfBkgFit(Int_t typeBkgFit){ftypeofbkg=typeBkgFit;}
  void     SetTypeOfv2BkgFit(Int_t typev2BkgFit){ftypeofv2bkg=typev2BkgFit;}
  void     SetInitialGaussianMean(Double_t mean) {fMass=mean;} // change the default value of the mean
  void     SetInitialGaussianSigma(Double_t sigma) {fSigmaSgn=sigma;} // change the default value of the sigma
  Bool_t   SetFixThisParam(Int_t thispar,Bool_t fixpar);
  void     SetFitInitialParameters(Double_t *initialPars) {fFitParsInitial=initialPars;}
  void     SetFitMinParameters(Double_t *minPars) {fFitParsMin=minPars;}
  void     SetFitMaxParameters(Double_t *maxPars) {fFitParsMax=maxPars;}
  void     SetFitParNames(TString *namesPars) {fparNames=namesPars;}
  void     SetfCounter(Int_t counter) {fCounter=counter;}
  void     SetIsPtFit(Bool_t isPtFit) {fIsPtFit=isPtFit;}

  //getters
  TH1F*    GetHistoClone(Int_t type) const; //return the histogram
  void     GetRangeFit(Double_t &minvalue, Double_t &maxvalue) const {minvalue=fminMass; maxvalue=fmaxMass;}
  Double_t GetMinRangeFit()const {return fminMass;}
  Double_t GetMaxRangeFit()const {return fmaxMass;}
  Int_t    GetTypeOfBkgFit() const {return ftypeofbkg;}
  Int_t    GetTypeOfv2BkgFit() const {return ftypeofv2bkg;}
  Int_t    GetBinN()       const {return fNbin;}
  Int_t    GetNFinalPars() const {return fNFinalPars;}
  Double_t* GetFitInitialParameters() const {return fFitParsInitial;}
  Double_t* GetFitMinParameters() const {return fFitParsMin;}
  Double_t* GetFitMaxParameters() const {return fFitParsMax;}
  Double_t GetMean() const {return fMass;}
  Double_t GetMeanError() const {return fMassError;}
  Double_t GetSigma()const {return fSigmaSgn;}
  Double_t GetSigmaError()const {return fSigmaError;}
  Double_t GetSignal() const {return fSignal;}
  Double_t GetSignalError() const {return fEsgn2;}
  Double_t GetSignificance() const {return fSignificance;}
  Double_t GetSignificanceError() const {return fESignificance;}
  Double_t GetAveragePt() const {return fAvgPt;}
  Double_t GetAveragePtError() const {return fAvgPtError;}
  Double_t Getv2() const {return fv2;}
  Double_t Getv2Error() const {return fv2Error;}
  Double_t GetChiSquareYield() const {return fLastChi2Yield;}
  Double_t GetChiSquareFlow() const {return fLastChi2Flow;}
  Double_t GetChiSquare() const {return fLastChi2;}
  Bool_t*  GetFixParam()const {return fFixPar;}
  Bool_t   GetFixThisParam(Int_t thispar)const;

  void     AddFunctionsToHisto();
  void     DrawFit(TString cname) const;

  void ComputeSignificance();
  void WriteHisto (TString path="./", TString suffix="") const;
    //  void WriteCanvas(TString path="./", TString suffix="", TCanvas *c) const;

  static     FlowFitter* GetInstance() {return fgInstance;}

 private:

  void     ComputeNFinalPars();
// * F I T   P A R T ****************************************
  double Sgn(double x, double *p) {   
    if(x<=mpi) return 0.0;
    if (fIsPtFit) return double( p[2]*TMath::Gaus(x,p[3],p[4]));
    else return double( p[3]*TMath::Gaus(x,p[4],p[5]));
  }
  double Bgr(double x, double *p) {
    if(x<=mpi) return 0.0;
    if (fIsPtFit) {
      if (ftypeofbkg==2) return double( p[5]*exp(p[6]*(x-p[3])));
      else if (ftypeofbkg==0) return double( p[5]*TMath::Power((x-mpi),p[6]));
      else return double( p[5]*sqrt(x-mpi)*exp(-1*p[6]*(x-mpi)) );
    }
    else { if (ftypeofbkg==0) return double( p[6]*TMath::Power((x-mpi),p[7]));
      else if (ftypeofbkg==1) return double( p[6]*sqrt(x-mpi)*exp(-1*p[7]*(x-mpi)) );
      else if (ftypeofbkg==2) return double( p[6]*exp(p[7]*(x-p[4])));
      else return double (p[6]*(x-p[4])+p[7]);
    }
    //return  p[6]*(2*TMath::Power(x-p[4],2)-1)+p[7]*(x-p[4])+p[8]; //too complex function??
    //return  p[6]*p[7]*(x-p[4])+p[8]; //too complex function??
  }


    double FlowBgr(double x, double *p) {
      if(x<=mpi) return 0.0;
      if (fIsPtFit) return double (p[1]);
      else {
	//if (ftypeofv2bkg==0) return double( p[1]*(x-p[4])+p[2] );
	if (ftypeofv2bkg==0) return double( p[1]*(x-mpi)+p[2] );
	//else return double( p[1]*(x-p[4])+p[2]+p[8]*(x-p[4])*(x-p[4]));
	else return double( p[1]*(x-mpi)+p[2]+p[8]*(x-mpi)*(x-mpi));
      }
    }

  double Yield(double x, double *p) {
    return double( Sgn(x,p)+Bgr(x,p) );
  }
  double YieldPointer(double *x, double *p) {
    return double( Yield(x[0],p) );
  }
  double SgnPointer(double *x, double *p) {
    return double( Sgn(x[0],p) );
  }
  double BgrPointer(double *x, double *p) {
    return double( Bgr(x[0],p) );
  }

  Double_t SgnFracPointer(Double_t *x, Double_t *p) {
    return Double_t( Sgn(x[0],p)/(Sgn(x[0],p)+Bgr(x[0],p)) );
  }
  double BgrFracPointer(double *x, double *p) {
    return double( Bgr(x[0],p)/(Sgn(x[0],p)+Bgr(x[0],p)) );
  }

  double Flow(double x, double *p, double xmin=-1, double xmax=-1) {
    double fracSgn = Sgn(x,p);
    double fracBgr = Bgr(x,p);
    double total = fracSgn+fracBgr;
    fracSgn/=total;
    fracBgr/=total;
    //  printf( " %f |%f %f | %f %f ===>  ",total, fracSgn,fracBgr, fracSgn/(fracSgn+fracBgr), fracBgr/(fracSgn+fracBgr) );
    if((xmin+xmax)>0){
      TF1 sgnf("sgnfrac",this,&FlowFitter::SgnFracPointer,fminMass,fmaxMass,8); sgnf.SetParameters(p);
      TF1 bgrf("bgrfrac",this,&FlowFitter::BgrFracPointer,fminMass,fmaxMass,8); bgrf.SetParameters(p);
      fracSgn = sgnf.Integral(xmin,xmax)/(xmax-xmin);
      fracBgr = bgrf.Integral(xmin,xmax)/(xmax-xmin);
    }
    /*  double total = fracSgn+fracBgr;
	if(total<1e-3) return 0.0;
	fracSgn/=total;
	fracBgr/=total;*/ //old version
    //Printf( "%f |%f %f | %f %f | %f %f ", fracSgn+fracBgr, fracSgn,fracBgr, fracSgn/(fracSgn+fracBgr), fracBgr/(fracSgn+fracBgr), xmax, xmin );
    return double( fracSgn*p[0]+fracBgr*FlowBgr(x,p) );
  }


  double FlowPointer(double *x, double *p) {
    return double( Flow(x[0],p) );
  }
  double FlowBgrPointer(double *x, double *p) {
    return double( FlowBgr(x[0],p) );
  }

  static void myFcn(Int_t& , Double_t* , Double_t &val, Double_t *par, Int_t );
  void FitMe();
    
  TH1F*     fhistoInvMass;     // histogram to fit
  TH1F*     fhistoFlow;     // histogram to fit
  Int_t ftypeofbkg;      //type of background for yield
  Int_t ftypeofv2bkg;    //0=linear, 1=quadratic
  Double_t  fminMass;          // lower mass limit
  Double_t  fmaxMass;          // upper mass limit
  Int_t fNbin; 
  Int_t fNFinalPars;          
  Double_t*  fFitParsInitial;      //[fNFinalPars] Initial param
  Double_t*  fFitParsMin;          //[fNFinalPars] Min param
  Double_t*  fFitParsMax;          //[fNFinalPars] Max param
  Double_t* fparValues;            //[fNFinalPars] 
  Double_t* fparErrors;            //[fNFinalPars] 
  TString* fparNames;             //[fNFinalPars] 
  Double_t fLastChi2Yield;
  Double_t fLastChi2Flow;
  Double_t fLastChi2;
  
  Double_t fMass;
  Double_t fMassError;
  Double_t fSigmaSgn;
  Double_t fSigmaError;
  Double_t fv2;
  Double_t fv2Error;
  Bool_t*   fFixPar;           //[fNFinalPars] for each par if kTRUE it is fixed in fit

  Double_t  fSignal;
  Double_t  fBackground;
  Double_t  fSignificance;
  Double_t  fEsgn2;
  Double_t  fEbgr2;
  Double_t  fESignificance;

  Double_t  fAvgPt;
  Double_t  fAvgPtError;

  Int_t fCounter;
  Bool_t fIsPtFit;
  static FlowFitter* fgInstance;         // global pointer on itself
  Double_t fChi2NDF;
  ClassDef(FlowFitter,1); // class for invariant mass fit
};

#endif
