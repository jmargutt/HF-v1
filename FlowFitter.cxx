#include <TH1F.h>
#include <TMath.h>
#include <TVirtualFitter.h>
#include <TList.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TVirtualFitter.h>
#include <TPaveText.h>
#include <TFile.h>
#include <iostream>

#include "FlowFitter.h"

using std::cout;
using std::endl;

FlowFitter* FlowFitter::fgInstance = 0; 

ClassImp(FlowFitter)

 //************** constructors
FlowFitter::FlowFitter(): 
  TNamed(),
  fhistoInvMass(0),
  fhistoFlow(0),
  ftypeofbkg(0),
  ftypeofv2bkg(0),
  fminMass(0),
  fmaxMass(0),
  fNbin(0),
  fNFinalPars(0),
  fFitParsInitial(0),
  fFitParsMin(0),
  fFitParsMax(0),
  fparValues(0),
  fparErrors(0),
  fparNames(0),
  fLastChi2Yield(0),
  fLastChi2Flow(0),
  fLastChi2(0),
  fMass(0),
  fMassError(0),
  fSigmaSgn(0),
  fSigmaError(0),
  fFixPar(0),
  fSignal(0),
  fBackground(0),
  fSignificance(0),
  fEsgn2(0),
  fEbgr2(0),
  fESignificance(0),
  fAvgPt(0),
  fAvgPtError(0),
  fCounter(0),
  fIsPtFit(kFALSE),
  fChi2NDF(0)
{
  // default constructor
  cout<<"Default constructor:"<<endl;
  cout<<"Remember to set the Histo, the Type, the FixPar"<<endl;
  fgInstance = this;
}

//___________________________________________________________________________
FlowFitter::FlowFitter (const TH1F *massHistoToFit, const TH1F *flowHistoToFit,  Double_t minvalue, Double_t maxvalue, Bool_t isPtFit, Int_t typeofbkg, Int_t typeofv2bkg): 
  TNamed(),
  fhistoInvMass(0),
  fhistoFlow(0),
  ftypeofbkg(0),
  ftypeofv2bkg(0),
  fminMass(0),
  fmaxMass(0),
  fNbin(0),
  fNFinalPars(1),
  fFitParsInitial(0),
  fFitParsMin(0),
  fFitParsMax(0),
  fparValues(0),
  fparErrors(0), 
  fparNames(0), 
  fLastChi2Yield(0),
  fLastChi2Flow(0),
  fLastChi2(0),
  fMass(1.865),
  fMassError(0),
  fSigmaSgn(0.012),
  fSigmaError(0),
  fFixPar(0),
  fSignal(0),
  fBackground(0),
  fSignificance(0),
  fEsgn2(0),
  fEbgr2(0),
  fESignificance(0),
  fAvgPt(0),
  fAvgPtError(0),
  fCounter(0),
  fIsPtFit(kFALSE),
  fChi2NDF(0)
{
  fhistoInvMass= (TH1F*)massHistoToFit->Clone("fhistoInvMass");
  fhistoFlow= (TH1F*)flowHistoToFit->Clone("fhistoFlow");	
  fhistoInvMass->SetDirectory(0);
  fminMass=minvalue;	
  fmaxMass=maxvalue;	
  fIsPtFit=isPtFit;
  ftypeofv2bkg=typeofv2bkg;
  ftypeofbkg=typeofbkg;
  ComputeNFinalPars();
  fFixPar=new Bool_t[fNFinalPars];
  for(Int_t i=0;i<fNFinalPars;i++){
    fFixPar[i]=kFALSE;
  }
  fFitParsInitial = new Double_t[fNFinalPars];
  fFitParsMin = new Double_t[fNFinalPars];
  fFitParsMax = new Double_t[fNFinalPars];
  fparValues = new Double_t[fNFinalPars];
  fparErrors = new Double_t[fNFinalPars];
  fparNames = new TString[fNFinalPars];
  for(Int_t i=0;i<fNFinalPars;i++){
    fFitParsInitial[i]=0.;
    fFitParsMin[i]=0.;
    fFitParsMax[i]=0.;
    fparValues[i]=0.;
    fparErrors[i]=0.;
    fparNames[i]="";
  }
  fgInstance = this;
}

FlowFitter::FlowFitter(const FlowFitter &mfit):
  TNamed(),
  fhistoInvMass(mfit.fhistoInvMass),
  fhistoFlow(mfit.fhistoFlow),
  ftypeofbkg(mfit.ftypeofbkg),
ftypeofv2bkg(mfit.ftypeofv2bkg),
  fminMass(mfit.fminMass),
  fmaxMass(mfit.fmaxMass),
  fNbin(mfit.fNbin),
  fNFinalPars(mfit.fNFinalPars),
  fFitParsInitial(0),
  fFitParsMin(0),
  fFitParsMax(0),
  fparValues(0),
  fparErrors(0),
  fparNames(mfit.fparNames),
  fLastChi2Yield(mfit.fLastChi2Yield),
  fLastChi2Flow(mfit.fLastChi2Flow),
  fLastChi2(mfit.fLastChi2),
  fMass(mfit.fMass),
  fMassError(mfit.fMassError),
  fSigmaSgn(mfit.fSigmaSgn),
  fSigmaError(mfit.fSigmaError),
  fv2(mfit.fv2),
  fv2Error(mfit.fv2Error),
  fFixPar(0),
  fSignal(mfit.fSignal),
  fBackground(mfit.fBackground),
  fSignificance(mfit.fSignificance),
  fEsgn2(mfit.fEsgn2),
  fEbgr2(mfit.fEbgr2),
  fESignificance(mfit.fESignificance),
  fAvgPt(mfit.fAvgPt),
  fAvgPtError(mfit.fAvgPtError),
  fCounter(mfit.fCounter),
  fIsPtFit(mfit.fIsPtFit),
  fChi2NDF(mfit.fChi2NDF)
{
  //copy constructor
  
  //if(mfit.fParsSize > 0){
    //fFitPars=new Float_t[fParsSize];
    fFixPar=new Bool_t[fNFinalPars];
    //memcpy(fFitPars,mfit.fFitPars,mfit.fParsSize*sizeof(Float_t));
    memcpy(fFixPar,mfit.fFixPar,mfit.fNFinalPars*sizeof(Bool_t));
    fFitParsInitial= new Double_t[fNFinalPars];
    fFitParsMin= new Double_t[fNFinalPars];
    fFitParsMax= new Double_t[fNFinalPars];
    fparValues = new Double_t[fNFinalPars];
    fparErrors = new Double_t[fNFinalPars];
    fparNames = new TString[fNFinalPars];
    memcpy(fFitParsInitial,mfit.fFitParsInitial,mfit.fNFinalPars*sizeof(Double_t));
    memcpy(fFitParsMin,mfit.fFitParsMin,mfit.fNFinalPars*sizeof(Double_t));
    memcpy(fFitParsMax,mfit.fFitParsMax,mfit.fNFinalPars*sizeof(Double_t));
    memcpy(fparValues,mfit.fparValues,mfit.fNFinalPars*sizeof(Double_t));
    memcpy(fparErrors,mfit.fparErrors,mfit.fNFinalPars*sizeof(Double_t));
    memcpy(fparNames,mfit.fparNames,mfit.fNFinalPars*sizeof(TString));
    //}
    fgInstance = this;
}

//_________________________________________________________________________

FlowFitter::~FlowFitter() {

  //destructor

  cout<<"FlowFitter destructor called"<<endl;

  delete fhistoInvMass;
  delete fhistoFlow;

  delete[] fFitParsInitial;
  delete[] fFitParsMin;
  delete[] fFitParsMax;
  delete[] fparValues;
  delete[] fparErrors;
  delete[] fparNames;

  delete[] fFixPar;
  fgInstance = 0;
}

//___________________________________________________________________________
void FlowFitter::SetHisto(const TH1F *massHistoToFit,const TH1F *flowHistoToFit ){

  fhistoInvMass = new TH1F(*massHistoToFit);
  fhistoInvMass->SetDirectory(0);
  fhistoFlow = new TH1F(*flowHistoToFit);
  fhistoFlow->SetDirectory(0);
  //cout<<"SetHisto pointer "<<fhistoInvMass<<endl;
}

//___________________________________________________________________________
Bool_t FlowFitter::SetFixThisParam(Int_t thispar,Bool_t fixpar){

  //set the value (kFALSE or kTRUE) of one element of fFixPar
  //return kFALSE if something wrong

  if(thispar>=fNFinalPars) {
    cout<<"Error! Parameter out of bounds! Max is "<<fNFinalPars-1<<endl;
    return kFALSE;
    }
  /*if(!fFixPar){
    cout<<"Initializing fFixPar...";
    SetDefaultFixParam();
    cout<<" done."<<endl;
    }*/

  fFixPar[thispar]=fixpar;
  if(fixpar)cout<<"Parameter "<<thispar<<" is now fixed"<<endl;
  else cout<<"Parameter "<<thispar<<" is now free"<<endl;
  return kTRUE;
}

//___________________________________________________________________________
void FlowFitter::ComputeNFinalPars() {

  //compute the number of parameters of the total (signal+bgk) function
  cout<<"Info:ComputeNFinalPars... ";
  if (fIsPtFit) fNFinalPars=2; 
  else {
    fNFinalPars =1;   //v2 signal
    switch (ftypeofv2bkg){
    case 0: 
      fNFinalPars+=2; //flow bkg linear parametrization
      break;
    case 1:
      fNFinalPars+=3; //flow bkg quadratic parametrization
      break;
    default:
      cout<<"Error in computing fNFinalPars: check ftypeOfFit4Bkg"<<endl;
      break;
    }
  }

  switch ( ftypeofbkg) {//npar background func
    case 0: //0 low power - 1  Power function conv. with exponential Fit
    fNFinalPars+=2;
    break;
  case 1:
    fNFinalPars+=2;
    break;
  case 2:
    fNFinalPars+=2;
    break;
  case 3:
    fNFinalPars+=2;
    break;
  default:
    cout<<"Error in computing fNFinalPars: check ftypeOfFit4Bkg"<<endl;
    break;
    }

  

  fNFinalPars+=3; //gaussian signal
  cout<<": "<<fNFinalPars<<endl;
}

//_________________________________________________________________________

TH1F* FlowFitter::GetHistoClone(Int_t typeHisto) const{
  TH1F* hout;
  switch (typeHisto) {
  case 0:
    hout=(TH1F*)fhistoInvMass->Clone(fhistoInvMass->GetName());
    break;
  case 1:
    hout=(TH1F*)fhistoFlow->Clone(fhistoFlow->GetName());
    break;
  default:
    cout<<"ERROR: wrong type of hist requested" <<endl;
    break;
   }
  return hout;
}

//___________________________________________________________________________
Bool_t FlowFitter::GetFixThisParam(Int_t thispar)const{
  //return the value of fFixPar[thispar]
  if(thispar>=fNFinalPars) {
    cout<<"Error! Parameter out of bounds! Max is "<<fNFinalPars-1<<endl;
    return kFALSE;
  }
  if(!fFixPar) {
    cout<<"Error! Parameters to be fixed still not set"<<endl;
    return kFALSE;
  }
  return fFixPar[thispar];
} 

//___________________________________________________________________________
void FlowFitter::myFcn(Int_t& , Double_t* , Double_t &val, Double_t *par, Int_t ) {
    double tmp;
    int ndf=0;
    
    int minBinFit, maxBinFit;
    Double_t lastChi2Yield=0, lastChi2Flow=0, lastChi2=0;
    // YIELD
    FlowFitter *f = FlowFitter::GetInstance();
    TH1F* histoInvMass = f->GetHistoClone(0);
    minBinFit = histoInvMass->FindBin(f->GetMinRangeFit()+1e-6);
    maxBinFit = histoInvMass->FindBin(f->GetMaxRangeFit()-1e-6);
    for(int mb=minBinFit; mb!=maxBinFit; ++mb) {
      if(histoInvMass->GetBinCenter(mb+1)<mpi)
     	continue;
      if (histoInvMass->GetBinContent(mb+1)==0) {Printf("Bin Content Mass =0"); continue;}
      tmp =histoInvMass->GetBinContent(mb+1);
      tmp-=f->Yield(histoInvMass->GetBinCenter(mb+1),par);
      tmp/=histoInvMass->GetBinError(mb+1);
      lastChi2Yield+=tmp*tmp;
      ++ndf;
    }
    // // FLOW
    TH1F* histoFlow = f->GetHistoClone(1);
    minBinFit = histoFlow->FindBin(f->GetMinRangeFit()+1e-6);
    maxBinFit = histoFlow->FindBin(f->GetMaxRangeFit()-1e-6);
    // //Printf("minBINFIT = %d, maxBINFIT = %d", minBinFit, maxBinFit);
    
    for(int mb=minBinFit; mb!=maxBinFit; ++mb) {
      if(histoFlow->GetBinCenter(mb+1)<mpi)
	continue;
      if(histoFlow->GetBinContent(mb+1)==0) {Printf("Bin Content Flow =0"); continue;}
      tmp =histoFlow->GetBinContent(mb+1);
      tmp-=f->Flow(histoFlow->GetBinCenter(mb+1),par,
		   histoFlow->GetBinLowEdge(mb+1),
		   histoFlow->GetBinLowEdge(mb+2));
      tmp/=histoFlow->GetBinError(mb+1);
      lastChi2Flow += tmp*tmp;
	++ndf;
    }
    lastChi2 = lastChi2Yield+lastChi2Flow;
    val = lastChi2;
}

//___________________________________________________________________________
void FlowFitter::FitMe() {

    TVirtualFitter::SetDefaultFitter("Minuit");
    TVirtualFitter *minuit = TVirtualFitter::Fitter(NULL,fNFinalPars); 
    for (Int_t i=0; i<fNFinalPars; i++) minuit->SetParameter(i,fparNames[i],fFitParsInitial[i], 1e-6, fFitParsMin[i], fFitParsMax[i]);
    /* minuit->SetParameter(0,"v2sgn",     fFitParsInitial[0], 1e-6, fFitParsMin[0], fFitParsMax[0]);
    minuit->SetParameter(1,"M_{v2bgr}", fFitParsInitial[1], 1e-6, fFitParsMin[1], fFitParsMax[1]);
    minuit->SetParameter(2,"B_{v2bgr}", fFitParsInitial[2], 1e-6, fFitParsMin[2], fFitParsMax[2]);
    minuit->SetParameter(3,"A_{Sgn}",   fFitParsInitial[3], 1e-2, fFitParsMin[3], fFitParsMax[3]);
    minuit->SetParameter(4,"meanMass",  fFitParsInitial[4], 1e-6, fFitParsMin[4], fFitParsMax[4]);
    minuit->SetParameter(5,"sigmaMass", fFitParsInitial[5], 1e-6, fFitParsMin[5], fFitParsMax[5]);
    minuit->SetParameter(6,"Int_{bkg}", fFitParsInitial[6], 1e-2, fFitParsMin[6], fFitParsMax[6]);
    minuit->SetParameter(7,"Slope",     fFitParsInitial[7], 1e-2, fFitParsMin[7], fFitParsMax[7]);*/

    minuit->SetFCN(FlowFitter::myFcn);
    double argList[100];
    argList[0] = 5000; // FUNCTION CALLS
    argList[1] = 1e-6; // TOLERANCE
    minuit->ExecuteCommand("MIGRAD",argList,2);
    for(int i=0; i!=fNFinalPars; ++i) {
      fparValues[i] = minuit->GetParameter(i);
      fparErrors[i] = minuit->GetParError(i);
    }

    Printf("-------------> %f",  fChi2NDF); 
    if (!fIsPtFit) {
      fMass= fparValues[4];
      fMassError = fparErrors[4];
      fSigmaSgn= fparValues[5];
      fSigmaError = fparErrors[5];
      fv2= fparValues[0];
      fv2Error= fparErrors[0];
      
      AddFunctionsToHisto();
      ComputeSignificance();
    } else {
      fAvgPt=fparValues[0];
      fAvgPtError=fparErrors[0];
      AddFunctionsToHisto();
    }      
//DrawFit();
}
//___________________________________________________________________________
void FlowFitter::AddFunctionsToHisto(){
  TF1 *massFit = new TF1("massFit", this,&FlowFitter::YieldPointer, fminMass, fmaxMass, fNFinalPars);
  massFit->SetParameters(fparValues);
  massFit->SetLineColor(kBlue);

  TList *listFunct = fhistoInvMass->GetListOfFunctions();
  TList *listFunctFlow = fhistoFlow->GetListOfFunctions();

  listFunct->Add(massFit);
  TF1 *bgrFunct = new TF1("bgrFunct", this, &FlowFitter::BgrPointer, fminMass, fmaxMass, fNFinalPars);
  bgrFunct->SetParameters(fparValues);
  listFunct->Add(bgrFunct);

  TF1 *flowFit = new TF1("flowFit", this, &FlowFitter::FlowPointer, fminMass, fmaxMass, fNFinalPars);
  flowFit->SetParameters(fparValues);
  listFunctFlow->Add(flowFit);   
 
}

//_________________________________________________________________________
void FlowFitter::DrawFit(TString cvname) const{
  //TCanvas FlowFitter::DrawFit(Int_t nsigma) const{
  //draws histogram together with fit functions with default nice colors
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBottomMargin(0.15);
 
  TString cvtitle="fit of ";
  cvtitle+=fhistoInvMass->GetName();
  cvtitle+=fCounter;
  // TString cvname="cc";
  //cvname+=fCounter;
  
 TCanvas* c= new TCanvas(cvname,cvtitle);
 Printf(">>>>>>>>>>>Creating canvas with name %s", c->GetName());
 c->Divide(1,2,0,0);
 c->cd(1); gPad->SetRightMargin(0.05);
 TH1F* hdraw=GetHistoClone(0);
  if(!hdraw->GetFunction("massFit")){
    cout<<"Probably fit failed and you didn't try to refit with background only, there's no function to be drawn"<<endl;
    return;
  }
 
  hdraw->SetMinimum(0);
  hdraw->GetXaxis()->SetRangeUser(1.71,2.055);
  hdraw->SetMarkerStyle(20);
  Printf(">>>>>>>>>>>Drawing histo with name %s in canvas %s", hdraw->GetName(), c->GetName());
  hdraw->Draw("PE");
  //if(hdraw->GetFunction("massFit")) hdraw->GetFunction("massFit")->DrawClone("same");

  /*TPaveText *pinfo=new TPaveText(0.6,0.86,1.,1.,"NDC");
  pinfo->SetBorderSize(0);
  pinfo->SetFillStyle(0);
  TString strI = "Pb-Pb, #sqrt{s_{NN}} = 2.76 TeV";
  pinfo->AddText(strI);
  strI = "Centrality 30-50%";
  pinfo->AddText(strI);
  strI = "7 * 10^{6} events";
  pinfo->AddText(strI);
  strI = "D*^{+} #rightarrow D^{0} #pi^{+}";
  pinfo->AddText(strI);
  strI = "4 < p_{T} < 6 GeV/c";
  pinfo->AddText(strI);
  pinfo->DrawClone();*/

  if (!fIsPtFit) {
    TPaveText *pinfoM=new TPaveText(0.6,0.7,1.,1.,"NDC");
    pinfoM->SetBorderSize(0);
    pinfoM->SetFillStyle(0);
    TString str=Form("Significance: %.1f #pm %.1f", fSignificance, fESignificance);
    pinfoM->AddText(str);
    str = Form("#mu = %.2f #pm %.2f",fMass*1000,fMassError*1000);
    pinfoM->AddText(str);
    str = Form("#sigma = %.3f #pm %.3f",fSigmaSgn*1000, fSigmaError*1000);
    pinfoM->AddText(str);
    str = Form("S/B = %.3f",fSignal/fBackground);
    pinfoM->AddText(str);  
    pinfoM->Draw();
  }
  c->cd(2); ;gPad->SetRightMargin(0.05); 
  gPad->SetTopMargin(0); 
  hdraw=GetHistoClone(1); hdraw->SetTitle(""); hdraw->GetYaxis()->SetRangeUser(-0.22,0.98); hdraw->Draw("E");  
  hdraw->SetMarkerStyle(20);
  hdraw->GetXaxis()->SetRangeUser(1.71,2.05);
  hdraw->Draw("PE");

  if (!fIsPtFit) {
    TPaveText *pinfoF=new TPaveText(0.6,0.8,1.,.87,"NDC");
    pinfoF->SetBorderSize(0);
    pinfoF->SetFillStyle(0);
    TString strF = Form("v_{2}^{sgn} = %.3f #pm %.3f",fv2,fv2Error);
    pinfoF->AddText(strF);
    pinfoF->Draw();
  }
}

//_________________________________________________________________________

void FlowFitter::ComputeSignificance() {
  TF1 *SgnFit = new TF1("SgnFit", this, &FlowFitter::SgnPointer, fminMass, fmaxMass, fNFinalPars);
  SgnFit->SetParameters(fparValues);
  TF1 *BgrFit = new TF1("BgrFit", this, &FlowFitter::BgrPointer, fminMass, fmaxMass, fNFinalPars);
  BgrFit->SetParameters(fparValues);
    double xmin = fparValues[4]-3*fparValues[5];
    double xmax = fparValues[4]+3*fparValues[5];
    int imin = fhistoInvMass->FindBin(xmin);
    int imax = fhistoInvMass->FindBin(xmax);
    double xfmin = fhistoInvMass->GetBinLowEdge(imin);
    double xfmax = fhistoInvMass->GetBinLowEdge(imax+1);
    fSignal = SgnFit->Integral(xfmin,xfmax)/fhistoInvMass->GetBinWidth(1);
    fBackground = BgrFit->Integral(xfmin,xfmax)/fhistoInvMass->GetBinWidth(1);
    double sgnInt=fparValues[3]*fparValues[5]*sqrt(2*TMath::Pi());
    double esgnInt=sqrt(2*TMath::Pi()*(fparErrors[5]*fparErrors[5]*fparValues[3]*fparValues[3]+ fparErrors[3]*fparErrors[3]*fparValues[5]*fparValues[5]));
    fEsgn2=(esgnInt/sgnInt)*fSignal;
    double xMinSideBand=fparValues[4]+4*fparValues[5];
    int binMinSideBand=fhistoInvMass->GetBinLowEdge(xMinSideBand);
    int binMaxSideBand=fhistoInvMass->GetNbinsX();
    double bgrSide = fhistoInvMass->Integral(binMinSideBand,binMaxSideBand);
    double sum2=0;
    for (int i=binMinSideBand; i<=binMaxSideBand;i++){
      sum2+=fhistoInvMass->GetBinError(i)*fhistoInvMass->GetBinError(i);
    }
    double ebgrSide=sqrt(sum2);
    fEbgr2=ebgrSide/bgrSide*fBackground;
    //double ebgr = sqrt(fBackground); //to be corrected
    fSignificance = fSignal/sqrt(fSignal+fBackground);
    //old calc
    //  esignificance = sqrt(+significance*significance/(signal+background)/(signal+background)*(1/4.*esgn*esgn+ebgr*ebgr)
    //			 +background/(signal+background)/(signal+background)*esgn*esgn);
    //new calc
    fESignificance = fSignificance * sqrt((fEsgn2*fEsgn2+fEbgr2*fEbgr2)/(4.*(fSignal+fBackground)*(fSignal+fBackground))+(fBackground/(fSignal+fBackground))*fEsgn2*fEsgn2/fSignal/fSignal);
    // Printf("=======SIGNIFICANCE");
    //Printf("%.3f %.3f", significance, esignificance);
}

//_________________________________________________________________________

void FlowFitter::WriteHisto(TString path, TString suffix) const {

  //Write the histogram in the output file 

  TH1F* hgetM=(TH1F*)fhistoInvMass->Clone();
  TH1F* hgetF=(TH1F*)fhistoFlow->Clone();
 
  path+= "FlowFitterOutput";
  path+= suffix.Data();
  path += ".root";
  TFile *output;
 
  output = new TFile(path.Data(),"recreate");
  output->cd();
  hgetM->Write();
  hgetF->Write();
  cout<<" "<<hgetM->GetName()<<" written in "<<path<<endl;
  cout<<" "<<hgetF->GetName()<<" written in "<<path<<endl;
  delete output;
  
}
