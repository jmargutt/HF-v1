#include <TSystem.h>
#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TF1.h>
#include <iostream>
#include <TCanvas.h> 
#include <TVirtualFitter.h>
#include <TMinuit.h>
#include <TGraphAsymmErrors.h>

//#include <AliFlowCommonHist.h>
//#include <AliFlowCommonHistResults.h>
//#include <AliHFMassFitter.h>
#include "FlowFitter.h"

using namespace::std;

void LoadLibraries();
TH1F* GetMassMerged(TFile *f, double &ptmin, double &ptmax, int iBinMin=3, int iBinMax=6);
TH1F* GetMassMergedOld(TFile *rooFile, double &ptmin, double &ptmax, int iBinMin=3, int iBinMax=5, int binWitdh=1, TString meth="QC", TString qvector="", TString gap="", TString refPart="TPC");
TH1F* GetV2Merged(TFile *f, int iBinMin=3, int iBinMax=6);
TH1F* GetV2MergedOld(TFile *rooFile,int iBinMin, int iBinMax, int binWitdh=1, TString meth="QC", TString qvector="", TString gap="",  TString refPart="TPC") ;
TH1F* GetAvgPtOld(TFile *rooFile,int iBinMin, int iBinMax, int binWitdh=1, TString meth="QC", TString qvector="", TString gap="") ;
double ExtractV2Nominal(TH1D *handler, int ptBinMin, int ptBinMax);
double ExtractV2Error(TH1D *handler, int ptBinMin, int ptBinMax);
TH1F* MergeFlowHistograms(TH1F *hist1, TH1F* hist2, char *name);
void FitMassPlotStandard(TH1F* massPlot, Double_t xMinFit, Double_t xMaxFit, int bkgFit);
TH1F* ReBinFlowHistogram(TH1F *hist1, Int_t nnewbins, Double_t *xnewBins);

TString outputFileName = "FinishAnalysisResults_2060MB_pPb_EtaGap.root";
TString finalPlotsFileName="prova2060SP_EtaqGap.root"; //file to store the final plots
 
TString postfixCuts = "Loose";
Int_t typeOfBkgFit = 1; //0 low power - 1  Power function conv. with exponential Fit
Int_t typeOfv2BkgFit =0;
Double_t valMin=0.1396, valMax=0.158; //for DStar
//Double_t valMin=1.71, valMax = 2.06;
//TString method = "SP TPC"; //name for plots
TString method = "SP TPC "; //name for plots
TString methUsed = "SP";
TString rfp = "TPC";
TString Dmeson = "DStar";
const Int_t nptBins=9;

int plotDSnew(Bool_t saveHisto=kTRUE){

  LoadLibraries();
  TFile *file = TFile::Open(outputFileName.Data(),"READ");

  //To be improved... this depend from the width used! 
 // Int_t ptBinMin[nptBins]={3,5,7,9}; 
  //Int_t ptBinMax[nptBins]={4,6,8,12};
  //Int_t ptBinMin[nptBins]={4,5,6,7}; //w=1, 3-4, 4-5, 5-6, 6-8
  //Int_t ptBinMax[nptBins]={4,5,6,8};
 
  Int_t ptBinMin[nptBins]={2,3,4,5,6,7,8,9,10};
  Int_t ptBinMax[nptBins]={2,3,4,5,6,7,8,9,10};
  Int_t width[nptBins]={1,1,1,1,1,1,1,1,1};
  //To be improved... should be calculated as avearge pt
  Double_t ptBin[nptBins]={1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5};
  Double_t ptBinLimit[nptBins+1]={1,2,3,4,5,6,7,8,9,10};
 
  Double_t ptBinErr[nptBins]={0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};

  Double_t ptMin[nptBins], ptMax[nptBins];

  Double_t arrayv2[nptBins], arrayv2Err[nptBins], arraymass[nptBins], arraymassErr[nptBins], arraysigma[nptBins], arraysigmaErr[nptBins];
  Double_t arraySign[nptBins], arraySignErr[nptBins];
  Double_t arrayAvgPt[nptBins], arrayAvgPtErrL[nptBins], arrayAvgPtErrR[nptBins];

  for (Int_t i=0; i<nptBins; i++){
    arrayv2[i]=0;
    arrayv2Err[i]=0;
    arraymass[i]=0;
    arraymassErr[i]=0;
    arraysigma[i]=0;
    arraysigmaErr[i]=0;
    arraySign[i]=0;
    arraySignErr[i]=0;
    arrayAvgPt[i]=0;
    arrayAvgPtErrL[i]=0;
    arrayAvgPtErrR[i]=0;
  }


  for (Int_t i=0; i<nptBins;i++){
    TH1F *m;
    TH1F *flow;
    TH1F *avgPt;

    if (methUsed=="QC"){
      m = GetMassMergedOld(file,ptMin[i],ptMax[i],ptBinMin[i],ptBinMax[i], width[i]);
      flow = GetV2MergedOld(file,ptBinMin[i],ptBinMax[i], width[i]);
      //avgPt = GetAvgPtOld(file,ptBinMin[i],ptBinMax[i], width[i]);
      // TCanvas *c = new TCanvas(); c->cd(); m->Draw();
      //TCanvas *cf = new TCanvas(); cf->cd(); flow->Draw();
      //      TCanvas *cpt = new TCanvas(); cpt->cd(); avgPt->Draw();
      //m = GetMassMerged(file,ptMin[i],ptMax[i],ptBinMin[i],ptBinMax[i]);
      cout<< ptMin[i]<<" "<< ptMax[i]<<endl;
      //flow = GetV2Merged(file,ptBinMin[i],ptBinMax[i]);
      //m->Rebin(2);
    } else if (methUsed=="SP" && rfp=="TPC"){
      TH1F* mass1 = GetMassMergedOld(file,ptMin[i],ptMax[i],ptBinMin[i],ptBinMax[i],width[i],methUsed, "Qa", "GAP1");
      mass1->Add( GetMassMergedOld(file,ptMin[i],ptMax[i],ptBinMin[i],ptBinMax[i],width[i],methUsed, "Qb", "GAP1"));
      m = mass1;
      TH1F* flow1 = GetV2MergedOld(file,ptBinMin[i],ptBinMax[i],width[i],methUsed, "Qa", "GAP1");
      TH1F* flow2 = GetV2MergedOld(file,ptBinMin[i],ptBinMax[i],width[i],methUsed, "Qb", "GAP1");
      flow = MergeFlowHistograms(flow1, flow2, "FlowMerged");
      //m = GetMassMerged(file,ptMin[i],ptMax[i],ptBinMin[i],ptBinMax[i]);
      cout<< ptMin[i]<<" "<< ptMax[i]<<endl;
      //flow = GetV2Merged(file,ptBinMin[i],ptBinMax[i]);
    } else if (methUsed=="SP" && rfp=="VZE"){
      m =GetMassMergedOld(file,ptMin[i],ptMax[i],ptBinMin[i],ptBinMax[i],width[i],methUsed, "Qa", "GAP1", rfp);
      TH1F* flow1 = GetV2MergedOld(file,ptBinMin[i],ptBinMax[i],width[i],methUsed, "Qa", "GAP1", rfp);
      TH1F* flow2 = GetV2MergedOld(file,ptBinMin[i],ptBinMax[i],width[i],methUsed, "Qb", "GAP1", rfp);
      flow = MergeFlowHistograms(flow1, flow2, "FlowMerged");
      flow = flow2;
    }
    //Double_t xbinnew[10] = {0.138, 0.139, 0.141, 0.143, 0.146, 0.148, 0.150, 0.152,  0.154, 0.158};
    
    //TH1F *histo = ReBinFlowHistogram(flow,9,xbinnew);
    //TCanvas *c = new TCanvas(); c->cd(); histo->Draw();
    //flow=histo;

    FlowFitter *f = new FlowFitter(m,flow, valMin, valMax,kFALSE,typeOfBkgFit, typeOfv2BkgFit);
    if (f) cout<< "I am happy"<<endl;
    else cout<<":("<<endl;
    
    //Initial parameters 
    TString parName[8] = {"v_{2}^{sgn}", "M_{v2}^{bgr}", "B_{v2}^{bgr}", "A^{sgn}", "#mu", "#sigma", "M2^{bgr}", "M1^{bgr}"};
    Double_t parIni[8] =  { +0.12,  +0.2,  +0.2,  50,   0.145,  0.005, 10,   0.5}; //v2POI - M_v2 - B_v2 - IntTot - Mean - Sigma - intBkg - slope
    //Double_t parIni[8]=  { +0.12,  +0.2,  +0.2,  50,   0.145,  0.0005, 10,   0.5}; //v2POI - M_v2 - B_v2 - IntTot - Mean - Sigma - intBkg - expcoeff
    Double_t parMin[8] =  { -2.0,  -10.0, -10.0, 1,   0.1445,  0.0001, 1,   -10.}; 
    Double_t parMax[8] =  { +2.0,  +10.0, +10.0, 1000,  0.1465,  0.010000, 7000, 50}; 
    
    f->SetFitInitialParameters(parIni);
    f->SetFitMinParameters(parMin);
    f->SetFitMaxParameters(parMax);
    f->SetFitParNames(parName);
    f->SetRangeFit(valMin,valMax);
    f->SetTypeOfBkgFit(typeOfBkgFit);
    Printf("The counter i is = %d", i);
    f->SetfCounter(i);
    Printf("type of v2 bkg fit is %d ====>", f->GetTypeOfv2BkgFit());
    f->FitMe();
//      f->DrawFit(Form("doubleFit%d",i),i);
      f->DrawFit(Form("doubleFit%d",i));
    Printf(">>>>>>fine del fit");
    //    if (saveHisto)  f->WriteHisto("./",Form("Pt%1.0f_%1.0f_w%d",ptMin[i], ptMax[i],width[i]));

    //retrieve values for the final plots
    arrayv2[i]=f->Getv2();
    arrayv2Err[i]=f->Getv2Error();
    arraymass[i]=f->GetMean();
    arraymassErr[i]=f->GetMeanError();
    arraysigma[i]=f->GetSigma();
    arraysigmaErr[i]=f->GetSigmaError();
    arraySign[i]=f->GetSignificance();
    arraySignErr[i]=f->GetSignificanceError();
    
    if (methUsed=="QC"){
      /*  FlowFitter *fpt = new FlowFitter(m,avgPt, valMin, valMax,kTRUE);
      if (fpt) cout<< "I am happy"<<endl;
      else cout<<":("<<endl;
      
      //Initial parameters 
         TString parNamefpt[7] = {"pt_{sgn}", "pt_{bkg}", "A^{sgn}", "#mu", "#sigma", "M2^{bgr}", "M1^{bgr}"};
      Double_t parInifpt[7] =  { +7,   7,    50, 0.145,  0.0005,   10,  0.5}; //v2POI - M_v2 - B_v2 - IntTot - Mean - Sigma - intBkg - expcoeff
      Double_t parMinfpt[7] =  { +1,   1,      1, 0.1445,  0.0001, 1, -10};
      Double_t parMaxfpt[7] =  { 12,  12,   1000, 0.1464,  0.01, +7000,    50};
      
      fpt->SetIsPtFit(kTRUE);
      fpt->SetFitInitialParameters(parInifpt);
      fpt->SetFitMinParameters(parMinfpt);
      fpt->SetFitMaxParameters(parMaxfpt);
      fpt->SetFitParNames(parNamefpt);
      fpt->SetRangeFit(valMin,valMax);
      fpt->SetTypeOfBkgFit(typeOfBkgFit);
      
      fpt->FitMe();
      fpt->DrawFit(Form("doubleFitPt%d",i));
      arrayAvgPt[i]=fpt->GetAveragePt();
      Printf("Average Pt = %f",arrayAvgPt[i]);
      arrayAvgPtErrL[i]=arrayAvgPt[i]-ptBinLimit[i];
      Printf("Average Pt Err L= %f",arrayAvgPtErrL[i]);
      arrayAvgPtErrR[i]=ptBinLimit[i+1]-arrayAvgPt[i];
      Printf("Average Pt Err R= %f",arrayAvgPtErrL[i]);*/
      arrayAvgPt[i]= ptBin[i];
      arrayAvgPtErrR[i]= ptBinErr[i];
      arrayAvgPtErrL[i]= ptBinErr[i]; 
    }
  
    if (methUsed=="SP") {
      arrayAvgPt[i]= ptBin[i];
      arrayAvgPtErrR[i]= ptBinErr[i];
      arrayAvgPtErrL[i]= ptBinErr[i]; 
    }
  
  }

  TGraphAsymmErrors *myV2 = new TGraphAsymmErrors(nptBins,arrayAvgPt,arrayv2,arrayAvgPtErrL,arrayAvgPtErrR,arrayv2Err,arrayv2Err);
  myV2->SetTitle(Form("%s v_{2} %s",Dmeson.Data(),method.Data()));
  myV2->SetName(Form("%s v_{2} %s",Dmeson.Data(),method.Data()));
  myV2->SetLineColor( kBlue );
  myV2->SetMarkerColor( kBlue );
  myV2->SetMarkerStyle( 20 );
  myV2->SetFillColor(kWhite);
  // myV2->GetYaxis()->SetRangeUser(-0.1,+0.8);
  myV2->GetXaxis()->SetRangeUser(0.0,12.0);
  myV2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  myV2->GetYaxis()->SetTitle("v_{2}");

  TGraphAsymmErrors *myMass = new TGraphAsymmErrors(nptBins,arrayAvgPt,arraymass,arrayAvgPtErrL,arrayAvgPtErrR,arraymassErr,arraymassErr);
  myMass->SetTitle(Form("%s mass %s",Dmeson.Data(), method.Data()));
  myMass->SetName(Form("%s mass %s",Dmeson.Data(), method.Data()));
  myMass->SetLineColor( kBlue );
  myMass->SetMarkerColor( kBlue );
  myMass->SetMarkerStyle( 20 );
  myMass->SetFillColor(kWhite);
  // myMass->GetYaxis()->SetRangeUser(-0.1,+0.8);
  myMass->GetXaxis()->SetRangeUser(0.0,12.0);
  myMass->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  myMass->GetYaxis()->SetTitle("mass");

  TGraphAsymmErrors *mySigma = new TGraphAsymmErrors(nptBins,arrayAvgPt,arraysigma,arrayAvgPtErrL,arrayAvgPtErrR,arraysigmaErr,arraysigmaErr);
  mySigma->SetTitle(Form("%s sigma %s",Dmeson.Data(),method.Data()));
  mySigma->SetName(Form("%s sigma %s",Dmeson.Data(),method.Data()));
  mySigma->SetLineColor( kBlue );
  mySigma->SetMarkerColor( kBlue );
  mySigma->SetMarkerStyle( 20 );
  mySigma->SetFillColor(kWhite);
  // mySigma->GetYaxis()->SetRangeUser(-0.1,+0.8);
  mySigma->GetXaxis()->SetRangeUser(0.0,12.0);
  mySigma->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  mySigma->GetYaxis()->SetTitle("sigma");

  TGraphAsymmErrors *mySignif = new TGraphAsymmErrors(nptBins,arrayAvgPt,arraySign,arrayAvgPtErrL,arrayAvgPtErrR,arraySignErr,arraySignErr);
  mySignif->SetTitle(Form("%s Significance %s",Dmeson.Data(),method.Data()));
  mySignif->SetName(Form("%s Significance %s",Dmeson.Data(),method.Data()));
  mySignif->SetLineColor( kBlue );
  mySignif->SetMarkerColor( kBlue );
  mySignif->SetMarkerStyle( 20 );
  mySignif->SetFillColor(kWhite);
  //mySignif->GetYaxis()->SetRangeUser(-0.1,+0.8);
  mySignif->GetXaxis()->SetRangeUser(0.0,12.0);
  mySignif->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  mySignif->GetYaxis()->SetTitle("Significance");

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  myV2->Draw("APS");

  TCanvas *c2 = new TCanvas("c2","c2");
  c2->cd();
  myMass->Draw("APS");

  TCanvas *c3 = new TCanvas("c3","c3");
  c3->cd();
  mySigma->Draw("APS");

  TCanvas *c4 = new TCanvas("c4","c4");
  c4->cd();
  mySignif->Draw("APS");

  if (saveHisto){
    cout << "----> Saving output in " << finalPlotsFileName.Data() <<endl;
    TFile *output = new TFile(finalPlotsFileName.Data(),"recreate");
    output->cd();
    myV2->Write();
    myMass->Write();
    mySigma->Write();
    mySignif->Write();
    output->Close();
  }

  printf("===========results \n");
  for (Int_t i=0; i<4;i++){
    printf("%f %f %f %f %f %f %f %f\n", arrayv2[i], arrayv2Err[i], arraymass[i], arraymassErr[i], arraysigma[i], arraysigmaErr[i], arraySign[i], arraySignErr[i]); 
  }
  return 0;
}


// * G E T   M A S S   M E R G E D **************************************************
TH1F* GetMassMerged(TFile *f, double &ptmin, double &ptmax, int iBinMin, int iBinMax) {
  TList* fileKeys = dynamic_cast<TList*>(f->GetListOfKeys());
  TList* listTemp = (TList*)f->Get(fileKeys->At(0)->GetName());
  TH2F *handler = (TH2F*)listTemp->FindObject("Control_Flow_Mass_POIAliFlowCommonHistQC");
  if (!listTemp->FindObject("Control_Flow_Mass_POIAliFlowCommonHistQC")) cout<<"problem"<<endl;
  int bins = handler->GetXaxis()->GetNbins();
  double binMin = handler->GetXaxis()->GetBinLowEdge(1);
  double binMax = handler->GetXaxis()->GetBinLowEdge(bins+1);
   ptmin = handler->GetYaxis()->GetBinLowEdge(iBinMin); //incluso
   ptmax = handler->GetYaxis()->GetBinLowEdge(iBinMax+1); //escluso
  double masswidth = handler->GetXaxis()->GetBinWidth(1)*1000;
  TH1F *mass = new TH1F(Form("Mass_%d%d", iBinMin, iBinMax),Form("(%.0f < p_{t} <%.0f) GeV/c;Mass (GeV/c^{2});Counts/%.1f MeV/c^{2}",
			ptmin,ptmax,masswidth),
			bins,binMin,binMax);
  for(int mb=0; mb!=13; ++mb) {   
      listTemp = (TList*)f->Get(fileKeys->At(mb)->GetName());
      handler = (TH2F*)listTemp->FindObject("Control_Flow_Mass_POIAliFlowCommonHistQC");
      mass->Add( (TH1D*) handler->ProjectionX(Form("MB%d",mb),iBinMin,iBinMax) );
  }
  return mass;
}

 // * G E T   M A S S   M E R G E D **************************************************
TH1F* GetMassMergedOld(TFile *rooFile, double &ptmin, double &ptmax,int iBinMin, int iBinMax, int binWidth, TString meth, TString qvector, TString gap, TString refPart) {
    TString objName;
    TString prefix;
    TString posfix;
    if (meth.EqualTo("QC")){
      objName = Form("AliFlowCommonHistQC");

      //prefix = Form("FlowD2H_QC_v2_Dstar3050w%d%s/DStarw%dQCTPCMB", binWidth, postfixCuts.Data(), binWidth);
        prefix = Form("FlowD2H_QC_v2_EtaGap05_DStar_v2_pPb/DStarw%dQCTPCMB", binWidth);

      posfix=Form("v2_EtaGap05_DStar_v2_pPb");
    } else   if (meth.EqualTo("SP")){
      objName = Form("AliFlowCommonHist_%s",meth.Data());

      //prefix = Form("FlowD2H_%s_v2_Dstar3050w%d%s/DStarw%d%s%sMB",meth.Data(),binWidth, postfixCuts.Data(), binWidth,meth.Data(),refPart.Data());
      prefix = Form("FlowD2H_SP_v2_EtaGap05_DStar_v2_pPb/DStarw%dSPTPCMB", binWidth);

      posfix=Form("SPv2%s%s_EtaGap05_DStar_v2_pPb",qvector.Data(),gap.Data());
    }

    /* if ( method.find("SP")!=std::string::npos) { // THE SP CASE
      objName = Form("AliFlowCommonHist_%s",  method.data() );
      prefix = Form("FlowD2H_%s_v%d_Dstar3050QCSPw%dLoose/%sw%d%s%sMB",
			       method.data(),
			       harmonic,
		               ptwidth,
			       dmeson.data(),
		               ptwidth,
			       subMethod.data(),
			       rfpSource.data());

      posfix=Form("%sv%d%s%s_Dstar3050QCSPw%dLoose", method.data(), harmonic, Qt, gap, ptwidth);
      }*/
    
    cout<<prefix.Data()<<endl;

    cout<<posfix.Data()<<endl;

    TList *listCC = (TList*)  rooFile->Get(Form("%s0%s", prefix.Data(), posfix.Data()));
    if (!listCC) printf("ERROR\n");
    AliFlowCommonHist *common = (AliFlowCommonHist*) listCC->FindObject( objName.Data() );
    if (!common) { printf("ERROR reading FILE\n"); return NULL;}
    TH2F *handler = common->GetHistMassPOI();
    int bins = handler->GetXaxis()->GetNbins();
    double binMin = handler->GetXaxis()->GetBinLowEdge(1);
    double binMax = handler->GetXaxis()->GetBinLowEdge(bins+1);
    ptmin = handler->GetYaxis()->GetBinLowEdge(iBinMin); //incluso
    ptmax = handler->GetYaxis()->GetBinLowEdge(iBinMax+1); //escluso
    double masswidth = handler->GetXaxis()->GetBinWidth(1)*1000;
    TH1F *mass = new TH1F(Form("MASS_%d%d",iBinMin,iBinMax),
			  Form(" (%.0f < p_{t} <%.0f GeV/c); M(k#pi#pi)-M(k#pi) (GeV/c^{2});Counts/%.1f MeV/c^{2}",
			       ptmin,ptmax,masswidth),
			  bins,binMin,binMax);
    //for(int mb=0; mb!= numberOfMB; ++mb) {
    for(int mb=0; mb!=13; ++mb) {   
      listCC = (TList*)  rooFile->Get(Form("%s%d%s",
					   prefix.Data(),
					   mb,
					   posfix.Data()));
      if(!listCC) { printf("ERROR in list %d\n",mb); continue; }
      common = (AliFlowCommonHist*) listCC->FindObject( objName.Data() );
      handler = (TH2F*) common->GetHistMassPOI();
      mass->Add( (TH1F*) handler->ProjectionX(Form("MB%d",mb),iBinMin,iBinMax) );
    }
    cout << "DONE!" <<endl;
    return mass;
  }

TH1F *GetAvgPtOld(TFile *rooFile,int iBinMin, int iBinMax, int binWidth, TString meth, TString qvector, TString gap){
  TString objReference;
  TString objName;
  TString prefix;
  TString posfix;

  if (meth.EqualTo("QC")){
    objName = Form("AliFlowCommonHistQC");
    
    prefix = Form("FlowD2H_QC_v2_Dstar3050w%d%s/DStarw%dQCTPCMB", binWidth, postfixCuts.Data(), binWidth);
      
    posfix=Form("v2_Dstar3050w%d%s", binWidth,postfixCuts.Data());
  } else   if (meth.EqualTo("SP")){
    objName = Form("AliFlowCommonHist_%s",meth.Data());
    
    prefix = Form("FlowD2H_%s_v2_Dstar3050w%d%s/DStarw%d%sTPCMB",meth.Data(),binWidth, postfixCuts.Data(), binWidth,meth.Data());
    
    posfix=Form("%sv2%s%s_Dstar3050w%d%s",meth.Data(),qvector.Data(),gap.Data(),binWidth,postfixCuts.Data());
  } 
  
  cout<<prefix.Data()<<endl;
  
  cout<<posfix.Data()<<endl;
 
  TList *listCC = (TList*)  rooFile->Get(Form("%s0%s", prefix.Data(), posfix.Data()));
  if (!listCC) printf("ERROR\n");
  AliFlowCommonHist *common = (AliFlowCommonHist*) listCC->FindObject( objName.Data() );
  if (!common) { printf("ERROR reading FILE\n"); return NULL;}
  TProfile *handler = common->GetHistProMeanPtperBin();
  if (!handler) {Printf("handler not exists"); return NULL;}

  Double_t binLimits[14] = {0.138, 0.139, 0.141, 0.143, 0.144, 0.145, 0.146, 0.147, 0.148, 0.150, 0.152,  0.154, 0.156, 0.158};

TH1F *avgpt = new TH1F(Form("AVGPT_%d%d",iBinMin,iBinMax),
			  Form(" (p_{t}); M(k#pi#pi)- M(k#pi) (GeV/c^{2});<pT>"),
			  13,binLimits);

 for(int mb=0; mb!=13; ++mb) {
      listCC = (TList*) rooFile->Get(Form("%s%d%s",prefix.Data(),mb,posfix.Data()));
      if(!listCC) { printf("ERROR in list %d\n",mb); continue; }
      common = (AliFlowCommonHist*) listCC->FindObject( objName.Data() );
      handler = (TProfile*) common->GetHistProMeanPtperBin();
      double avg=0, err=0;
      int ent=0;
      for(int j=iBinMin; j!=iBinMax+1; ++j) {
	
	avg+=handler->GetBinContent(j)*handler->GetBinEntries(j);
	err+=handler->GetBinError(j)*handler->GetBinError(j);
	ent+=handler->GetBinEntries(j);
	Printf("mb =%d, j=%d, Bin content = %f, Entries = %f, Error = %f",mb,j, handler->GetBinContent(j), handler->GetBinEntries(j), handler->GetBinError(j));
      }
      if (ent!=0) avg /= ent;
      else avg=0;
      if (err!=0) err = sqrt(err);
      else err=0;
      avgpt->SetBinContent( mb+1, avg );
      avgpt->SetBinError( mb+1, err );
 }
    /*   TF1 *myAvg = new TF1("myAvgPt","[0]",0.140,0.158);
    mass->Fit(myAvg);
    avgPt = myAvg->GetParameter(0);
    eLowPt = avgPt - ptmin;
    eHigPt = ptmax - avgPt;*/
 return avgpt;
}


// * G E T  V2   M E R G E D **************************************************
TH1F* GetV2Merged(TFile *f, int iBinMin, int iBinMax) {
  // const Int_t nBins= numberOfMB;
  Double_t binLimits[14] = {0.138, 0.139, 0.141, 0.143, 0.144, 0.145, 0.146, 0.147, 0.148, 0.150, 0.152,  0.154, 0.156, 0.158};
  TH1D *handler;
  //Getting mass bins
    TList* fileKeys = dynamic_cast<TList*>(f->GetListOfKeys());
    TList* listTemp = (TList*)f->Get(fileKeys->At(0)->GetName());
    TH2F *ref2 = (TH2F*)listTemp->FindObject("Control_Flow_Mass_POIAliFlowCommonHistQC");
    double ptmin = ref2->GetYaxis()->GetBinLowEdge(iBinMin); 
    double ptmax = ref2->GetYaxis()->GetBinLowEdge(iBinMax+1);
    TH1F*merged = new TH1F(Form("V2_%d%d",iBinMin,iBinMax), Form("(%.0f < p_{t} <%.0f GeV/c); M(k#pi#pi)-M(k#pi)(GeV/c^{2});v_{2}", 
								   ptmin,ptmax), 13, binLimits);

    for(int mb=0; mb!=13; ++mb) {
      listTemp = (TList*)f->Get(fileKeys->At(mb)->GetName());
      handler = (TH1D*)listTemp->FindObject(Form("diffFlow_2p_FlowD2H_QC_v2_Dstar3050QCSPw4Loo_DStarw4QCTPCMB%dv2_Dstar3050QCSPw4Loo",mb));
      if (!handler) {continue;};
      double v2nom = ExtractV2Nominal(handler,iBinMin,iBinMax);
      double v2err = ExtractV2Error(handler,iBinMin,iBinMax);   
      merged->SetBinContent(mb+1,v2nom);
      merged->SetBinError(mb+1,v2err);
    }
  
    return merged;
  }

// * G E T   V 2   M E R G E D ********************************************
TH1F* GetV2MergedOld(TFile *rooFile,int iBinMin, int iBinMax, int binWidth, TString meth, TString qvector, TString gap, TString refPart) {
    TString objReference;
    TString objName;
    TString prefix;
    TString posfix;

if (meth.EqualTo("QC")){
      objName = Form("AliFlowCommonHistResults2ndOrderQC");

      objReference = Form("AliFlowCommonHistQC");      
        prefix = Form("FlowD2H_QC_v2_EtaGap05_DStar_v2_pPb/DStarw%dQCTPCMB", binWidth);
    
        posfix=Form("v2_EtaGap05_DStar_v2_pPb");

 } else if (meth.EqualTo("SP")){
      objName = Form("AliFlowCommonHistResults_SP");

      objReference = Form("AliFlowCommonHist_SP");      
      
 //     prefix = Form("FlowD2H_%s_v2_Dstar3050w%d%s/DStarw%d%s%sMB", meth.Data(), binWidth, postfixCuts.Data(), binWidth, meth.Data(),refPart.Data());

  //    posfix=Form("%sv2%s%s_Dstar3050w%d%s",meth.Data(), qvector.Data(), gap.Data(), binWidth, postfixCuts.Data());
  
     prefix = Form("FlowD2H_SP_v2_EtaGap05_DStar_v2_pPb/DStarw%dSPTPCMB", binWidth);
     
     posfix=Form("SPv2%s%s_EtaGap05_DStar_v2_pPb",qvector.Data(),gap.Data());
     
 }

    cout<<prefix.Data()<<endl;

    cout<<posfix.Data()<<endl;

    //    const Int_t nBins= numberOfMB;
    Double_t binLimits[14] = {0.138, 0.139, 0.141, 0.143, 0.144, 0.145, 0.146, 0.147, 0.148, 0.150, 0.152,  0.154, 0.156, 0.158};
    AliFlowCommonHistResults *common;
    TH1D *handler;
    //Getting mass bins
    TList *listCC = (TList*)  rooFile->Get( Form("%s0%s",prefix.Data(),posfix.Data()));
    AliFlowCommonHist *ref1 = (AliFlowCommonHist*) listCC->FindObject( objReference.Data() );
    if (!ref1) { printf("ERROR reading FILE\n"); return NULL;}
    TH2F *ref2 = (TH2F*) ref1->GetHistMassPOI();
    //int bins = ref2->GetXaxis()->GetNbins();
    //double binMin = ref2->GetXaxis()->GetBinLowEdge(1);
    //double binMax = ref2->GetXaxis()->GetBinLowEdge(bins+1);    
    double ptmin = ref2->GetYaxis()->GetBinLowEdge(iBinMin); 
    double ptmax = ref2->GetYaxis()->GetBinLowEdge(iBinMax+1);
    // TH1D *merged = new TH1D(Form("V2_%d%d",iBinMin,iBinMax),
    //Form("%s %s (%.0f,%.0f) GeV/c;Mass (GeV/c^{2});v_{2}",
    //				  dmeson.data(),
    //				  ccName.data(),
    //				 ptmin,ptmax),
    //                           bins,binMin,binMax);
    TH1F *merged = new TH1F(Form("V2_%d%d",iBinMin,iBinMax), Form(" (%.0f < p_{t} <%.0f GeV/c); M(k#pi#pi)-M(k#pi)(GeV/c^{2});v_{2}", 
								   ptmin,ptmax), 13, binLimits);

    // for(int mb=0; mb!= numberOfMB; ++mb) {
    for(int mb=0; mb!=13; ++mb) {
      listCC = (TList*)  rooFile->Get( Form("%s%d%s",prefix.Data(),mb,posfix.Data()));
      common = (AliFlowCommonHistResults*) listCC->FindObject( objName.Data() );
      handler = (TH1D*) common->GetHistDiffFlowPtPOI();
      double v2nom = ExtractV2Nominal(handler,iBinMin,iBinMax);
      cout<< "mb " << mb << "v2 = "<< v2nom <<endl;
      double v2err = ExtractV2Error(handler,iBinMin,iBinMax);   
      merged->SetBinContent(mb+1,v2nom);
      merged->SetBinError(mb+1,v2err);
    }
  cout << "DONE!" <<endl;
    return merged;
  }

// * E X T R A C T   V 2   N O M I N A L ***********************************
double ExtractV2Nominal(TH1D *handler, int ptBinMin, int ptBinMax) {
  double dSum1=0, dSum2=0;
  for(int i=ptBinMin; i!=ptBinMax+1; ++i) {
    double value = handler->GetBinContent(i);
    double error = handler->GetBinError(i);
    if(error>0.) {
      dSum1+=value/(error*error);
      dSum2+=1./(error*error);
    }
  }
  if(dSum2>0) return dSum1/dSum2;
  return 0;
}
// * E X T R A C T   V 2   E R R O R ****************************************
double ExtractV2Error(TH1D *handler, int ptBinMin, int ptBinMax) {
  double dSum2=0;
  for(int i=ptBinMin; i!=ptBinMax+1; ++i) {
    double error = handler->GetBinError(i);
    if(error>0.)
      dSum2+=1./(error*error);
  }
  if(dSum2>0) return pow(1./dSum2,0.5);
  return 0;
}
// * M E R G E  V2  P L O T ****************************************
 TH1F* MergeFlowHistograms(TH1F *hist1, TH1F* hist2, char *name) {
    TH1F *ret = (TH1F*) hist1->Clone( name );
    for( int i=1; i!=ret->GetNbinsX()+1; ++i ) {
      double dSum1=0, dSum2=0;
      double value = hist1->GetBinContent(i);
      double error = hist1->GetBinError(i);
      if(error>0.) {
	dSum1+=value/(error*error);
	dSum2+=1./(error*error);
      }
      value = hist2->GetBinContent(i);
      error = hist2->GetBinError(i);
      if(error>0.) {
	dSum1+=value/(error*error);
	dSum2+=1./(error*error);
      }
      if (dSum2>0){
	ret->SetBinContent(i, dSum1/dSum2 );
	ret->SetBinError(i, pow(1./dSum2,0.5) );
      }
    }
    return ret;
  }
  

/*
*/
/*double Flow(double x, double *p, double xmin=-1, double xmax=-1) {
  double fracSgn = Sgn(x,p);
  double fracBgr = Bgr(x,p);
  double total = fracSgn+fracBgr;
  fracSgn/=total;
  fracBgr/=total;
  //  printf( " %f |%f %f | %f %f ===>  ",total, fracSgn,fracBgr, fracSgn/(fracSgn+fracBgr), fracBgr/(fracSgn+fracBgr) );
  if((xmin+xmax)>0){
    TF1 sgnf("sgnfrac",SgnFracPointer,xMinFit,xMaxFit,8); sgnf.SetParameters(p);
    TF1 bgrf("bgrfrac",BgrFracPointer,xMinFit,xMaxFit,8); bgrf.SetParameters(p);
    fracSgn = sgnf.Integral(xmin,xmax)/(xmax-xmin);
    fracBgr = bgrf.Integral(xmin,xmax)/(xmax-xmin);
    }*/
  /*  double total = fracSgn+fracBgr;
    if(total<1e-3) return 0.0;
    fracSgn/=total;
    fracBgr/=total;*/ //old version
    //Printf( "%f |%f %f | %f %f | %f %f ", fracSgn+fracBgr, fracSgn,fracBgr, fracSgn/(fracSgn+fracBgr), fracBgr/(fracSgn+fracBgr), xmax, xmin );
/*  return double( fracSgn*p[0]+fracBgr*FlowBgr(x,p) );
    }*/
/*double FlowPointer(double *x, double *p) {
  return double( Flow(x[0],p) );
}
double FlowBgrPointer(double *x, double *p) {
  return double( FlowBgr(x[0],p) );
}
void myFcn(Int_t& , Double_t* , Double_t &val, Double_t *par, Int_t ) {
  double tmp;
  int ndf=0;
  double xMinFit=0.140;
  double xMaxFit=0.158;
  int nPars=8;

  int minBinFit, maxBinFit;
  // YIELD
  minBinFit = mass->FindBin(xMinFit+1e-6);
  maxBinFit = mass->FindBin(xMaxFit-1e-6);
  double lastChi2Yield=0;
  for(int mb=minBinFit; mb!=maxBinFit; ++mb) {
    if(mass->GetBinCenter(mb+1)<mpi)
      continue;
    if (mass->GetBinContent(mb+1)==0) {Printf("Bin Content =0"); continue;}
    tmp =mass->GetBinContent(mb+1);
    tmp-=Yield(mass->GetBinCenter(mb+1),par);
    tmp/=mass->GetBinError(mb+1);
    lastChi2Yield+=tmp*tmp;
    ++ndf;
  }
  // FLOW
  minBinFit = flow->FindBin(xMinFit+1e-6);
  maxBinFit = flow->FindBin(xMaxFit-1e-6);
  //Printf("minBINFIT = %d, maxBINFIT = %d", minBinFit, maxBinFit);
  double lastChi2Flow=0;
  for(int mb=minBinFit; mb!=maxBinFit; ++mb) {
    if(flow->GetBinCenter(mb+1)<mpi)
      continue;
    if(flow->GetBinContent(mb+1)==0) {Printf("Bin Content =0"); continue;}
    tmp =flow->GetBinContent(mb+1);
    tmp-=Flow(flow->GetBinCenter(mb+1),par,
	      flow->GetBinLowEdge(mb+1),
	      flow->GetBinLowEdge(mb+2));
    tmp/=flow->GetBinError(mb+1);
    lastChi2Flow += tmp*tmp;
    ++ndf;
  }
  double lastChi2 = lastChi2Yield+lastChi2Flow;
  val = lastChi2;
}

  void FitMe() {
    double parIni[8] =  { +0.12,  +0.2,  +0.2,  5,   0.145,  0.0006, 3300,   4};
    double parMin[8] =  { -2.0,  -10.0, -10.0, 1,   0.1445,  0.0001, 1,   -10.}; 
    double parMax[8] =  { +2.0,  +10.0, +10.0, 1000,  0.1465,  0.01000, 7000, 50}; 

    TVirtualFitter::SetDefaultFitter("Minuit");
    TVirtualFitter *minuit = TVirtualFitter::Fitter(NULL,nPars); //???
    minuit->SetParameter(0,"v2sgn",     parIni[0], 1e-6, parMin[0], parMax[0]);
    minuit->SetParameter(1,"M_{v2bgr}", parIni[1], 1e-6, parMin[1], parMax[1]);
    minuit->SetParameter(2,"B_{v2bgr}", parIni[2], 1e-6, parMin[2], parMax[2]);
    minuit->SetParameter(3,"A_{Sgn}",   parIni[3], 1e-2, parMin[3], parMax[3]);
    minuit->SetParameter(4,"meanMass",  parIni[4], 1e-6, parMin[4], parMax[4]);
    minuit->SetParameter(5,"sigmaMass", parIni[5], 1e-6, parMin[5], parMax[5]);
    minuit->SetParameter(6,"M_{Bgr}",   parIni[6], 1e-2, parMin[6], parMax[6]);
    minuit->SetParameter(7,"B_{Bgr}",   parIni[7], 1e-2, parMin[7], parMax[7]);
    minuit->SetFCN(myFcn);
    double argList[100];
    argList[0] = 5000; // FUNCTION CALLS
    argList[1] = 1e-6; // TOLERANCE
    minuit->ExecuteCommand("MIGRAD",argList,2);
    for(int i=0; i!=8; ++i) {
      parValues[i] = minuit->GetParameter(i);
      parErrors[i] = minuit->GetParError(i);
    }
    TF1 *massFit = new TF1("massFit", YieldPointer, xMinFit, xMaxFit, 8);
    massFit->SetParameters(parValues);
    massFit->SetLineColor(kBlue);
    mass->GetListOfFunctions()->Add(massFit);
    TF1 *bgrFunct = new TF1("bgrFunct", BgrPointer, xMinFit, xMaxFit, 8);
    bgrFunct->SetParameters(parValues);
    mass->GetListOfFunctions()->Add(bgrFunct);
    TF1 *flowFit = new TF1("flowFit", FlowPointer, xMinFit, xMaxFit, 8);
    flowFit->SetParameters(parValues);
    flow->GetListOfFunctions()->Add(flowFit);   
  }
*/

void LoadLibraries(){
    gSystem->Load("libGeom");
    gSystem->Load("libVMC");
    gSystem->Load("libXMLIO");
    gSystem->Load("libPhysics");
  
    //load needed libraries:
    gSystem->AddIncludePath("-I$ROOTSYS/include");
    //gSystem->Load("libTree");

    // for AliRoot
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libPWGflowBase");
}

void FitMassPlotStandard(TH1F* massPlot, Double_t xMinFit, Double_t xMaxFit, int bkgFit){
  TCanvas *canvas2 = new TCanvas("c","c");
    AliHFMassFitter *fitter;
    fitter = new AliHFMassFitter(massPlot, massPlot->GetBinLowEdge(3), massPlot->GetBinLowEdge(massPlot->GetNbinsX()-2),1,bkgFit);
    fitter->SetRangeFit(xMinFit, xMaxFit);
    fitter->SetInitialGaussianSigma(0.01);
    fitter->SetInitialGaussianMean(1.86);
    Bool_t ok=fitter->MassFitter(kFALSE);
    if (ok) fitter->DrawHere(canvas2->cd());
    canvas2->DrawClone();
    //  massFromFit[ptBin]=fitter->GetMean();
    // sigmaFromFit[ptBin]=fitter->GetSigma();
    //  fitter->Signal(3,signal[ptBin],esignal[ptBin]);
    //  fitter->Background(3,bkg[ptBin],ebkg[ptBin]);
    // fitter->Significance(3,significance[ptBin],esignificance[ptBin]);
    }

TH1F* ReBinFlowHistogram(TH1F *src, Int_t nbins, Double_t *bin ) {
    TH1F *ret = new TH1F( Form("%s_rebin",src->GetName()),
			  Form("%s_rebin",src->GetTitle()),
			  nbins, bin );
    for( int j=1; j!=nbins+1; ++j ) {
      double dSum1=0, dSum2=0;
      printf("%d in [%f,%f] (",j,bin[j-1],bin[j]);
      for(int i=1; i!=src->GetNbinsX()+1; ++i) {
	if( ( src->GetBinLowEdge(i)+1e-5 > bin[j-1] ) &&
	    ( src->GetBinLowEdge(i+1)-1e-5 < bin[j] ) ) {
	  printf("%d,",i);
	  double value = src->GetBinContent(i);
	  //	  if(fabs(value)>0.5) continue; /// TO BE REMOVED: TO EXCLUDE OUTLIERS FOUND IN SP
	  double error = src->GetBinError(i);
	  if(error>0.) {
	    dSum1+=value/(error*error);
	    dSum2+=1./(error*error);
	  }
	}
      }
      printf("): %f +- %f \n",dSum1/dSum2,pow(1./dSum2,0.5));
      if ( dSum1 !=0 ) ret->SetBinContent(j, dSum1/dSum2 ); else ret->SetBinContent(j, 0);
      if ( dSum2 != 0 ) ret->SetBinError(j, pow(1./dSum2,0.5) ); else ret->SetBinError(j, 0 );
    }
    return ret;
  }
