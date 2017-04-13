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

#include <THnSparse.h>
//#include <AliFlowCommonHist.h>
//#include <AliFlowCommonHistResults.h>
//#include <AliHFMassFitter.h>
#include "FlowFitter.h"

using namespace::std;

void LoadLibraries();


TString outputFileName = "AnalysisResults.root";
TString finalPlotsFileName="plot.root"; //file to store the final plots

TString postfixCuts = "Loose";
Int_t typeOfBkgFit = 1; //0 low power - 1  Power function conv. with exponential Fit
Int_t typeOfv2BkgFit =0;
Double_t valMin=0.1396, valMax=0.158; //for DStar
//Double_t valMin=1.71, valMax = 2.06;
//TString method = "SP TPC"; //name for plots
TString method = "SP TPC "; //name for plots
TString methUsed = "SP";
TString rfp = "TPC";
TString Dmeson = "D0";
const Int_t nptBins=4;

void plotEverything(Bool_t saveHisto=kTRUE){
    
    Double_t binsPt[40] = {1.66484,1.68792,1.71099,1.73407,1.75715,1.78022,1.8033,1.82638,1.83407,1.84176,1.84946,1.85715,1.86484,1.87253,1.88022,1.88792,1.89561,1.9033,1.91099,1.92253,1.94561,1.96869,1.99176,2.01484,2.03792,2.06099+0.00384615};
    
    LoadLibraries();
    TFile *file = TFile::Open(outputFileName.Data(),"READ");
    // file->ls();
    TDirectoryFile *dir = (TDirectoryFile*)file->Get("PWGHF_D2H_HFvn_Dzero_3050_step2_QoverM_cent_SP");
    //dir->ls();
    TList *list = (TList*)dir->Get("coutputv2Dzero_3050_step2_QoverM_cent_SP");
    
    
    THnSparseD *cc2030 = (THnSparseD*)list->FindObject("hMassScalProduCQA_pt0centr300_400");
//    THnSparseD *cc30401 = (THnSparseD*)list->FindObject("hMassScalProduCQA_pt0centr300_400");
//    THnSparseD *cc2030 = (THnSparseD*)cc20301->Clone("sum");
//    cc2030->Add(cc30401);
    Double_t arrayv2[nptBins], arrayv2Err[nptBins], arraymass[nptBins], arraymassErr[nptBins], arraysigma[nptBins], arraysigmaErr[nptBins];
    Double_t arraySign[nptBins], arraySignErr[nptBins];
    Double_t arrayAvgPt[nptBins], arrayAvgPtErrL[nptBins], arrayAvgPtErrR[nptBins];
    
    for (Int_t i=0; i<4; i++){
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
    

    
  //    cc2030->GetAxis(3)->SetRangeUser(1.1,1.9);//d0-d0bar
    Double_t Etaranges[5] = {-0.8,-0.4,0,0.4,0.8};
    TH1D *mass[4];
    TH1F* hv2cpfnStandpxtot2[4];
    TH1F* hv2cpfnStandpxtot2noreb[4];
    TProfile* hv2cpfStandtot[4];
    TProfile* hv2cpfnStandtot[4];
    TH2F* hv2cStandtot[4];
    
    for(int i=0;i<4;i++){
        cc2030->GetAxis(2)->SetRangeUser(Etaranges[i],Etaranges[i+1]);//eta
        mass[i] = (TH1D*)cc2030->Projection(1);
        mass[i]->Rebin(2);
        mass[i]->SetTitle(Form("eta_range_%f_%f",Etaranges[i],Etaranges[i+1]));
        hv2cStandtot[i] = (TH2F*)cc2030->Projection(0,1);
        hv2cStandtot[i]->SetName(Form("v1%d",i));
        hv2cpfStandtot[i] = (TProfile*)hv2cStandtot[i]->ProfileX(Form("hv2cpfStandtot%d",i));
        hv2cpfnStandtot[i] = (TProfile*)hv2cpfStandtot[i]->Rebin(25, Form("hv2cpfnStandtot%d",i), binsPt);
        hv2cpfnStandpxtot2[i] = (TH1F*)hv2cpfnStandtot[i]->ProjectionX(Form("hv2cpfnStandpxtot2%d",i));
        hv2cpfnStandpxtot2[i]->SetLineColor(kRed);
        hv2cpfnStandpxtot2[i]->SetLineWidth(2);
        hv2cpfnStandpxtot2noreb[i] = (TH1F*)hv2cpfStandtot[i]->ProjectionX(Form("hv2cpfnStandpxtot2noreb%d",i));
    }
    
    
    TCanvas *ciccio[4];
    for(int i=0;i<4;i++){
        ciccio[i] = new TCanvas();
        ciccio[i]->Divide(1,2);
        ciccio[i]->cd(1);
        mass[i]->Draw();
        ciccio[i]->cd(2);
        hv2cpfnStandpxtot2noreb[i]->Draw("");
        hv2cpfnStandpxtot2[i]->Draw("same");
    }
    
    
    
    //
    //    TH1F *resoh = (TH1F*)list->FindObject("hScalProdQAQB_centr200_300");
    //    TProfile *reso = (TProfile*)resoh->Clone();
    //    Double_t meanreso = reso->GetMean();
    //    //  hv2cpfnStandpxtot2->Scale(1./TMath::Sqrt(meanreso));
    //
    //
    //
    
    //
    Int_t typeOfBkgFit = 2; //
    Int_t typeOfv2BkgFit =0;
    Double_t valMin=1.66484, valMax=2.06099+0.00384615; //for D0
    FlowFitter *f[4];
    
    for(int i=0;i<4;i++)
    {
        TH1F *maas1;
        TH1F *v11;
        
        
        //Initial parameters
        TString parName[8] = {"v_{2}^{sgn}", "M_{v2}^{bgr}", "B_{v2}^{bgr}", "A^{sgn}", "#mu", "#sigma", "M2^{bgr}", "M1^{bgr}"};
        Double_t parIni[8] =  { +0.12,  +0.2,  +0.2,  50,   1.865,  0.005, 10,   0.5}; //v2POI - M_v2 - B_v2 - IntTot - Mean - Sigma - intBkg - slope
        Double_t parMin[8] =  { -2.0,  -10.0, -10.0, 1,   1.86,  0.0001, 1,   -10.};
        Double_t parMax[8] =  { +2.0,  +10.0, +10.0, 1000,  1.869,  0.010000, 7000, 50};
        
        
        maas1 = (TH1F*)mass[i]->Clone(Form("mass%d",i));
        v11 = (TH1F*)hv2cpfnStandpxtot2[i]->Clone(Form("v1_%d",i));
        
        f[i] = new FlowFitter(maas1,v11, valMin, valMax,kFALSE,typeOfBkgFit, typeOfv2BkgFit);
        if (f[i]) cout<< "I am happy"<<endl;
        else cout<<":("<<endl;
        
        f[i]->SetFitInitialParameters(parIni);
        f[i]->SetFitMinParameters(parMin);
        f[i]->SetFitMaxParameters(parMax);
        f[i]->SetFitParNames(parName);
        f[i]->SetRangeFit(valMin,valMax);
        f[i]->SetTypeOfBkgFit(typeOfBkgFit);
        Printf("The counter i is = %d", i);
        f[i]->SetfCounter(i);
        Printf("type of v1 bkg fit is %d ====>", f[i]->GetTypeOfv2BkgFit());
        f[i]->FitMe();
        //      f->DrawFit(Form("doubleFit%d",i),i);
        f[i]->DrawFit(Form("doubleFit%d",i));
        Printf(">>>>>>fine del fit");
        if (saveHisto)  f[i]->WriteHisto("./",Form("eta%f_%f_w%d",Etaranges[i], Etaranges[i+1],1));
        
       
        arrayv2[i]=f[i]->Getv2();
        arrayv2Err[i]=f[i]->Getv2Error();
        arraymass[i]=f[i]->GetMean();
        arraymassErr[i]=f[i]->GetMeanError();
        arraysigma[i]=f[i]->GetSigma();
        arraysigmaErr[i]=f[i]->GetSigmaError();
        arraySign[i]=f[i]->GetSignificance();
        arraySignErr[i]=f[i]->GetSignificanceError();
        
        
        arrayAvgPt[i]= Etaranges[i]+0.2;
        arrayAvgPtErrR[i]= 0.2;
        arrayAvgPtErrL[i]= 0.2;
        
    }
    
    TGraphAsymmErrors *myV2 = new TGraphAsymmErrors(nptBins,arrayAvgPt,arrayv2,arrayAvgPtErrL,arrayAvgPtErrR,arrayv2Err,arrayv2Err);
    myV2->SetTitle(Form("%s v_{1} %s",Dmeson.Data(),method.Data()));
    myV2->SetName(Form("%s v_{1} %s",Dmeson.Data(),method.Data()));
    myV2->SetLineColor( kBlue );
    myV2->SetMarkerColor( kBlue );
    myV2->SetMarkerStyle( 20 );
    myV2->SetFillColor(kWhite);
    // myV2->GetYaxis()->SetRangeUser(-0.1,+0.8);
   // myV2->GetXaxis()->SetRangeUser(0.0,12.0);
    myV2->GetXaxis()->SetTitle("#eta");
    myV2->GetYaxis()->SetTitle("v_{1}");
    
    TGraphAsymmErrors *myMass = new TGraphAsymmErrors(nptBins,arrayAvgPt,arraymass,arrayAvgPtErrL,arrayAvgPtErrR,arraymassErr,arraymassErr);
    myMass->SetTitle(Form("%s mass %s",Dmeson.Data(), method.Data()));
    myMass->SetName(Form("%s mass %s",Dmeson.Data(), method.Data()));
    myMass->SetLineColor( kBlue );
    myMass->SetMarkerColor( kBlue );
    myMass->SetMarkerStyle( 20 );
    myMass->SetFillColor(kWhite);
    // myMass->GetYaxis()->SetRangeUser(-0.1,+0.8);
   // myMass->GetXaxis()->SetRangeUser(0.0,12.0);
    myMass->GetXaxis()->SetTitle("#eta");
    myMass->GetYaxis()->SetTitle("mass");
    
    TGraphAsymmErrors *mySigma = new TGraphAsymmErrors(nptBins,arrayAvgPt,arraysigma,arrayAvgPtErrL,arrayAvgPtErrR,arraysigmaErr,arraysigmaErr);
    mySigma->SetTitle(Form("%s sigma %s",Dmeson.Data(),method.Data()));
    mySigma->SetName(Form("%s sigma %s",Dmeson.Data(),method.Data()));
    mySigma->SetLineColor( kBlue );
    mySigma->SetMarkerColor( kBlue );
    mySigma->SetMarkerStyle( 20 );
    mySigma->SetFillColor(kWhite);
    // mySigma->GetYaxis()->SetRangeUser(-0.1,+0.8);
 //   mySigma->GetXaxis()->SetRangeUser(0.0,12.0);
    mySigma->GetXaxis()->SetTitle("#eta");
    mySigma->GetYaxis()->SetTitle("sigma");
    
    TGraphAsymmErrors *mySignif = new TGraphAsymmErrors(nptBins,arrayAvgPt,arraySign,arrayAvgPtErrL,arrayAvgPtErrR,arraySignErr,arraySignErr);
    mySignif->SetTitle(Form("%s Significance %s",Dmeson.Data(),method.Data()));
    mySignif->SetName(Form("%s Significance %s",Dmeson.Data(),method.Data()));
    mySignif->SetLineColor( kBlue );
    mySignif->SetMarkerColor( kBlue );
    mySignif->SetMarkerStyle( 20 );
    mySignif->SetFillColor(kWhite);
    //mySignif->GetYaxis()->SetRangeUser(-0.1,+0.8);
  //  mySignif->GetXaxis()->SetRangeUser(0.0,12.0);
    mySignif->GetXaxis()->SetTitle("#eta");
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
    
}



//
void LoadLibraries(){
    gSystem->Load("libGeom");
    gSystem->Load("libVMC");
    gSystem->Load("libXMLIO");
    gSystem->Load("libPhysics");
    
    //load needed libraries:
    gSystem->AddIncludePath("-I$ROOTSYS/include");
    //gSystem->Load("libTree");
    
    // for AliRoot
    gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libPWGflowBase");
}
//

