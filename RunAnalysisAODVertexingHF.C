
class AliAnalysisGrid;
class AliAnalysisAlien;

void RunAnalysisAODVertexingHF()
{
    //
    // Test macro for AliAnalysisTaskSE's for heavy-flavour candidates
    // It has the structure of a Analysis Train:
    // - in this macro, change things related to running mode
    //   and input preparation
    // - add your task using a AddTaskXXX macro
    //
    // A.Dainese, andrea.dainese@lnl.infn.it
    // "grid" mode added by R.Bala, bala@to.infn.it
    //
    
    
    gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/base -I$ALICE_PHYSICS/PWGHF/vertexingHF -I$ALICE_PHYSICS/PWG/FLOW/Base -I$ALICE_PHYSICS/PWG/FLOW/Tasks -I$ALICE_PHYSICS/PWG -I$ALICE_PHYSICS/PWGPP -g");
    
    //
    TString trainName = "D2H";
    TString analysisMode = "grid"; // "local", "grid", or "proof"
    TString inputMode    = "list"; // "list", "xml", or "dataset"
    Long64_t nentries=123567890,firstentry=0;
    Bool_t useParFiles=kFALSE;
    Bool_t useAlienPlugin=kTRUE;
    TString pluginmode="test";
    Bool_t saveProofToAlien=kFALSE;
    TString proofOutdir = "";
    TString loadMacroPath="$ALICE_PHYSICS/PWGHF/vertexingHF/macros/";
    //TString loadMacroPath="./"; // this is normally needed for CAF
    //
    
    if(analysisMode=="grid") {
        // Connect to AliEn
        TGrid::Connect("alien://");
    } else if(analysisMode=="proof") {
        // Connect to the PROOF cluster
        if(inputMode!="dataset") {printf("Input mode must be dataset, for proof analysis\n"); return;}
        gEnv->SetValue("XSec.GSI.DelegProxy","2");
        TProof::Open("alicecaf");
        //TProof::Reset("alicecaf");
        if(saveProofToAlien) {
            TGrid::Connect("alien://");
            if(gGrid) {
                TString homedir = gGrid->GetHomeDirectory();
                TString workdir = homedir + trainName;
                if(!gGrid->Cd(workdir)) {
                    gGrid->Cd(homedir);
                    if(gGrid->Mkdir(workdir)) {
                        gGrid->Cd(trainName);
                        ::Info("VertexingTrain::Connect()", "Directory %s created", gGrid->Pwd());
                    }
                }
                gGrid->Mkdir("proof_output");
                gGrid->Cd("proof_output");
                proofOutdir = Form("alien://%s", gGrid->Pwd());
            }
        }
    }
    
    
    // AliRoot libraries
    if(analysisMode=="local" || analysisMode=="grid" || analysisMode=="test") {
        TString loadLibraries="LoadLibraries.C"; loadLibraries.Prepend(loadMacroPath.Data());
        gROOT->LoadMacro(loadLibraries.Data());
        LoadLibraries(useParFiles);
    } else if (analysisMode=="proof") {
        gSystem->Load("libTree.so");
        gSystem->Load("libGeom.so");
        gSystem->Load("libPhysics.so");
        gSystem->Load("libVMC.so");
        gSystem->Load("libMinuit.so");
        // Enable the needed packages
        //gProof->ClearPackages();
        TString parDir="/afs/cern.ch/user/d/dainesea/code/";
        TString parFile;
        if(!useParFiles) {
            gProof->UploadPackage("AF-v4-17");
            gProof->EnablePackage("AF-v4-17");
            // --- Enable the PWGHFvertexingHF Package
            parFile="PWGHFvertexingHF.par"; parFile.Prepend(parDir.Data());
            gProof->UploadPackage(parFile.Data());
            gProof->EnablePackage("PWGHFvertexingHF");
        } else {
            // --- Enable the STEERBase Package
            parFile="STEERBase.par"; parFile.Prepend(parDir.Data());
            gProof->UploadPackage(parFile.Data());
            gProof->EnablePackage("STEERBase");
            // --- Enable the ESD Package
            parFile="ESD.par"; parFile.Prepend(parDir.Data());
            gProof->UploadPackage(parFile.Data());
            gProof->EnablePackage("ESD");
            // --- Enable the AOD Package
            parFile="AOD.par"; parFile.Prepend(parDir.Data());
            gProof->UploadPackage(parFile.Data());
            gProof->EnablePackage("AOD");
            // --- Enable the ANALYSIS Package
            parFile="ANALYSIS.par"; parFile.Prepend(parDir.Data());
            gProof->UploadPackage(parFile.Data());
            gProof->EnablePackage("ANALYSIS");
            // --- Enable the ANALYSISalice Package
            parFile="ANALYSISalice.par"; parFile.Prepend(parDir.Data());
            gProof->UploadPackage(parFile.Data());
            gProof->EnablePackage("ANALYSISalice");
            // --- Enable the CORRFW Package
            parFile="CORRFW.par"; parFile.Prepend(parDir.Data());
            gProof->UploadPackage(parFile.Data());
            gProof->EnablePackage("CORRFW");
            // --- Enable the PWGHFbase Package
            parFile="PWGHFbase.par"; parFile.Prepend(parDir.Data());
            gProof->UploadPackage(parFile.Data());
            gProof->EnablePackage("PWGHFbase");
            // --- Enable the PWGHFvertexingHF Package
            parFile="PWGHFvertexingHF.par"; parFile.Prepend(parDir.Data());
            gProof->UploadPackage(parFile.Data());
            gProof->EnablePackage("PWGHFvertexingHF");
        }
        gProof->ShowEnabledPackages(); // show a list of enabled packages
    }
    
    
    // Create Alien plugin, if requested
    if(useAlienPlugin) {
        if(analysisMode!="grid") {printf("Analysis mode must be grid, to use alien plugin\n"); return;}
        AliAnalysisGrid *alienHandler = CreateAlienHandler(pluginmode,useParFiles);
        if(!alienHandler) return;
    }
    
    
    //-------------------------------------------------------------------
    // Prepare input
    TChain *chainAOD = 0;
    TString dataset; // for proof
    
    if(!useAlienPlugin) {
        TString makeAODInputChain="../MakeAODInputChain.C"; makeAODInputChain.Prepend(loadMacroPath.Data());
        if(inputMode=="list") {
            // Local files
            gROOT->LoadMacro(makeAODInputChain.Data());
            chainAOD = MakeAODInputChain();// with this it reads ./AliAOD.root and ./AliAOD.VertexingHF.root
            //chainAOD = MakeAODInputChain("alien:///alice/cern.ch/user/r/rbala/newtrain/out_lhc08x/180100/",1,1);
            printf("ENTRIES %d\n",chainAOD->GetEntries());
        } else if(inputMode=="xml") {
            // xml
            gROOT->LoadMacro(makeAODInputChain.Data());
            chainAOD = MakeAODInputChain("collection_aod.xml","collection_aodHF.xml");
        } else if(inputMode=="dataset") {
            // CAF dataset
            //gProof->ShowDataSets();
            dataset="/ITS/dainesea/AODVertexingHF_LHC08x_180100";
        }
    }
    
    // Create the analysis manager
    AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
    mgr->SetDebugLevel(10);
    // Connect plug-in to the analysis manager
    if(useAlienPlugin) mgr->SetGridHandler(alienHandler);
    
    // Input
    AliAODInputHandler *inputHandler = new AliAODInputHandler("handler","handler for D2H");
    if(analysisMode=="proof" ) {
        inputHandler->AddFriend("./AliAOD.VertexingHF.root");
        //inputHandler->AddFriend("deltas/AliAOD.VertexingHF.root");
        if(saveProofToAlien) mgr->SetSpecialOutputLocation(proofOutdir);
    }
    mgr->SetInputEventHandler(inputHandler);
    //-------------------------------------------------------------------
    
    
    //-------------------------------------------------------------------
    // Analysis tasks (wagons of the train)
    //
    // First add the task for the PID response setting
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskSE *setupTask = AddTaskPIDResponse(kFALSE,kTRUE);
    
    //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
    //AliAnalysisTaskSE *PIDqaTask = AddTaskPIDqa();
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AddTaskPhysicsSelection(kFALSE,kTRUE);
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask * taskMu = AddTaskMultSelection(kFALSE);
  
//    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/AddTaskFlowQnVectorCorrectionsToLegoTrainNewDetConfig.C");
//    AddTaskFlowQnVectorCorrectionsToLegoTrainNewDetConfig("alien:///alice/cern.ch/user/p/pwg_hf/common/QnConfig/LHC15o/pass1");
  
    TString taskName;
    
    ////// ADD THE FULL D2H TRAIN
    /*taskName="../AddD2HTrain.C"; taskName.Prepend(loadMacroPath.Data());
     gROOT->LoadMacro(taskName.Data());
     Bool_t readMC=kFALSE;
     AddD2HTrain(readMC);//,1,0,0,0,0,0,0,0,0,0,0);*/
    
    ////// OR ADD INDIVIDUAL TASKS
    
    //taskName="AddTaskHFQA.C"; taskName.Prepend(loadMacroPath.Data());
    //gROOT->LoadMacro(taskName.Data());
    //AliAnalysisTaskSEHFQA *QATask = AddTaskHFQA(AliAnalysisTaskSEHFQA::kDplustoKpipi, "DplustoKpipiCuts_kINT7.root", kFALSE, kFALSE, 1, "");
    
    //D0 Spectra
//    gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/vertexingHF/charmFlow/AddTaskHFvn.C");
  gROOT->LoadMacro("AliAnalysisTaskZDCEP.cxx++");
  gROOT->LoadMacro("AddTaskZDCEP.C");
  AddTaskZDCEP("alien:///alice/cern.ch/user/j/jmargutt/15oHI_ZDCcalibVar_CenVtxCen_VtxRbR_Ecom.root");
  
  gROOT->LoadMacro("AliAnalysisTaskHFv1.cxx++");
  gROOT->LoadMacro("AddTaskHFv1.C");
  //  AliAnalysisTaskFlowD2H *DZero = AddTaskD0Mass("alien:///alice/cern.ch/user/a/adubla/D0toKpiCutsPbPb1040NoRecVtxPileupRej.root", "FirstTest_D0_PbPb", AliRDHFCuts::kD0Cuts, 2,kTRUE,kTRUE,kFALSE,kFALSE,kFALSE,1,0,0,0,0,0);
//  AliAnalysisTaskSEHFvn *DZero = AddTaskHFvn(1,"alien:///alice/cern.ch/user/j/jmargutt/HF/D0toKpi2011RefCutsPbPb1040NoRecVtxNoPileupRej_cent.root",AliAnalysisTaskSEHFvn::kD0toKpi, "D0toKpiCuts",kFALSE, "_test",AliAnalysisTaskSEHFvn::kVZERO,10.,40.,kTRUE,AliAnalysisTaskSEHFvn::kSP,"QoverM");
  AliAnalysisTaskHFv1 *DZero = AddTaskHFv1(1,kTRUE,"alien:///alice/cern.ch/user/j/jmargutt/HF/D0toKpi2011RefCutsPbPb1040NoRecVtxNoPileupRej_cent.root",AliAnalysisTaskHFv1::kD0toKpi, "D0toKpiCuts",kFALSE, "_test",AliAnalysisTaskHFv1::kZDC,10.,40.,kFALSE,AliAnalysisTaskHFv1::kSP,"QoverM");
  
    //HERE SET THE DATA MEMBRS........
    
    //taskName="AddTaskD0Mass.C"; taskName.Prepend(loadMacroPath.Data());
    //gROOT->LoadMacro(taskName.Data());
    //AliAnalysisTaskSED0Mass *D0mass = AddTaskD0Mass( 0, kFALSE, kTRUE, kFALSE, 1, 0, 0, 100, "", "D0toKpiCutsppRecVtxPileupRej_pPb.root", "D0toKpiCuts" );
    
    /*
     // attach a private task (not committed)
     // (the files MyTask.h MyTask.cxx AddMyTask.C have to be declared in plugin
     // configuration, see below)
     
     if(analysisMode.Data()=="proof") {
     gProof->LoadMacro("MyTask.cxx++g");
     } else {
     gROOT->LoadMacro("MyTask.cxx++g");
     }
     gROOT->LoadMacro("AddMyTask.C");
     MyTask *myTask = AddMyTask();
     
     
     if(analysisMode.Data()=="proof") {
     gProof->LoadMacro("AliDStarJets.cxx++g");
     } else {
     gROOT->LoadMacro("AliDStarJets.cxx++g");
     }
     gROOT->LoadMacro("AddTaskDStarJets.C");
     AliDStarJets *myTask = AddTaskDStarJets();
     */
    //-------------------------------------------------------------------
    
    //
    // Run the analysis
    //
    if(chainAOD) printf("CHAIN HAS %d ENTRIES\n",(Int_t)chainAOD->GetEntries());
    
    if(!mgr->InitAnalysis()) return;
    mgr->PrintStatus();
    if(analysisMode=="grid" && !useAlienPlugin) analysisMode="local";
    if(analysisMode!="proof") {
        mgr->StartAnalysis(analysisMode.Data(),chainAOD,nentries,firstentry);
    } else {
        // proof
        mgr->StartAnalysis(analysisMode.Data(),dataset.Data(),nentries,firstentry);
    }
    
    return;
}
//_____________________________________________________________________________
//
AliAnalysisGrid* CreateAlienHandler(TString pluginmode="test",Bool_t useParFiles=kFALSE)
{
    // Check if user has a valid token, otherwise make one. This has limitations.
    // One can always follow the standard procedure of calling alien-token-init then
    //   source /tmp/gclient_env_$UID in the current shell.
    AliAnalysisAlien *plugin = new AliAnalysisAlien();
    // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
    plugin->SetRunMode(pluginmode.Data());
//    plugin->SetUser("adubla");
    // Set versions of used packages
    plugin->SetAPIVersion("V1.1x");
   // plugin->SetROOTVersion("v5-34-30");
    plugin->SetAliPhysicsVersion("vAN-20170406-1");
    plugin->SetNtestFiles(1);
    //   gROOT->LoadMacro("$ALICE_ROOT/PWGHF/vertexingHF/AddGoodRuns.C");
    
    // Declare input data to be processed.
    //************************************************
    // Set data search pattern for DATA
    //************************************************
    //Method 1: To create automatically xml through plugin
    plugin->SetGridDataDir("/alice/data/2015/LHC15o/"); // specify LHC period
    plugin->SetDataPattern("/pass1/AOD/*AOD.root"); // specify reco pass and AOD set
    plugin->SetFriendChainName("./AliAOD.VertexingHF.root");
    // OR plugin->SetFriendChainName("deltas/AliAOD.VertexingHF.root");
    // Adds only the good runs from the Monalisa Run Condition Table
    // More than one period can be added but the period name has to be removed from GridDataDir (to be tested)
    Int_t totruns=15;
    plugin->SetRunPrefix("000");
//    plugin->AddRunNumber(195529);
//    plugin->AddRunNumber(195531);
//    plugin->AddRunNumber(195532);
//    plugin->AddRunNumber(195566);
//    plugin->AddRunNumber(195567);
//    plugin->AddRunNumber(195568);
//    plugin->AddRunNumber(195592);
//    plugin->AddRunNumber(195593);
//    plugin->AddRunNumber(195596);
//    plugin->AddRunNumber(195633);
//    plugin->AddRunNumber(195635);
//    plugin->AddRunNumber(195644);
//    plugin->AddRunNumber(195673);
//    plugin->AddRunNumber(195675);
//    plugin->AddRunNumber(195677);
  //  plugin->SetNrunsPerMaster(totruns);
    
    
    Int_t cyclenumber = 1;//5
    Int_t runcycle[] = {0,1,90};
  
    Int_t runArray[] = {246994, 246991, 246989, 246984, 246982, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246676, 246675, 246495, 246493, 246488, 246487, 246434, 246431, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246115, 246113, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245692, 245683}; // bad runs are commented out - 89 entries
    for (Int_t i =  runcycle[cyclenumber - 1]; i < runcycle[cyclenumber] ; i++)
    {
        if (i == sizeof(runArray) / sizeof(runArray[1])) break;
        plugin->AddRunNumber(runArray[i]);
    }
    plugin->SetNrunsPerMaster(totruns);

    
    //plugin->SetGridDataDir("/alice/data/2013/LHC13b"); // specify LHC period
    //   plugin->SetDataPattern("ESDs/pass3/AOD/*/AliAOD.root"); // specify reco pass and AOD set
    /*    plugin->SetFriendChainName("./AliAOD.VertexingHF.root");
     //    // OR plugin->SetFriendChainName("deltas/AliAOD.VertexingHF.root");
     //    // Adds only the good runs from the Monalisa Run Condition Table
     //    // More than one period can be added but the period name has to be removed from GridDataDir (to be tested)
     Int_t totruns=12;
     plugin->SetRunPrefix("000");
     plugin->AddRunNumber(195344);
     plugin->AddRunNumber(195346);
     plugin->AddRunNumber(195351);
     plugin->AddRunNumber(195389);
     plugin->AddRunNumber(195390);
     plugin->AddRunNumber(195391);
     plugin->AddRunNumber(195478);
     plugin->AddRunNumber(195479);
     plugin->AddRunNumber(195480);
     plugin->AddRunNumber(195481);
     plugin->AddRunNumber(195482);
     plugin->AddRunNumber(195483);
     //    //totruns += AddGoodRuns(plugin,"LHC10b"); // specify LHC period
     //    //totruns += AddGoodRuns(plugin,"LHC10c"); // specify LHC period
     //    //totruns += AddGoodRuns(plugin,"LHC10d"); // specify LHC period
     plugin->SetNrunsPerMaster(totruns);
     
     */
    // Method 2: Declare existing data files (e.g xml collections)
    
    //plugin->AddDataFile("/alice/cern.ch/user/r/rbala/000168068_000170593.xml");
    //  plugin->SetDataPattern("*AliAOD.root");
    //  plugin->SetFriendChainName("./AliAOD.VertexingHF.root");
    
    //************************************************
    // Set data search pattern for MONTECARLO
    //************************************************
    /*
     plugin->SetGridDataDir("/alice/sim/LHC10d3"); // specify MC sample
     plugin->SetDataPattern("AOD005/*AliAOD.root"); // specify AOD set
     plugin->SetFriendChainName("./AliAOD.VertexingHF.root");
     // OR plugin->SetFriendChainName("deltas/AliAOD.VertexingHF.root");
     // Adds only the good runs from the Monalisa Run Condition Table
     // More than one period can be added!
     Int_t totruns=0;
     totruns += AddGoodRuns(plugin,"LHC10b","LHC10d3"); // specify LHC period for anchor runs; and the name of the MC production
     //totruns += AddGoodRuns(plugin,"LHC10c","LHC10f7"); // specify LHC period for anchor runs;  and the name of the MC production
     //totruns += AddGoodRuns(plugin,"LHC10d","LHC10f7"); // specify LHC period for anchor runs;  and the name of the MC production
     plugin->SetNrunsPerMaster(totruns);
     */
    //
    // Define alien work directory where all files will be copied. Relative to alien $HOME.
    plugin->SetGridWorkingDir("D0v1/test0");
    // Name of executable
    plugin->SetExecutable("myHFanalysis.sh");
    // Declare alien output directory. Relative to working directory.
    plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
    // Declare the analysis source files names separated by blancs. To be compiled runtime
    // using ACLiC on the worker nodes.
    plugin->SetAnalysisSource("AliAnalysisTaskHFv1.cxx AliAnalysisTaskZDCEP.cxx");
    // Declare all libraries (other than the default ones for the framework. These will be
    // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
    plugin->SetAdditionalLibs("libPWGflowBase.so libPWGflowTasks.so libPWGHFbase.so libPWGHFvertexingHF.so libGui.so libRAWDatabase.so libCDB.so libSTEER.so libTRDbase.so libPWGTRD.so AliAnalysisTaskHFv1.cxx AliAnalysisTaskHFv1.h AliAnalysisTaskZDCEP.cxx AliAnalysisTaskZDCEP.h");
    // use par files
    if(useParFiles) {
        plugin->EnablePackage("STEERBase.par");
        plugin->EnablePackage("ESD.par");
        plugin->EnablePackage("AOD.par");
        plugin->EnablePackage("ANALYSIS.par");
        plugin->EnablePackage("OADB.par");
        plugin->EnablePackage("ANALYSISalice.par");
        plugin->EnablePackage("CORRFW.par");
        plugin->EnablePackage("PWGHFbase.par");
        plugin->EnablePackage("PWGHFvertexingHF.par");
    }
    plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/base -I$ALICE_PHYSICS/PWGHF/vertexingHF -I$ALICE_PHYSICS/PWG/FLOW/Base -I$ALICE_PHYSICS/PWG/FLOW/Tasks -I$ALICE_PHYSICS/PWG -I$ALICE_PHYSICS/PWGPP -g");

    plugin->SetSplitMaxInputFileNumber(5);
    plugin->SetNtestFiles(1);
    plugin->SetDefaultOutputs(kTRUE);
    // merging via jdl
    plugin->SetMergeViaJDL(kTRUE);
    plugin->SetOneStageMerging(kFALSE);
    plugin->SetMaxMergeStages(2);
    
    // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
    plugin->SetAnalysisMacro("AnalysisHF.C");
    // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
    // Optionally modify the name of the generated JDL (default analysis.jdl)
    plugin->SetJDLName("TaskHF.jdl");
    
    return plugin;
}
