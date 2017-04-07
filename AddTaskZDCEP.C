AliAnalysisTask * AddTaskZDCEP(TString ZDCCalibFileName,
                               TString Label="",
                               const char* suffix="")
{
  // the manager is static, so get the existing manager via the static method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    printf("No analysis manager to connect to!\n");
    return NULL;
  }
  
  // just to see if all went well, check if the input event handler has been connected
  if (!mgr->GetInputEventHandler()) {
    printf("This task requires an input event handler!\n");
    return NULL;
  }
  
  // get the default name of the output file ("AnalysisResults.root")
  TString file = "AnalysisResults";
  
  // get the common input container from the analysis manager
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  
  // create the flow analysis tasks
  TString AnalysisTaskName = "AnalysisTaskZDCEP";
  AnalysisTaskName += Label;
  AnalysisTaskName += suffix;
  AliAnalysisTaskZDCEP *taskZDC = new AliAnalysisTaskZDCEP(AnalysisTaskName);
  taskZDC->SelectCollisionCandidates(AliVEvent::kINT7);
  
  // add list for ZDC towers gain equalization
  TString ZDCTowerEqFileName = "alien:///alice/cern.ch/user/j/jmargutt/15oHI_EZDCcalib.root";
  TFile* ZDCTowerEqFile = TFile::Open(ZDCTowerEqFileName,"READ");
  gROOT->cd();
  TList* ZDCTowerEqList = (TList*)(ZDCTowerEqFile->FindObjectAny("EZNcalib"));
  if(ZDCTowerEqList) {
    taskZDC->SetTowerEqList(ZDCTowerEqList);
    cout << "ZDCTowerEq set (from " <<  ZDCTowerEqFileName.Data() << ")" << endl;
  } else {
    cout << "ERROR: ZDCTowerEqList not found!" << endl;
    exit(1);
  }
  delete ZDCTowerEqFile;
  
  // add list for ZDC Q-vector re-centering
  TFile* ZDCCalibFile = TFile::Open(ZDCCalibFileName,"READ");
  if(!ZDCCalibFile) {
    cout << "ERROR: ZDC calibration not found!" << endl;
    exit(1);
  }
  TList* ZDCCalibList = dynamic_cast<TList*>(ZDCCalibFile->FindObjectAny("Q Vectors"));
  if(ZDCCalibList) {
    taskZDC->SetZDCCalibList(ZDCCalibList);
    cout << "ZDC calibration set (from " <<  ZDCCalibFileName.Data() << ")" << endl;
  }
  delete ZDCCalibFile;
  
  // connect the task to the analysis manager
  mgr->AddTask(taskZDC);
  
  // create a data container for the output of the flow event task
  TString taskZDCEPname = "ZDCEPExchangeContainer";
  AliAnalysisDataContainer *coutputFE = mgr->CreateContainer(taskZDCEPname,
                                                             AliFlowEventSimple::Class(),
                                                             AliAnalysisManager::kExchangeContainer);
  // connect the input data to the flow event task
  mgr->ConnectInput(taskZDC,0,cinput);
  // and connect the output to the flow event task
  mgr->ConnectOutput(taskZDC,1,coutputFE);
  
  return taskZDC;
}

