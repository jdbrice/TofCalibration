#include "iostream.h"

class     StChain;
StChain  *chain=0;
class     St_db_Maker;
St_db_Maker *dbMk =0;
//MAKE SURE THIS IS INDEED STARTLESS
Int_t iEvt=0,istat=0,nEvents=0;
void t0Calibration(const Char_t *fileList = "test.lis",
				  const Char_t *histname = "test_hist_calibFromTable.root")
//                   const Char_t *ntuplename = "test.ntuple.root")
{
  Int_t nEvents = 10000000;
  Int_t nfiles = 1000;
  //
  // First load some shared libraries we need
  //
	cout<<"Loading Filelist: "<< fileList<<endl;
  if (gClassTable->GetID("TTable") < 0) {
    gSystem->Load("libStar");
    gSystem->Load("libPhysics");
  }  
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("StarMagField");
  gSystem->Load("StMagF");
  gSystem->Load("StDetectorDbMaker");
  gSystem->Load("StTpcDb");
  //  gSystem->Load("StDbUtilities");
  gSystem->Load("StDaqLib");
  gSystem->Load("StDbBroker");
//  gSystem->Load("StDetectorDbMaker");
  gSystem->Load("StDbUtilities");
  gSystem->Load("St_db_Maker");

  // rafael - begin
  gSystem->Load("StEvent");
  gSystem->Load("StEventMaker");
  gSystem->Load("StarMagField");
 
  gSystem->Load("libtpc_Tables");
  gSystem->Load("libGeom");
  gSystem->Load("St_g2t");
  gSystem->Load("xgeometry");
  gSystem->Load("St_geant_Maker"); 
  // rafael - end

  gSystem->Load("StBTofUtil");
  //  gSystem->Load("StBTofMaker");
  gSystem->Load("StBTofMatchMaker");
  gSystem->Load("StVpdCalibMaker");
  gSystem->Load("StBTofCalibMaker");
  gSystem->Load("StBTofQAMaker");


  // Handling depends on whether file is a ROOT file or XDF file
  //
  chain  = new StChain("StChain");
  
  StMuDstMaker *muDstMaker = new StMuDstMaker(0,0,"",fileList,"MuDst.root",nfiles);
  // rafael - begin
  muDstMaker->SetStatus("*",0);
  muDstMaker->SetStatus("MuEvent*",1);
  muDstMaker->SetStatus("PrimaryVertices*",1);
  muDstMaker->SetStatus("PrimaryTrack*",1);
  muDstMaker->SetStatus("GlobalTrack*",1);
  muDstMaker->SetStatus("BTof*",1);
  // rafael - end  
  
  cout<<endl<<"============  Data Base ========="<<endl;
  dbMk = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","StarDb");
  
  // rafael - begin
  St_geant_Maker *geantMk = new St_geant_Maker();
  geantMk->LoadGeometry("detp geometry y2014");
  geantMk->SetActive(kFALSE);
  
  StBTofMatchMaker *matchMaker = new StBTofMatchMaker("btofMatch");

   matchMaker->SetDebug(0);
   matchMaker->setMuDstIn(kTRUE);
   matchMaker->setAlignFileName("./geometry.txt");
   matchMaker->setIdealGeometry(kFALSE);
   matchMaker->setCalculateAlign(kTRUE);


  // rafael - end
  
  StVpdCalibMaker *vpdCalib = new StVpdCalibMaker("vpdCalib");
  vpdCalib->SetDebug(1);
  vpdCalib->setMuDstIn();
  vpdCalib->setUseVpdStart(kFALSE);//should be false for startless
  vpdCalib->setInitFromFile(kTRUE);
  vpdCalib->setCalibFilePvpd("./db/pvpdCali_4DB.dat");
  
  StBTofCalibMaker *btofCalib = new StBTofCalibMaker("btofCalib");
  btofCalib->setMuDstIn(kTRUE);
  btofCalib->setInitFromFile(kTRUE);
  btofCalib->setCalibFileTot( "./db/totCali_4DB.dat");
  btofCalib->setCalibFileZhit("./db/zCali_4DB.dat");
  btofCalib->setCalibFileT0(  "./db/t0_4DB.dat");
  

  StBTofQAMaker *qaMaker = new StBTofQAMaker(muDstMaker,histname);
  qaMaker->setHistogramOn(kTRUE);
  
  // Initialize chain
  //
  Int_t iInit = chain->Init();
  if (iInit) chain->Fatal(iInit,"on init");
  chain->PrintInfo();
  //
  // Event loop
  //
  int istat = 0, i = 1;
 EventLoop: if (i <= nEvents && istat != 2) {
    
    cout << endl << "============================ Event " << i
	 << " start ============================" << endl;
    
    chain->Clear();
    istat = chain->Make(i);
    if (istat == 2) 
      {cout << "Last  event processed. Status = " << istat << endl; break;}
    if (istat == 3) 
      {cout << "Error event processed. Status = " << istat << endl; break;}
    
    //   gObjectTable->Print();
    i++;
    goto EventLoop;
  }
  
  i--;
  cout<<endl<<"============================ Event "<<i<<" finish ============================"<<endl;
  
  //
  // Chain Finish
  if (nEvents > 1) {
    chain->Finish();
  }
  
  delete chain;
}



