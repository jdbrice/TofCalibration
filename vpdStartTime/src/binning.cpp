
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include <sys/types.h>
#include "math.h"
#include "string.h"

#ifndef __CINT__
#include "TROOT.h"
#include "Riostream.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"  
#include "TRandom.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPostScript.h"
#include "TString.h"
 #include "TLeaf.h"

using namespace std;
using std::cout;   
using std::endl;  

#endif

#include "binning.h"


Int_t binning(TChain* chain, Int_t choice)
{
  Int_t nevents = (int)chain->GetEntries();

  if(choice==0) { 
    TOFrPicoDst *t1 = new TOFrPicoDst(chain);
    cout<< " woking on binning! choice = "<<choice<<endl;
 
    for(Int_t i=0;i<nevents;i++) {
      t1->GetEntry(i);
      if ( skipBadRun( t1 ) )
        continue;

      if(i%(nevents/10)==0) cout<< i <<" events processed!"<<endl;     
      Int_t ieast = t1->numberOfVpdEast;
      Int_t iwest = t1->numberOfVpdWest;
      
      if(ieast<minpVPDHits&&iwest<minpVPDHits) continue;
     
      if(iwest>=minpVPDHits){
        for(Int_t j=0; j<nPVPDChannel/2; j++) {
		
	

	  Double_t tot = t1->vpdTotWest[j];
	  if(tot>minTOT_pVPD && tot<maxTOT_pVPD) pvpdtots[j].push_back(tot);
        }
      }
      if(ieast>=minpVPDHits){
        for(Int_t j=nPVPDChannel/2; j<nPVPDChannel; j++) {
          Double_t tot = t1->vpdTotEast[j-nPVPDChannel/2];
          if(tot>minTOT_pVPD && tot<maxTOT_pVPD) pvpdtots[j].push_back(tot);
        }
      }
    }
    
    for(Int_t i=0; i<nPVPDChannel; i++) {
      Int_t size = pvpdtots[i].size();
      cout<<"channel "<<i<<": "<< size << " hits"<<endl;
      if(size<maxpvpdNBin) {
	Double_t step = (maxTOT_pVPD-minTOT_pVPD)/maxpvpdNBin;
	for(Int_t j=0; j<=maxpvpdNBin; j++) {
          pvpdTOTBins[i][j] = step*j+minTOT_pVPD; 


	}
	cout << " pVPD Channel " << i << " seems dead ! " << "( " << size << " hits)" <<endl;
      } else {
	Int_t step = size/maxpvpdNBin; 
	std::sort(pvpdtots[i].begin(), pvpdtots[i].end());
	pvpdTOTBins[i][0] = pvpdtots[i].at(0);
	pvpdTOTBins[i][maxpvpdNBin] = pvpdtots[i].at(size-1);
	Double_t tmp;
	for(Int_t j=1; j<maxpvpdNBin; j++) {
	  pvpdTOTBins[i][j] = (pvpdtots[i].at(step*j)+pvpdtots[i].at(step*j-1))/2.0;
	  tmp = 0.;
	  for(Int_t k=step*(j-1); k<step*j; k++) {
	    tmp += pvpdtots[i].at(k);
	  }
	  pvpdTOT[i][j-1] = tmp/step;
	}
	tmp = 0;
	for(Int_t k=step*(maxpvpdNBin-1); k<size; k++) {
	  tmp += pvpdtots[i].at(k);
	}
	pvpdTOT[i][maxpvpdNBin-1] = tmp/(size-step*(maxpvpdNBin-1));
      }
    } 
    
  } else { 
    TOFrPicoDst *t2 = new TOFrPicoDst(chain);

    char rootfile[500];

    for(Int_t i=0;i<nevents;i++) {
      t2->GetEntry(i);


      Double_t t0 = t2->T0;
      Int_t hits = t2->nTofHits;
      if(hits<=0) continue;

      for(Int_t j=0; j<hits; j++) {
	Int_t tray = t2->tray[j]-1-nTray*trayStart;
	Int_t module = (t2->module[j]-1)%4;
	Int_t cell = t2->cell[j]-1;
	Int_t tcell = module*6+cell;

	if(tray>=nTray || tray<0) continue;
	if(module>=nModule || module<0) continue;
	if(cell>=nCell || cell<0) continue;	
	
	Double_t p = t2->pt[j]*TMath::CosH(t2->eta[j]);
	Double_t dedx = t2->dedx[j];
	if(p<0.3 || p>0.6 || dedx > dedxcut) continue;
        if(sqrt(pow(t2->vertexX,2)+pow(t2->vertexY,2))>rCut) continue;    

         if(TMath::Abs(t2->vertexZ-t2->vzVpd)>vzCut) continue; 


	Double_t tracklength = t2->length[j];
	Double_t betagamma = p/piMass;
	Double_t velocity = TMath::Sqrt(1.0/(1.0/betagamma/betagamma+1.0))*c_light;
	Double_t tofPi = tracklength/velocity;

	Double_t tdc = t2->leTime[j];
        tdc -= phaseDiff;
        while(tdc>51200) tdc-=51200;
	Double_t dt = tdc-t0-tofPi-tofrDelay[tray][0][tcell];

	if(TMath::Abs(dt)>toftlimit) continue;
 
	if(choice==1) { 
	  Double_t tot = t2->tot[j];

	  if(tot>minTOT_TOFr && tot<maxTOT_TOFr) tofrtots[tray][0][tcell].push_back(tot);
	  
	}
	if(choice==2) { 
	  Double_t z = t2->zLocal[j];
	  if(z>minZ_TOFr && z<maxZ_TOFr) tofrzs[tray][0][tcell].push_back(z);
	}
      }
    }
    
    for(Int_t i=0; i<nTray; i++) {
      for(Int_t j=0; j<nBoard; j++) {
	for(Int_t k=0; k<ntCell; k++) {
	  if(choice==1) {
	    Int_t size = tofrtots[i][j][k].size();
            cout << "TOFr "<<trayStart*nTray<<"th Tray "<<TDIGStart*nBoard<<"th Board "<<k<<"th Cell ( "<<size<< " hits)"<<endl;
	    if(size<maxNBin) {
	      Double_t step = (maxTOT_TOFr-minTOT_TOFr)/maxNBin;
	      for(Int_t l=0; l<=maxNBin; l++) {
		tofrTOTBins[i][j][k][l] = step*l; 
	      }
	      cout << " TOFr Channel at Tray " << trayStart*nTray << " Board " << TDIGStart*nBoard <<" Cell "<<k<< " seems dead ! (" << size << " hits)"<<endl;
	    } else {
	      Int_t step = (int)(size*1./maxNBin);
	      Int_t binl = (int)(size*0./maxNBin);
	      Int_t binh = size-binl;
	      std::sort(tofrtots[i][j][k].begin(), tofrtots[i][j][k].end());
	      tofrTOTBins[i][j][k][0] = tofrtots[i][j][k].at(binl);
	      tofrTOTBins[i][j][k][maxNBin] = tofrtots[i][j][k].at(binh-1);
	      Double_t tmp;
	      for(Int_t l=1; l<maxNBin; l++) {
		tofrTOTBins[i][j][k][l] = (tofrtots[i][j][k].at(binl+step*l)+tofrtots[i][j][k].at(binl+step*l-1))/2.0;
		tmp = 0.;
		for(Int_t m=binl+step*(l-1); m<binl+step*l; m++) {
		  tmp += tofrtots[i][j][k].at(m);
		}
		tofrTOT[i][j][k][l-1] = tmp/step;
	      }
	      tmp = 0;
	      for(Int_t m=binl+step*(maxNBin-1); m<binh; m++) {
		tmp += tofrtots[i][j][k].at(m);
	      }
	      tofrTOT[i][j][k][maxNBin-1] = tmp/(size-2*binl-step*(maxNBin-1));
	    }
	  }
	  if(choice==2) {
	    Int_t size = tofrzs[i][j][k].size();
	    if(size<maxNZBin) {
	      Double_t step = (maxZ_TOFr-minZ_TOFr)/maxNZBin;
	      for(Int_t l=0; l<=maxNZBin; l++) {
		tofrZBins[i][j][k][l] = step*l; 
	      }
	      cout << " TOFr Channel at Tray " << i << " Board " << TDIGStart*nBoard << " Cell "<<k<<" seems dead ! " << endl;
	    } else {
	      Int_t step = size/maxNZBin;
	      std::sort(tofrzs[i][j][k].begin(), tofrzs[i][j][k].end());
	      tofrZBins[i][j][k][0] = tofrzs[i][j][k].at(0);
	      tofrZBins[i][j][k][maxNZBin] = tofrzs[i][j][k].at(size-1);
	      Double_t tmp;
	      for(Int_t l=1; l<maxNZBin; l++) {
		tofrZBins[i][j][k][l] = (tofrzs[i][j][k].at(step*l)+tofrzs[i][j][k].at(step*l-1))/2.0;
		tmp = 0.;
		for(Int_t m=step*(l-1); m<step*l; m++) {
		  tmp += tofrzs[i][j][k].at(m);
		}
		tofrZ[i][j][k][l-1] = tmp/step;
	      }
	      tmp = 0;
	      for(Int_t m=step*(maxNZBin-1); m<size; m++) {
		tmp += tofrzs[i][j][k].at(m);
	      }
	      tofrZ[i][j][k][maxNZBin-1] = tmp/(size-step*(maxNZBin-1));
	    }
	  }
	}
      }
    } 
  } 

  
  if(choice==0) {
    for(Int_t i=0; i<nPVPDChannel; i++) {
      pvpdtots[i].clear();
    }
  } else {
    for(Int_t i=0; i<nTray; i++) {      
      for(Int_t j=0; j<nBoard; j++) {
	for(Int_t k=0; k<ntCell; k++) {
	  tofrtots[i][j][k].clear();
	  tofrzs[i][j][k].clear();	
	}
      }
    }  
  }
  cout<<"TOFr tot	"<<endl;
for(int tmpk=0; tmpk<maxNBin; tmpk++)
	cout<<setw(10)<<tofrTOT[0][0][0][tmpk];
  return 1;
}