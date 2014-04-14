#include <stdio>
#include <stdlib>
#include <fstream>
#include <iostream>
#include "iomanip.h"

void convertVpd()
{
  gROOT->Reset();
  gStyle->Reset("plain");
  
  const Int_t maxN = 50;
  const Int_t nVPDChannel = 38;
  const Float_t c_light = 29.9792458;

  const Float_t vzOffset = 9.67581530999990633e+01; // run14 auau15GeV

  Int_t nChannel;
  Float_t tot[maxN], corr[maxN];  
  Float_t totbound[maxN+1], corrbound[maxN+1];
  ifstream indata;
  indata.open("pvpdCali.dat");

  ofstream outdata;
  outdata.open("pvpdCali_4DB.dat");
  for(int i=0;i<nVPDChannel;i++) {
    indata>>nChannel;
	
    for(int j=0;j<=maxN;j++){
      indata>>totbound[j];
	
    }
	
    for(int j=0;j<=maxN;j++){
      indata>>corrbound[j];  
	
      //vzOffset
      if(i>=nVPDChannel/2) {
        corrbound[j] -= 2.*vzOffset/c_light;
      }
    }

    for(int j=0;j<maxN;j++){
      indata>>tot[j];

    }

    for(int j=0;j<maxN;j++){
      indata>>corr[j];  

      //vzOffset
      if(i>=nVPDChannel/2) {
        corr[j] -= 2.*vzOffset/c_light;
      }
    } 
	
    //output
    outdata<<setw(4) <<i+1<<endl;// change to tubeID instead of index
    outdata<<setw(4) << maxN<<endl;
    for(int j=0;j<=maxN;j++){
      outdata << setw(16) << totbound[j] << flush;
    }
    outdata<<endl;
    for(int j=0;j<=maxN;j++){
      outdata<<setw(16) << corrbound[j] << flush;
    }
    outdata << endl;
  }
  outdata.close();

  cout << "done!" << endl;
  
  return;
}
  

