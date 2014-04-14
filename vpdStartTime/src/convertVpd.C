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
  //const Float_t vzOffset = 2454.21696153 + 0.114676;
  
  //const Float_t vzOffset = 2454.21696153-6.7258;
  //const Float_t vzOffset = 2454.21696153-6.77748;
  //const Float_t vzOffset = 2454.21696153-6.9333-0.0732;
  //const Float_t vzOffset = 2454.21696153-7.015-0.064; run 11051
  //const Float_t vzOffset = 2454.21696153-7.18157; //run 11065
//  const Float_t vzOffset = 2454.21696153-9.68563; //run 11087
//  const Float_t vzOffset = 2454.21696153-9.45167; //39earlyproduction(testsample)
//  const Float_t vzOffset = 2454.21696153-9.36829; //39earlyproduction(kefengs)
//  const Float_t vzOffset = 2454.21696153-9.43839;//39earlyhalflist

//  const Float_t vzOffset = 4542.37-2.02792-0.0030382;//run11 hlt
//  const Float_t vzOffset = 4540.34208-0.019722;//run11 hlt
//-6.65938e-03
//const Float_t vzOffset = 4542.37-2.02792-4.59159-0.0259989+9.59631e-02+6.65938e-03;//run11 auau19 calibration
//  const Float_t vzOffset = 4542.37-2.02792-4.59159+9.59631e-02+3.88363e+00;//auau200rff
//  const Float_t vzOffset = 4542.37-2.02792-4.59159+9.59631e-02+9.91213;//run 12 pp with run11 information. 
    const Float_t vzOffset = 4542.37-2.02792-4.59159+9.59631e-02 - 4487 - 0.5076;//run 13 - daniel. 
  const Float_t c_light = 29.9792458;

  Int_t nChannel;
  Float_t tot[maxN], corr[maxN];  
  Float_t totbound[maxN+1], corrbound[maxN+1];
  ifstream indata;
  indata.open("pvpdCali.dat");

  ofstream outdata;
  outdata.open("pvpdCali_4DB.dat");
  for(int i=0;i<nVPDChannel;i++) {
    indata>>nChannel;
	//cout<<nChannel<<endl;
    for(int j=0;j<=maxN;j++){
      indata>>totbound[j];
	//  cout<<totbound[j]<<"	";
    }
	//cout<<endl;
    for(int j=0;j<=maxN;j++){
      indata>>corrbound[j];  
	//  cout<<corrbound[j]<<"	";	  
      //vzOffset
      if(i>=nVPDChannel/2) {
        corrbound[j] -= 2.*vzOffset/c_light;
      }
    }
	//cout<<endl;
    for(int j=0;j<maxN;j++){
      indata>>tot[j];
	//	  cout<<tot[j]<<"	";	  
    }
	//	cout<<endl;
    for(int j=0;j<maxN;j++){
      indata>>corr[j];  
	//  cout<<corr[j]<<"	";	 	  
      //vzOffset
      if(i>=nVPDChannel/2) {
        corr[j] -= 2.*vzOffset/c_light;
      }
    } 
	//cout<<endl;
    //output
    outdata<<setw(4) <<i+1<<endl;// change to tubeID instead of index
    outdata<<setw(4) << maxN<<endl;
    for(int j=0;j<=maxN;j++){
      outdata<<setw(16)<<totbound[j]<<flush;
    }
    outdata<<endl;
    for(int j=0;j<=maxN;j++){
      outdata<<setw(16)<<corrbound[j]<<flush;
    }
    outdata<<endl;
  }
  outdata.close();

  cout<<"done!"<<endl;
  return;
}
  

