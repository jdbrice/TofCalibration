/*************************************************************************
*
*	Name:	fitTHubs.C
*	Description:
*		Fits each of the thub ranges to a 0 order polynomial
*		and reports the value 
*
*
*************************************************************************/

#include "iostream"
#include "sstream"
using namespace std;

void thubs( char* infile = "allhists.root" ) {


	gStyle->SetOptFit( 1111111);
	
	TCanvas* c1 = new TCanvas( "c1", "c1", 800, 800);

	// open the PDF for plot output
	c1->Print( "thub.pdf[" );

	TFile * input = new TFile( infile );

	TH2F * tHub = (TH2F*)BTofTimeRes_vs_Tray->Clone( "tHub" );	
	tHub->FitSlicesY();


	/*
	*	THUB 1
	*/
	TH1F * th1 = new TH1F( "th1", "", 120, 0, 120);

	for ( int bi = 0; bi < 120; bi ++ ){

		if ( ((bi >= 0 && bi <= 20 ) || (bi >= 51 && bi <= 60 )) && 
			!( bi == 8 || bi == 23 || bi == 93 || bi == 108 || bi == 193 ))
			th1->SetBinContent( bi, tHub_1->GetBinContent( bi ) );

	}

	cout << "1-20 & 51-60" << endl << endl;
	TFitResultPtr th1Fit = th1->Fit( "pol0", "W");

	th1->SetTitle( "1-20 & 51-60");
	th1->Draw();
	c1->Print( "thub.pdf" );

	/*
	*	THUB 2
	*/
	TH1F * th2 = new TH1F( "th2", "", 120, 0, 120);

	for ( int bi = 0; bi < 120; bi ++ ){

		if ( ((bi >= 21 && bi <= 50 ) ) && 
			!( bi == 8 || bi == 23 || bi == 93 || bi == 108 || bi == 193 ))
			th2->SetBinContent( bi, tHub_1->GetBinContent( bi ) );

	}

	cout << "21-50" << endl << endl;
	TFitResultPtr th2Fit = th2->Fit( "pol0", "W");

	th2->SetTitle( "21-50");
	th2->Draw();
	c1->Print( "thub.pdf" );

	/*
	*	THUB 3
	*/
	TH1F * th3 = new TH1F( "th3", "", 120, 0, 120);

	for ( int bi = 0; bi < 120; bi ++ ){

		if ( ((bi >= 66 && bi <= 95 ) ) && 
			!( bi == 8 || bi == 23 || bi == 93 || bi == 108 || bi == 193 ))
			th3->SetBinContent( bi, tHub_1->GetBinContent( bi ) );

	}

	cout << "66-95" << endl << endl;
	TFitResultPtr th3Fit = th3->Fit( "pol0", "W");

	th3->SetTitle( "66-95");
	th3->Draw();
	c1->Print( "thub.pdf" );

	/*
	*	THUB 4
	*/
	TH1F * th4 = new TH1F( "th4", "", 120, 0, 120);

	for ( int bi = 0; bi < 120; bi ++ ){

		if ( ((bi >= 96 && bi <= 120 ) || (bi >= 61 && bi <= 65 )) &&
			!( bi == 8 || bi == 23 || bi == 93 || bi == 108 || bi == 119 ))
			th4->SetBinContent( bi, tHub_1->GetBinContent( bi ) );

	}

	cout << "96-120 & 61-65" << endl << endl;
	TFitResultPtr th4Fit = th4->Fit( "pol0", "W");

	th4->SetTitle( "96-120 & 61-65");
	th4->Draw();
	c1->Print( "thub.pdf" );

	//close the pdf 
	c1->Print( "thub.pdf]" );

}