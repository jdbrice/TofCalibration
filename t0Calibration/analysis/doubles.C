/*******************************************************************
*
*	Name:	Doubles.C
*	Description:
*		A simple script to fit trays with two or more peaks in 
*		their t0 timing resolution.
*
*
*
********************************************************************/

#include "iostream"
#include "sstream"
using namespace std;

void doubles( char* infile = "allhists.root" ) {


	gStyle->SetOptFit( 1111111);
	
	TCanvas* c1 = new TCanvas( "c1", "c1", 800, 800);
	TFile * input = new TFile( infile );

	TH2F * tHub = (TH2F*)BTofTimeRes_vs_Tray->Clone( "tHub" );

	// open the pdf file for saving the QA plots
	c1->Print( "doubleTrays.pdf[" );

	/*
	*
	* Get Tray 8
	*
	*/
	TH1D* tray8A = BTofTimeRes_vs_Tray->ProjectionY( "tray8A", 8, 8);
	TH1D* tray8B = BTofTimeRes_vs_Tray->ProjectionY( "tray8B", 8, 8);
	
	// Set the axis to a reasonable range
	tray8A->GetXaxis()->SetRangeUser( -10, 10);

	// define a fit for the desired range
	// must be set to about the correct range for the peak
	tray8A->Fit( "gaus", "", "", -10, 1);
	tray8A->Draw();

	// save the plot to the PDF
	c1->Print( "doubleTrays.pdf" );

	// set the axis to a reasonable range
	tray8B->GetXaxis()->SetRangeUser( -10, 10);

	// define a fit for the desired range
	// must be set to about the correct range for the peak
	tray8B->Fit( "gaus", "", "", 1, 10);
	tray8B->Draw();

	// print to the PDF
	c1->Print( "doubleTrays.pdf" );

	// close the pdf
	c1->Print( "doubleTrays.pdf]" );
}