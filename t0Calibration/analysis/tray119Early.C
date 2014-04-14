
#include "iostream"
#include "sstream"
using namespace std;


static const int TRAYS = 120;
static const int MODULES = 32;
static const int BOARDS = 6;

static const int NUMBER_OF_OFFSETS = 1;

struct t0Offset {

	public:
		double t0;
		int firstTray;	// inclusive
		int lastTray;	// inclusive
		
		int firstModule = 1; // inclusive
		int lastModule = 32; // inclusive

		std::vector<int> excludeTrays;
};

t0Offset newOffsets[ NUMBER_OF_OFFSETS ];


double applyOffsets( int tray, int module, double timing ){

	for ( int i = 0; i < NUMBER_OF_OFFSETS; i ++ ){

		if ( tray >= newOffsets[i].firstTray && tray <= newOffsets[i].lastTray ){
			if ( module >= newOffsets[i].firstModule && module <= newOffsets[i].lastModule ){
				return timing + newOffsets[i].t0;
			}
		}

	}

	return timing;
}

void tray119Early() {


	gStyle->SetOptFit( 1111111);
	
	TCanvas* c1 = new TCanvas( "c1", "c1", 800, 800);

	TFile * input = new TFile( "allhists.root" );

	TH2F * tHub = (TH2F*)BTofTimeRes_vs_Tray->Clone( "tHub" );

	c1->Print( "doubleTrays.pdf[" );

	/*
	*
	* Tray 119
	*
	*/
	TH1D* tray119A = BTofTimeRes_vs_Tray->ProjectionY( "tray119A", 119, 119);
	TH1D* tray119B = BTofTimeRes_vs_Tray->ProjectionY( "tray119B", 119, 119);
	
	tray119A->GetXaxis()->SetRangeUser( -10, 10);
	tray119A->Fit( "gaus", "", "", -1.30, -.60);
	tray119A->Draw();
	c1->Print( "doubleTrays.pdf" );

	tray119B->GetXaxis()->SetRangeUser( -10, 10);
	tray119B->Fit( "gaus", "", "", -.5, 2);
	tray119B->Draw();
	c1->Print( "doubleTrays.pdf" );



	/*
	*
	*	Extract the mean from each fit and print as a double check
	*/


	TF1 * fitTray119A = tray119A->GetFunction( "gaus" );
	cout << "Tray119 Modules 1-16: " << fitTray119A->GetParameter(1) << endl;
	TF1 * fitTray119B = tray119B->GetFunction( "gaus" );
	cout << "Tray119B Modules 17-32: " << fitTray119B->GetParameter(1) << endl;


	/*
	*
	*	Create the offset for each of the trays' module subreagion
	*
	*/

	/*
	*	Tray 119 FIT TO SECOND PEAK s
	*/ 
	newOffsets[ 0 ].t0 = fitTray119A->GetParameter(1);
	newOffsets[ 0 ].firstTray = 119;
	newOffsets[ 0 ].lastTray = 119;
	newOffsets[ 0 ].firstModule = 1;
	newOffsets[ 0 ].lastModule = 32;

	// close the pdf file to ensure it gets flushed out.
	c1->Print( "doubleTrays.pdf]" );

	cout << "Applying offsets " << endl;
	applyOffsets();


	
}



void applyOffsets(){

	double t0[ TRAYS ][ MODULES ][ BOARDS ];

	/*
	*	Should have already applied the tray by tray offsets
	*/	
	ifstream infile;
    infile.open("t0_4DB.dat");
	if ( infile.is_open() ){

		while ( !infile.eof()){
			int tray, board, module;
			infile >> tray >> module >> board;
			double timing;
			infile >> timing;
			
			if ( tray <= TRAYS && module <= MODULES && board <= BOARDS ){
				//cout << "[ " << tray << " ] " << "[ " << module << " ] " << "[ " << board << " ]  = " << timing << endl;
				t0[ tray-1 ][ module-1 ][board-1] = timing;
			}

		}

	} else {
		cout << "Could not open file \n";
	}


	// now apply the offsets and write out a new file:
	ofstream out( "tray119Early_t0_4DB.dat");
	for ( int it = 0; it < TRAYS; it ++ ){
		for ( int im = 0; im < MODULES; im ++ ){
			for ( int ib = 0; ib < BOARDS; ib ++ ){
				
				double ct = t0[it][im][ib];

				ct = applyOffsets( it+1, im+1, ct );

				out << (it+1) << " " << (im+1) << " " << (ib+1) << endl;
				out << ct << endl;

			}		
		}
	}

	out.close();




}