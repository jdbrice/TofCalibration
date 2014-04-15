/*******************************************************************
*
*	Name:	t0Modify.C
*	Description:
*		Basic script for applying t0 offsets to a DB file
*
*
********************************************************************/


#include <iostream>
#include <fstream>
using namespace std;

static const int TRAYS = 120;
static const int MODULES = 32;
static const int BOARDS = 6;

static const int NUMBER_OF_OFFSETS = 6;

struct t0Offset {

	public:
		double t0;
		int firstTray;	// inclusive
		int lastTray;	// inclusive
		
		int firstModule = 1; // inclusive
		int lastModule = 32; // inclusive
};

t0Offset newOffsets[ NUMBER_OF_OFFSETS ];

double applyOffsets( int tray, double timing ){

	for ( int i = 0; i < NUMBER_OF_OFFSETS; i ++ ){

		if ( tray >= newOffsets[i].firstTray && tray <= newOffsets[i].lastTray ){
			return timing + newOffsets[i].t0;
		}

	}

	return timing;
}

double applyOffsets( int tray, int module, double timing ){
	for ( int i = 0; i < NUMBER_OF_OFFSETS; i ++ ){

		if ( tray >= newOffsets[i].firstTray && tray <= newOffsets[i].lastTray ){
			if ( module >= newOffsets[i].firstModule && module <= newOffsets[i].lastModule ){
				if ( tray == 8){
					cout << "Tray: " << tray << " : module: " << module << " t0 : " << newOffsets[i].t0 << endl;
				}
				return timing + newOffsets[i].t0;
			}
		}

	}

	return timing;
}


void t0Modify( char* t0In = "t0_4DB.dat", char* t0Out = "out_t0_4DB.dat") {

	for ( int i = 0; i < NUMBER_OF_OFFSETS; i++ ){
		newOffsets[ i ].firstModule = 1;
		newOffsets[ i ].lastModule = 32;
	}

	// THub 1 A
	newOffsets[ 0 ].t0 = -0.07259;
	newOffsets[ 0 ].firstTray = 1;
	newOffsets[ 0 ].lastTray = 20;
	

	// THub 1 B
	newOffsets[ 1 ].t0 = -0.07259;
	newOffsets[ 1 ].firstTray = 51;
	newOffsets[ 1 ].lastTray = 60;

	// THub 2
	newOffsets[ 2 ].t0 = -0.08924;
	newOffsets[ 2 ].firstTray = 21;
	newOffsets[ 2 ].lastTray = 50;

	// THub 3 
	newOffsets[ 3 ].t0 = -0.06149;
	newOffsets[ 3 ].firstTray = 66;
	newOffsets[ 3 ].lastTray = 95;

	// THub 3 A
	newOffsets[ 4 ].t0 = -0.03179;
	newOffsets[ 4 ].firstTray = 96;
	newOffsets[ 4 ].lastTray = 120;

	// THub 3 B
	newOffsets[ 5 ].t0 = -0.03179;
	newOffsets[ 5 ].firstTray = 61;
	newOffsets[ 5 ].lastTray = 65;

	// allocate space for the t0 offsets 
	double t0[ TRAYS ][ MODULES ][ BOARDS ];

	// open the input t0 database file and read in the values for each tray, board, module
	ifstream infile;
    infile.open( t0In );
	if ( infile.is_open() ){

		while ( !infile.eof()){
			int tray, board, module;
			infile >> tray >> module >> board;
			double timing;
			infile >> timing;
			
			if ( tray <= TRAYS && module <= MODULES && board <= BOARDS ){
				cout << "[ " << tray << " ] " << "[ " << module << " ] " << "[ " << board << " ]  = " << timing << endl;
				t0[ tray-1 ][ module-1 ][board-1] = timing;
			}

		}

	} else {
		cout << "Could not open file \n";
		cout << " New t0 values will be with respect to 0\n"; 
	}


	// now apply the offsets and write out a new file:
	ofstream out( t0Out );
	for ( int it = 0; it < TRAYS; it ++ ){
		for ( int im = 0; im < MODULES; im ++ ){
			for ( int ib = 0; ib < BOARDS; ib ++ ){
				
				double ct = t0[it][im][ib];

				// array is zero indexed but tray, module, board are indexed at 1
				// apply the offsets
				ct = applyOffsets( it+1, im+1, ct );

				// output the tray, module, board
				out << (it+1) << " " << (im+1) << " " << (ib+1) << endl;

				// output the timing offset for this tray, module, board
				out << ct << endl;

			}		
		}
	}

	out.close();

	return;
}

