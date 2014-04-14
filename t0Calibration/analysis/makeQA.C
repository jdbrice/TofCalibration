



void makeQA(){

	gStyle->SetOptStat( 1111 );
	gStyle->SetStatX( 0.35 );
	TCanvas* can = new TCanvas( "QA Plots", "", 800, 800 );

	TFile * input = new TFile( "allhists.root" );

	can->Print( "t0QA.pdf[" );


	TF1* gaus = new TF1( "g", "gaus" );
	gaus->SetRange( -.5, .5 );

	BTofTimeRes_vs_Tray->FitSlicesY( gaus );

	gPad->SetLeftMargin( .15 );
	gPad->SetRightMargin( .15 );
	BTofTimeRes_vs_Tray->SetTitle( "BTof Time Resolution vs Tray" );
	BTofTimeRes_vs_Tray->GetXaxis()->SetTitle( "Tray" );
	BTofTimeRes_vs_Tray->GetYaxis()->SetTitleOffset( 2 );
	BTofTimeRes_vs_Tray->GetYaxis()->SetTitle( "Time (ns) " );
	BTofTimeRes_vs_Tray->GetYaxis()->SetRangeUser( -.5, .5);
	BTofTimeRes_vs_Tray->Draw("colz");
	can->Print( "t0QA.pdf" );

	gStyle->SetOptStat(0);

	gPad->SetRightMargin( .1 );
	BTofTimeRes_vs_Tray_1->SetTitle( "Mean For Each Tray" );
	BTofTimeRes_vs_Tray_1->GetXaxis()->SetTitle( "Tray" );
	BTofTimeRes_vs_Tray_1->GetYaxis()->SetTitleOffset( 2 );
	BTofTimeRes_vs_Tray_1->GetYaxis()->SetTitle( "Time (ns) " );
	BTofTimeRes_vs_Tray_1->GetYaxis()->SetRangeUser( -.004, .004 );
	BTofTimeRes_vs_Tray_1->Draw("ep");
	can->Print( "t0QA.pdf" );


	BTofTimeRes_vs_Tray_2->SetTitle( "Sigma For Each Tray" );
	BTofTimeRes_vs_Tray_2->GetXaxis()->SetTitle( "Tray" );
	BTofTimeRes_vs_Tray_2->GetYaxis()->SetTitleOffset( 2 );
	BTofTimeRes_vs_Tray_2->GetYaxis()->SetTitle( "Time (ns) " );
	BTofTimeRes_vs_Tray_2->GetYaxis()->SetRangeUser( .12, .16 );
	BTofTimeRes_vs_Tray_2->Draw("ep");
	can->Print( "t0QA.pdf" );



	can->Print( "t0QA.pdf]" );
}