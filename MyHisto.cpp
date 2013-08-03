#include "MyHisto.h"
#include "pythia.h"
#include "fjClustering.h"
#include "root.h"

//open canvas
TCanvas *mycanvas = new TCanvas("mycanvas","pseudotop",200,10,700,500);
mycanvas->SetFillColor(4000);
mycanvas->SetGrid();

mycanvas->Divide(2,3)


//make tpads for drawing and displaying on canvas
Tpad * t1 = new Tpad("pseudo_pt", " pt_of_psedotop", 0,100, 0,100 );
Tpad * t2 = new Tpad("pseudo_y", " y_of_psedotop", 0,100, 0,100 );
Tpad * t3 = new Tpad("pseudo_m", " m_of_psedotop", 0,100, 0,100 );
Tpad * t4 = new Tpad("pseudo_eta", " eta of pseudotop", 0,100, 0,100 );
Tpad * t5 = new Tpad("top_pt", "top pt ", 0,100, 0,100 );
Tpad * t6 = new Tpad("top_y", " top y", 0,100, 0,100 );


//get histograms from another .root files
 
TFile *flMC_file = new TFile("ttbar_partonMode.root");
TFile* file = TFile::Open("ttbar_partonMode.root", "RECREATE"); 
TH1F* h1dMC_NN = (TH1F*)flMC_file->Get("h1dNN");
TH1F* h1dMC_NN = (TH1F*)flMC_file->Get("h1dNN");
TH1F* h1dMC_NN = (TH1F*)flMC_file->Get("h1dNN");
TH1F* h1dMC_NN = (TH1F*)flMC_file->Get("h1dNN");
TH1F* h1dMC_NN = (TH1F*)flMC_file->Get("h1dNN");
TH1F* h1dMC_NN = (TH1F*)flMC_file->Get("h1dNN");
TH1F* h1dMC_NN = (TH1F*)flMC_file->Get("h1dNN");




flMC_file->GetList()Write(); //write in memory objects of 1st file to current file





//write to and close root file
outFile->Write();
outFile->Close();


