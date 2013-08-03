#ifndef MYHISTO__HH
#define MYHISTO__HH 1

#include <string>
#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector>
#include <cstdio>
#include "Pythia.h"


class Histogram{

 public:

  //open canvas
  TCanvas *mycanvas = new TCanvas("mycanvas","pseudotop",200,10,700,500);
  mycanvas->SetFillColor(4000);
  mycanvas->SetGrid();

  //divides canvas to a 2X3 set of pads
  mycanvas->Divide(2,3);

  //Book histograms 
  TH1F* top_pmass = new TH1F("m_ptop"," top mass from pythia ;mass [GeV] ",100,150.,190.);

  //tops

  TH1F* top_y   = new TH1F("y_top",";#it{y}^{top}", 100, -5., 5.);

  //pseudotop
  TH1F* top_reconmass = new TH1F("reconmass_top","reconstructed mass of pseudotop", 100,130,190.);
  TH1F* trecon_pt = new TH1F("pT_tr","pseudotop  pt", 100, 0, 300.0);  
  TH1F* trecon_eta = new TH1F("eta_tr","pseudotop eta", 100, -5., 5.0);
  TH1F* trecon_y = new TH1F("y_tr","pesudotop y;it{p}_{T} [GeV]",100,-5.,5.);
 
  //W 
  TH1F* whad = new TH1F("whadronic"," hadronic pseudo W; mass [GeV]", 80,40.,116.);
  TH1F* wlep = new TH1F("wleptonic"," leptonic w mass from decay products; mass [GeV]", 100,40.,160.);
  TH1F* wlje = new TH1F("wlje","leptonic PseudoW; mass [GeV]", 100,40.,190.);
  


  //make tpads for drawing and displaying on canvas
  Tpad * t1 = new Tpad("pseudo_pt", " pt_of_psedotop", 0,100, 0,100 );
  Tpad * t1 = new Tpad("pseudo_pt", " pt_of_psedotop", 0,100, 0,100 );
  Tpad * t1 = new Tpad("pseudo_pt", " pt_of_psedotop", 0,100, 0,100 );
  Tpad * t1 = new Tpad("pseudo_pt", " pt_of_psedotop", 0,100, 0,100 );
  Tpad * t1 = new Tpad("pseudo_pt", " pt_of_psedotop", 0,100, 0,100 );
  Tpad * t1 = new Tpad("pseudo_pt", " pt_of_psedotop", 0,100, 0,100 );

  TH1F* top_mass = new TH1F("m_top"," top mass from decay partons ;mass [GeV] ",100,150.,190.);
  top_mass->SetMarkerColor(kBlue); 
  top_mass->SetFillColor(0);
  top_mass->SetLineColor(4);
  top_mass->Draw("AC*");
 
  top_mass->cd(1);
  
  TH1F* top_pt = new TH1F("pT_top","top quark pT;#it{p}_{T} [GeV];Frequency", 100, 0.  , 300.);



  TH1F* top_eta = new TH1F("eta_top","pseudrapitity of t;#eta^{top}", 100, -5., 5.);
  TGraph* top_y = new TGraph("topy", " top rapidity vs matching ration", 100, 
   


 
 

 void fill_it( ) {


}

 


 private:
 
 

 
 
   
};
#endif





