// Stdlib header file for input and output.
#include <iostream>

// Header file to access Pythia 8 program elements.
#include "Pythia.h"

//root 
#include "root.h"

//fastjet 
#include "fjClustering.h"

//Jets
#include "JetSeparation.h"
#include "MyTopEvent.h"

using namespace Pythia8; 


//functions 
bool isBhadron(int);
bool isElectron(const Particle &ptcl) { return ptcl.idAbs()==11;}
bool isTau(const Particle &ptcl) { return ptcl.idAbs() ==15;}
bool isMu(const Particle &ptcl) { return ptcl.idAbs()==13;}
bool isNu(const Particle &ptcl) { return ptcl.idAbs()==12 || ptcl.idAbs()==14 || ptcl.idAbs()==16; }
void fatal(TString msg) { printf("ERROR:\n\n  %s\n\n",msg.Data()); abort(); }
void PrintPtcl(const Particle &ptcl, TString comment="") {
  printf("  (pT,eta,phi,m) = (%6.1f GeV,%5.2f,%5.2f,%5.1f GeV) pdgID %4d : %s\n",
    ptcl.pT(),ptcl.eta(),ptcl.phi(),ptcl.m(),ptcl.id(),comment.Data());
} 
void PrintPseudojets( const vector<fastjet::PseudoJet> &);

//main function
int  main(int argc, char* argv[]) {


  bool partonMode=true;
  bool debug= false;

  // Create the ROOT application environment. 
  TApplication theApp("hist", &argc, argv);

  //generator without banner
  Pythia pythia("",false);
  
  //read settings from command file
  // if(partonMode){pythia.readFile("partonMode.cmnd");} else { pythia.readFile("particleMode.cmnd");}
  
  pythia.initTuneEE(2);
  pythia.readString("SoftQCD:minBias = off");
  pythia.readString("SoftQCD:singleDiffractive = off");
  pythia.readString("SoftQCD:doubleDiffractive = off");
  pythia.readString("HardQCD:all = off");
  pythia.readString("PartonLevel:FSR = off");
  pythia.readString("PartonLevel:ISR = on ");
  pythia.readString("HadronLevel:all = on");
  pythia.readString("Top:qqbar2ttbar = on");
  pythia.readString("Top:gg2ttbar = on");
  if (!partonMode) {
    pythia.readString("SoftQCD:minBias = on ");
    pythia.readString("HardQCD:all = on");
    pythia.readString("PartonLevel:FSR = on");
    pythia.readString("PartonLevel:ISR = on ");
    pythia.readString("HadronLevel:all = on ");
  } 


  //beam initialization (proton,proton, 7000 GeV)
  pythia.init(2212,2212,7000);
  
  
  //Create file on which histograms will be saved.
  TString of_name = partonMode ? "ttbar_partonLevel_histos_"+string(argv[1])+".root" : "ttbar_hadronLevel_histograms.root";
  TFile* outFile = new TFile("rootplots/"+of_name, "RECREATE");

  //_____________________________________________________________________________________

  
  //truth 
  TH1F* top_mass  = new TH1F("m_top"," truth mass ; mass [GeV] ",20,130.,190.);
  TH1F* top_pt    = new TH1F("pt_top","ttruth pT ;#it{p}_{T} [GeV];Frequency", 50, 0.  , 300.);
  TH1F* top_eta   = new TH1F("eta_top","truth #eta ;#eta^{top}", 20, -5., 5.);
  TH1F* top_y     = new TH1F("y_top","truth rapidity ;#it{y}^{top}", 20, -5., 5.);

  //not matched pseudos
  TH1F* pseudo_nomatch_mass  = new TH1F("pseudo_mass_unmatched"," unmatched leptonic pseudotop mass; mass [GeV]", 50 ,100,250.);
  TH1F* pseudo_nomatch_eta   = new TH1F("pseudo_eta_unmatched"," unmatched leptonic pseudotop eta; #eta", 50 ,-6,6.);
  TH1F* pseudo_nomatch_pt    = new TH1F("pseudo_pt_unmatched"," unmatched leptonic pseudotop pt;#it{p}_{T} [GeV]", 50 ,0,330.);

  //matched pseudos
  TH1F* pseudo_match_mass  = new TH1F("pseudo_mass_matched"," matched leptonic pseudotop; mass [GeV]", 50 ,100,250.);
  TH1F* pseudo_match_pt    = new TH1F("pseudo_pt_matched",  " matched leptonic pseudotop; #it{p}_{T} [GeV]", 50 ,0,330.);
  TH1F* pseudo_match_eta   = new TH1F("pseudo_eta_matched"," matched leptonic pseudotop ; #eta ", 50 ,-6,6.);

  //different Bjet-W pairs
  TH1F* bestbtop      = new TH1F("pseudo_btop_mass"," best b with lep W pseudotop mass; mass [GeV]", 50 ,100,230.);
  TH1F* withhad       = new TH1F("pseudo_HW_top"," pseudotop mass using had w and best b jet; mass [GeV]", 50 ,100,230.);
  
  //hadronic and leptonic Ws
  TH1F* whad = new TH1F("whadronic"," hadronic pseudo W; mass [GeV]", 20,40.,120.);
  TH1F* wlep = new TH1F("wleptonic"," leptonic pseudo w ; mass [GeV]", 20,40.,120.);
  
  //matching efficiency plots  
  TProfile *matchingEff_vs_topPt   =  new TProfile("eff_vs_topPt"," Matching efficiency vs pT #it{p}_{T} ; p_{T} [Gev]; efficiency",50,0,300);
  TProfile *nmatchingEff_vs_topPt  =  new TProfile("neff_vs_topPt","Matching inefficiency vs pT #it{p}_{T} ; p_{T} [Gev] ;  miss-match ",50,0,300);
  TProfile *matchingEff_vs_topmass =  new TProfile("eff_vs_topMass","Matching Efficiency vs mass_{top} ; mass [GeV]; efficieny ",20,130,190);
  TProfile *matchingEff_vs_topeta  =  new TProfile("eff_vs_topeta","Matching efficiency vs #eta ; #eta ; efficiency",20,-5,5);
  TProfile *matchingEff_vs_bmass   =  new TProfile("eff_vs_bmass", "Matching efficiency vs mass_{b}; mass [GeV]; efficiency", 5, 0, 20);

  TProfile *matchingEff_vs_tbmass  =  new TProfile("match_ratio_vs_mass", "Matching ratio top/tbar; mass [GeV]; efficiency", 5, 0, 180);

  TProfile *tbar_matching          =  new TProfile("eff_vs_ttbarmass", "Matching efficiency vs mass_{tbar}; mass [GeV]; efficiency", 5, 0, 20);



  
  //___________________________________________________________________________

  //argumetns that are sent in: clustering radius, matching deltaR, etc
  int nEv=atol(argv[2]);  
  double  R = .3 + .1*double( atof( argv[1]) );
  double delta_R = double( atof(   argv[2]) );

  //class instances
  fjClustering *clustering         = new fjClustering(fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best); 
  JetMatching  *jetcuts            = new JetMatching();
  MyTopEvent   *reconstruction     = new MyTopEvent();

  
  if(debug) cout<<"There are "<<nEv<<" events! and R="<<R<<endl;
  
  // begin event loop 
  for (int iEv = 0; iEv < nEv; ++iEv) { 
    Event &event = pythia.event;

    //structure objets 
    struct first{
      vector<Particle> tops;
      vector<Particle> top;
      vector<Particle> antitop;
      vector<Particle> B;
      vector<Particle> W;
      vector<Particle> Wplus;
      vector<Particle> Wminus;
  
    };
    
    struct second{
      vector<Particle> nutrinos;
      vector<Particle> muons;
      vector<Particle> electrons;
    };

    struct minestruct{

      first  partons;
      second leptons;

    } data;

     //minestruct* date = new minestruct;

    //reset objects    
    clustering->ClearJets();
    reconstruction->Clear();
    
    if (!pythia.next()) continue;//generate events.skip if necessary
    if (debug && iEv==0) {pythia.info.list(); pythia.event.list();} 

    // storage for B Hadrons
    vector<Particle> Bhadrons;
    
    // storage for neutrinos, leptons, and _all_ final particles (parton or hadron level particles depending on mode)
    vector<Particle>  all_stbl_ptcls, MET; 
    

    //stable/final particles--------------------------------------------------------
    
    bool skipevent1(false),skipevent2(false),skip(false);
    for( size_t ithpt(0); ithpt < event.size(); ++ithpt){
      const Particle &particle = event[ithpt];
      const Particle &m1 = event[particle.mother1()];
      const Particle &m2 = event[particle.mother2()];
      const Particle &m3 = event[m1.mother1()];
      const Particle &m4 = event[m2.mother2()];


    //skip dilepton event
      if ( particle.isLepton() && (m1.id() == -24 || m2.id() == -24) &&  ( m3.id() == -6 || m4.id()== -6) )
       { skipevent1 = true; break; }

    //find tops
     if(particle.idAbs() == 6   &&   particle.daughter1()!= particle.daughter2()  &&  ( particle.daughter1()!=6 || particle.daughter2() != 6 ) ) {
       data.partons.tops.push_back(particle);
       if(particle.id()==6)
         data.partons.top.push_back(particle);
       if(particle.id()==-6)
         data.partons.antitop.push_back(particle);
     }

     if (isBhadron(particle.id())) 
       Bhadrons.push_back(particle);    

     if( isNu(particle) ) 
       MET.push_back(particle);

      //final particles
     if(! particle.isFinal()) continue;

      //Let's cluster particles that are not leptons from topw 
     if ( particle.isLepton() && (m1.id()== 24 || m2.id() ==24) &&  (m3.id()==6 || m4.id()==6)) {
       munus++;
       skipevent2 = true;
       if (skipevent1 && skipevent2) break;

       if(isTau( particle ) ){ skip = true; break;}

       if( isNu(particle))
         data.leptons.nutrinos.push_back(particle);

       if( isMu(particle))
         data.leptons.muons.push_back(particle);

       if( isElectron(particle) )
         data.leptons.electrons.push_back(particle);

       continue;
     }

      //final particle are stored in the jetcluster object for clustering later
     clustering->push_back(particle);

   }

    //skip if all hadronic decay or tau from wplus
   if( skipevent1 && skipevent2 ){ cout<<"dileptonic decay from t-tbar"<<endl;  continue; }
   if( skip ){ cout<<"we found a tau! "<<endl; continue;} 

   if ( (data.partons.top.size() + data.partons.antitop.size()) > 2) {
    pythia.event.list();
    fatal(Form("Event %d has too many good tops!?",iEv)); }


    //skip this event
    if(  data.leptons.electrons.size() == 0 &&  data.leptons.muons.size() == 0) {cout<<"\n no electrons and no  muons: All hadronic decay?!"<< endl; continue; }



    if( data.leptons.nutrinos.size() == 0 ){ cout<<"where are my nutrinos!"<<endl; continue;}



    //find w decay products-----------------------------------------------


    //here we look for the decay products of the W 
    //and check the decay channels
    std::pair <int,int> daughters;
    for ( size_t itop=0; itop < data.partons.tops.size(); ++itop) {
      const Particle& top = data.partons.tops[itop];

      // We have a top! it should go to a W and a b (most of the time)
      int W_index = top.daughter1(), b_index = top.daughter2();
      if (event[W_index].idAbs()!= 24) { W_index = top.daughter2(); b_index = top.daughter1(); } 
      const Particle &Wboson = event[W_index];
      const Particle &bquark = event[b_index];

      if(bquark.idAbs()!=5 || Wboson.idAbs()!=24) { printf(" We have top -> %d %d\n",Wboson.id(),bquark.id()); continue;}

      //store particles
      data.partons.W.push_back(Wboson);
      data.partons.B.push_back(bquark);

      //check for decay channels here by looking for sisters
      //sisterList();

      if(debug) {
       TString tup = top.id()==6 ? " top" : " tbar";
       cout <<"daughters of"<<tup<< " are: "<<event[W_index].id() <<" and  " <<event[b_index].id()<<"\n"<<endl;
     }

   } 

   if (debug  && iEv < 10) {
    printf("\nEvent %d:\n",iEv);
    printf("  We found %d b-quarks\n",int(data.partons.B.size()));
    for (size_t i=0;i < data.partons.B.size();++i) PrintPtcl(data.partons.B[i],Form("b-quark %d",i));
      clustering->PrintJets(); }


    //cluster particles and then do matching--------------------------------------------------


    //cluster particles
    clustering->doClustering();

    vector<fastjet::PseudoJet> all_jets = clustering->GetJets();   
    vector<fastjet::PseudoJet> btags, lightjets;

    //fastjet::sorted_by_pt(matchedjets);

    bool skipoverlapremoval = false;

    if(skipoverlapremoval) goto step_2;

    if(data.leptons.muons.size()==0) goto step_1;
    jetcuts->SetParam(.4,2.4,3.0);
    all_jets= jetcuts->OverlapRemoval(data.leptons.muons, all_jets);

    step_1:
    if(data.leptons.electrons.size() == 0) goto step_2;
    jetcuts->SetParam(.4, 2.4, 3.0);
    all_jets= jetcuts->OverlapRemoval(data.leptons.electrons, all_jets);

    step_2:
      //find b jets/ delta R, etamax, ptmin  
    jetcuts->SetParam(.4, 2.4, 10.0); 
    btags = jetcuts->Match(data.partons.B, all_jets);
    lightjets= jetcuts->RemoveSubset(btags, all_jets);



    //pseudotop construction-------------------------------------------------------------- 



    //check if reuiremements are met
    if (! reconstruction->SelectedEvent(btags, lightjets, 30.0)) { continue;}

    if( debug && iEv < 20){
      cout<<"Bjets are: "<<endl;
      PrintPseudojets(btags);
      cout<<" Lightjets: \n ";
      PrintPseudojets(lightjets);
    }

    //call function that returns the pseudojet of the leptonic pseudotop
    //We send in bjet, non-bjets, nutrinos, muons, and electrons 
    fastjet::PseudoJet leptonpseudotop = reconstruction->Recon_Mass_Method_1(btags, lightjets, 
                                                                             data.leptons.nutrinos, data.leptons.muons,data.leptons.electrons );


    //W:s
    whad->Fill(           reconstruction->Returnhadronicw());
    wlep->Fill(           reconstruction->Returnleptonicw() );
    bestbtop->Fill(       reconstruction->Returnbestbtop() );
    withhad ->Fill(       reconstruction->Returnlastbtop() );

    //information about truth parton
    top_mass->Fill(        data.partons.top[0].m() );
    top_eta->Fill(         data.partons.top[0].eta() );
    top_y->Fill(           data.partons.top[0].y() );
    top_pt->Fill(          data.partons.top[0].pT() );

    //check is our pseudo lpeton top matches a tbar

    vector<bool> matches = reconstruction-> TopMatch( data.partons.tops, leptonpseudotop, hadronicpseudotop) ;
    bool tmatch = matches[0], tbarmatch=matches[1];


    // keep track on the matching efficiency
    // this is a profile (TPofile) - fill with (x,y), and it will keep track
    // of the y-mean and plot it in bins of x
    matchingEff_vs_topPt->Fill(      leptonpseudotop.pt(),  tmatch);
    matchingEff_vs_topmass->Fill(    leptonpseudotop.m(),   tmatch);
    matchingEff_vs_topeta->Fill(     leptonpseudotop.eta(), tmatch);

    matchingEff_vs_tbmass->Fill(  data.partons.antitop[0].m(), tbarmatch );
    tbar_matching->Fill( leptonpseudotop.m(),tbarmatch );

    
    //separate histograms of matched and unmatched pseudo top
    if(tmatch){
      pseudo_top_match_mass->Fill(   leptonpseudotop.m()   );
      pseudo_top_match_eta->Fill(    leptonpseudotop.eta() );
      pseudo_top_match_pt ->Fill(    leptonpseudotop.pt()  );
    }
    else if(tbarmatch) {
       pseudo_tbar_nomatch_mass->Fill(     leptonpseudotop.m());
       pseudo_tbar_nomatch_eta->Fill(      leptonpseudotop.m());
       pseudo_tbar_nomatch_pt->Fill(       leptonpseudotop.m());
     } 
     else 
     {
      nomatch_mass->Fill(      leptonpseudotop.m()   );
      nomatch_eta ->Fill(      leptonpseudotop.eta() );
      nomatch_pt  ->Fill(      leptonpseudotop.pt()  );
    }
    

  } //end of event loop




  //_____________________________________________________

  if(debug) cout<<" Times reached the filling of leptons: " <<munus<<endl;

  // Statistical summary
  if(debug) pythia.stat(); 

  //write to file and close it
  outFile->Write();
  outFile->Close();

  return 0;

  } //end main function 


//check for b hadrons
  bool isBhadron(int pdg) {
    if (pdg<0) pdg *=-1;
    if (pdg<500) return false;
    if (pdg/100%5 == 0 || pdg/1000==5) return true; 
    return false;
  }

  void PrintPseudojets( const vector<fastjet::PseudoJet> &myjets) {


    for(size_t n(0); n < myjets.size(); ++n) {
      printf(" (pT,eta,phi,m) = (%6.1f GeV,%5.2f,%5.2f,%5.1f GeV) : \n",
        myjets[n].pt(),
        myjets[n].eta(),
        myjets[n].phi(),
        myjets[n].m());	     
    }  
  }
