//       Stdlib header file for input and output.
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
void PrintPtcl(const Particle &ptcl, TString comment="")
{
	printf("  (pT,eta,phi,m) = (%6.1f GeV,%5.2f,%5.2f,%5.1f GeV) pdgID %4d : %s\n",
	       ptcl.pT(),ptcl.eta(),ptcl.phi(),ptcl.m(),ptcl.id(),comment.Data());
} 
void PrintPseudojets( const vector<fastjet::PseudoJet> &);

//main function
int  main(int argc, char* argv[]) 
{


	bool partonMode=true;
	bool debug= false;


  // Create the ROOT application environment. 
	TApplication theApp("hist", &argc, argv);

  //generator without banners
	Pythia pythia("",false);

  //read settings from command file
  // if(partonMode){pythia.readFile("partonMode.cmnd");} else { pythia.readFile("particleMode.cmnd");}

  //pythia.initTuneEE(2);
	pythia.readString("SoftQCD:minBias = off");
	pythia.readString("SoftQCD:singleDiffractive = off");
	pythia.readString("SoftQCD:doubleDiffractive = off");
	pythia.readString("HardQCD:all = off");
	pythia.readString("PartonLevel:FSR = off");
	pythia.readString("PartonLevel:ISR = on ");
	pythia.readString("HadronLevel:all = on");
	pythia.readString("Top:qqbar2ttbar = on");
	pythia.readString("Top:gg2ttbar = on");
	if (!partonMode) 
	{
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
	TH1D* top_mass  = new TH1D("m_top"," truth mass ; mass [GeV] ",20,130.,190.);
	TH1D* top_pt    = new TH1D("pt_top","truth pT ;#it{p}_{T} [GeV];Frequency", 50, 0.  , 300.);
	TH1D* top_eta   = new TH1D("eta_top","truth #eta ;#eta^{top}", 20, -5., 5.);
	TH1D* top_y     = new TH1D("y_top","truth rapidity ;#it{y}^{top}", 20, -5., 5.);

  //not matched pseudos
	TH1D* pseudo_tbar_match_mass  = new TH1D("pseudo_mass_tbar_matched"," tbarmatch leptonic pseudotop mass; mass [GeV]", 100,130,190.);
	TH1D* pseudo_tbar_match_eta   = new TH1D("pseudo_eta_tbar_matched"," tbarmatch leptonic pseudotop eta; #eta", 100 ,-6,6.);
	TH1D* pseudo_tbar_match_pt    = new TH1D("pseudo_pt_tbar_matched"," tbarmatch leptonic pseudotop pt;#it{p}_{T} [GeV]", 100 ,0,300.);

  //matched pseudos
	TH1D* pseudo_top_match_mass  = new TH1D("pseudo_mass_matched"," tmatched leptonic pseudotop; mass [GeV]", 100 ,110,200.);
	TH1D* pseudo_top_match_pt    = new TH1D("pseudo_pt_matched",  " tmatched leptonic pseudotop; #it{p}_{T} [GeV]", 100 ,0,300.);
	TH1D* pseudo_top_match_eta   = new TH1D("pseudo_eta_matched","  tmatched leptonic pseudotop ; #eta ", 100 ,-6,6.);


  //different Bjet-W pairs
	TH1D* leptonicw_bestb_pair        = new TH1D("pseudo_btop_mass","hadronic W pseudotop mass; mass [GeV]", 100 ,100,210.);
	TH1D* hadronicw_bestb_pair        = new TH1D("pseudo_HW_top"," pseudotop mass using had w and best b jet; mass [GeV]", 100 ,100,210.);

  //hadronic and leptonic Ws
	TH1D* whad = new TH1D("whadronic"," hadronic pseudo W; mass [GeV]", 100,0.,200.);

	TH1D* wlep_eta = new TH1D("wleptoniceta"," leptonic pseudo W; #eta ", 100,-6,6.);
	TH1D* wlep_pt = new TH1D("wleptonicpt","   leptonic pseudo W; #{p}_{T}", 100,0, 100.);

	TH1D* wparton_pt = new TH1D("whadrosdfnic"," true lep W pt ;   pt [GeV]", 20,0.,250.);
	TH1D* wparton_eta = new TH1D("whadsc"," true lep W eta ; mass [GeV]", 20,-6.,6.);

	TH2D* dr_vs_pseudomass   = new TH2D("cores"," delata R vs mass of pseudo top; mass [GeV]; delta_R", 100 , 130, 190., 20, 0., .7);
	TH2D* dr_vs_wlep_mass    = new TH2D("cosdres"," mass of w vs mass of pseudo W ;mass [GeV]; delta_R", 100 ,50, 90 , 50, 0., .7);
	TH2D* drpartons_vs_drpseudos = new TH2D("deltarspseudo"," Delta_R ; del_recon; del_partons ", 100, 0., 7 , 50 , 0., 7. );

  //matching efficiency plots  
	TProfile*  matchingEff_vs_topPt   =  new TProfile("eff_vs_topPt"," Matching efficiency vs pT #it{p}_{T} ; p_{T} [Gev]; efficiency",50,0,300);
	TProfile*  matchingEff_vs_topmass =  new TProfile("eff_vs_topMass","Matching Efficiency vs mass_{top} ; mass [GeV]; efficieny ",20,130,190);
	TProfile*  matchingEff_vs_topeta  =  new TProfile("eff_vs_topeta","Matching efficiency vs #eta ; #eta ; efficiency",20,-5,5);
	TProfile*  matchingEff_vs_bmass   =  new TProfile("eff_vs_bmass", "Matching efficiency vs mass_{b}; mass [GeV]; efficiency", 5, 0, 20);


  //3D matching plots 
	TH3D* deltas = new TH3D("deltars"," delR ratios; del_ratio ; del_partons; del_pseudos ", 10, 0., .6 , 100 , 140, 190. , 100, 60, 90.);



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


  //----------------------------------------------------------


  //things for event logging
	int nCandidateEvent=0, nSelectedevents=0;

  // begin event loop 
	std::pair<Pythia8::Particle, double> *checkmatch = new std::pair<Pythia8::Particle, double>;
	for (int iEv = 0; iEv < nEv; ++iEv)
	{ 
		Event &event = pythia.event;

      //structure objets. move these outside of loop later. 
		struct first
		{
			vector<Particle> tops;
			vector<Particle> top;
			vector<Particle> antitop;
			vector<Particle> B;
			vector<Particle> W;
			vector<Particle> Wplus;
			vector<Particle> Wminus;
		};

		struct second
		{
			vector<Particle> neutrinos;
			vector<Particle> muons;
			vector<Particle> electrons;
		};

		struct minestruct
		{
			first  partons;
			second leptons;

		} data;


      //reset objects    
		clustering->ClearJets();
		reconstruction->Clear(); 
		jetcuts->Clear();

		if(debug) cout<<"clear objects reached\n"<<endl;

      if (!pythia.next()) continue;//generate events.skip if necessary
      if (debug && iEv==0) {pythia.info.list(); pythia.event.list();} 
      
      // storage for B Hadrons
      vector<Particle> stble_Bhadrons;

      // storage for neutrinos, leptons, and _all_ final particles (parton or hadron level particles depending on mode)
      vector<Particle>  all_stbl_ptcls, MET; 

      if(debug) cout<<"break pie \n "<<endl;
      
      //stable/final particles--------------------------------------------------------

      
      
      bool skipevent1(false),skipevent2(false),skip(false);
      for( size_t ithpt(0); ithpt < event.size(); ++ithpt)
      {
      	const Particle &particle = event[ithpt];
      	const Particle &m1 = event[particle.mother1()];
      	const Particle &m2 = event[particle.mother2()];
      	const Particle &m3 = event[m1.mother1()];
      	const Particle &m4 = event[m2.mother2()];
      	


      	//check for dileptonic event and skip
      	if ( particle.isLepton() && (m1.idAbs() == 24 || m2.idAbs() == 24) &&  ( m3.idAbs() == 6 || m4.idAbs()== 6) )
      	{
      		if(m3.id() == 6 || m4.id() == 6 )
      			skipevent1 = true; 
      		if(m3.id() == -6 || m4.id() == -6)
      			skipevent2 = true;
      	} 

      	if (skipevent1 && skipevent2) break;


      	if( fabs(particle.eta()) < 4.5 ) 
      		MET.push_back(particle);

	  //find tops
      	if(particle.idAbs() == 6   &&  particle.daughter1()!= particle.daughter2()  &&  ( particle.daughter1()!=6 || particle.daughter2() != 6 ) ) 
      	{

      		data.partons.tops.push_back(particle);

      		if(particle.id()==6) 
      			data.partons.top.push_back(particle);

      		if(particle.id()==-6) 
      			data.partons.antitop.push_back(particle);
      	}


	  //final particles
      	if(! particle.isFinal()) continue;

      	//Let's get the stable Bhadrons. 
      	if (isBhadron(particle.id())) 
      		stble_Bhadrons.push_back(particle);    

	  //Let's get the particles used for reconstruction
      	if ( particle.isLepton() && (m1.idAbs() == 24 || m2.idAbs() == 24) &&  (m3.idAbs() == 6 || m4.idAbs() == 6)) 
      	{
      		if(isTau( particle ) ) { skip = true; break;}

      		if( isNu(particle))
      			data.leptons.neutrinos.push_back(particle);

      		if( isMu(particle))
      			data.leptons.muons.push_back(particle);

      		if( isElectron(particle) )
      			data.leptons.electrons.push_back(particle);
      		continue;
      	}

	  //final particle are stored in the jetcluster object for clustering 
      	clustering->push_back(particle);
      }

      if(debug) cout<<"break pizza \n"<<endl;

      //skip if dileptonic decay or tau from top's w
      if( skipevent1 && skipevent2 ){ cout<<"dileptonic decay from t-tbar \n "<<endl;  continue; }
      if( skip ){ cout<<"we found a tau!\n "<<endl; continue;} 

      if ( (data.partons.top.size() + data.partons.antitop.size()) > 2) 
      {
      	pythia.event.list();
      	fatal(Form("Event %d has too many good tops!?",iEv)); 
      }


      //skip this event
      if(  data.leptons.electrons.size() ==0 && data.leptons.muons.size() == 0  )  
      	{ cout<<"\n no electrons and no  muons: All hadronic decay?!"<< endl; continue; }

      if( data.leptons.neutrinos.size() == 0 )
      	{ cout<<"where are my neutrinos!"<<endl; continue; }



      //find top decay products-----------------------------------------------



      //here we look for the decay products of the W 
      //and  check the decay channels
      for ( size_t itop=0; itop < data.partons.tops.size(); ++itop)
      {
      	const Particle &top = data.partons.tops[itop];

	  // We have a top! it should go to a W and a b (most of the time)
      	int W_index = top.daughter1(), b_index = top.daughter2();
      	if (event[W_index].idAbs()!= 24) { W_index = top.daughter2(); b_index = top.daughter1(); } 
      	const Particle &Wboson = event[W_index];
      	const Particle &bquark = event[b_index];

      	if(bquark.idAbs()!=5 || Wboson.idAbs()!=24) { printf(" We have top -> %d %d\n",Wboson.id(),bquark.id()); continue;}

	  //store particles
      	data.partons.W.push_back(Wboson);
      	data.partons.B.push_back(bquark);

      	if(Wboson.id() == 24) 
      		data.partons.Wplus.push_back(Wboson);

	  //check for decay channels here by looking for sisters
	  //sisterList();

      	if(debug) 
      	{
      		TString tup = top.id()==6 ? " top" : " tbar";
      		cout <<"daughters of"<<tup<< " are: "<<event[W_index].id() <<" and  " <<event[b_index].id()<<"\n"<<endl; 
      	}

      } 


      if (debug  && iEv < 10)
      {
      	printf("\nEvent %d:\n",iEv);
      	printf("  We found %d b-quarks\n",int(data.partons.B.size()));
      	for (size_t i=0;i < data.partons.B.size();++i) PrintPtcl(data.partons.B[i],Form("b-quark %d",i) );
      		clustering->PrintJets(); 
      }




      //cluster particles and then do matching--------------------------------------------------




      //cluster particles
      clustering->doClustering();

      vector<fastjet::PseudoJet> all_jets = *(clustering->GetJets() );   
      vector<fastjet::PseudoJet> btagjets, lightjets;


      bool skipoverlapremoval = true;
      if(skipoverlapremoval) goto step_2;

      if(data.leptons.muons.size()==0) goto step_1;
      jetcuts->SetParam(.4,2.4,3.0);
      all_jets= jetcuts->OverlapRemoval(data.leptons.muons, all_jets);

      step_1:
      if(data.leptons.electrons.size() == 0) goto step_2;
      jetcuts->SetParam(.4, 2.4, 3.0);
      all_jets= jetcuts->OverlapRemoval(data.leptons.electrons, all_jets);

      step_2: 
      jetcuts->SetParam(.4, 2.4, 10.0); 
      jetcuts->Match(data.partons.B, all_jets, &btagjets);
      jetcuts->RemoveSubset(btagjets, all_jets, &lightjets);



      //pseudotop construction-------------------------------------------------------------- 




      //keep track of  to number of candiate pseudotops
      nCandidateEvent++;

      if(debug) cout<<"almost there!"<<endl;

      //check if reuiremements are met ( size of bjets, size of lightjets, size of all jets, MET )
      if (! jetcuts->SelectedEvent(2, 2, 4)) 
      {
      	nSelectedevents++;
      	continue;
      }


      if( debug && iEv < 20)
      {
      	cout<<"Bjets are: "<<endl;
      	PrintPseudojets(btagjets);
      	cout<<" Lightjets: \n ";
      	PrintPseudojets(lightjets);
      }



      //call function that returns the pseudojet of the leptonic pseudotop
      //We send in bjet, non-bjets, nutrinos, muons, and electrons 
      fastjet::PseudoJet leptonpseudotop = reconstruction->Recon_Mass_Method_1(btagjets, lightjets,  
                                                                               data.leptons.neutrinos, data.leptons.muons,data.leptons.electrons );

      //check is our pseudo lpeton top matches a tbar
      bool tmatch(false), tbarmatch(false);
      
      if(debug) cout<<"nearly there"<<endl;
      
      reconstruction->Closest_Match(data.partons.tops, leptonpseudotop, checkmatch);
      if(checkmatch->first.id() == 6 && checkmatch->second < .6) tmatch = true;
      if(checkmatch->first.id() == -6 && checkmatch->second < .6) tbarmatch = true;

      
      //get w to find scalling factor 
      std::pair<Pythia8::Particle, double> W_jj;
      fastjet::PseudoJet tempw = reconstruction->Return_W4vec();
      reconstruction->Closest_Match(data.partons.W, tempw , &W_jj);

      double tempass= reconstruction->Returnleptonicw() ; 
      double m_w_jj = checkmatch->first.m();
      double wjj_w = tempass/ m_w_jj;
      
      if(debug) cout<<"hang in there"<<endl;

      //W:s
      whad->Fill(           reconstruction->Returnhadronicw());
      //wlep->Fill(           reconstruction->Returnleptonicw() );
      //bestbtop->Fill(       reconstruction->Returnbestbtop() );
      //withhad ->Fill(       reconstruction->Returnlastbtop() );
      fastjet::PseudoJet lpw = reconstruction->Return_W4vec();
      wparton_pt->Fill(   data.partons.Wplus[0].pT() );
      wparton_eta->Fill(  data.partons.Wplus[0].eta() );
      wlep_pt->Fill(      lpw.pt() );
      wlep_eta->Fill(     lpw.eta());

      //retrieve the pseudojets for the b amd w jets
      std::pair<fastjet::PseudoJet, fastjet::PseudoJet> mypair =  reconstruction->Return_delta();

      //delta r of w and b
      pair<Particle, double> read, write;
      double del1_r = jetcuts->Return_DR(mypair.first, mypair.second);
      reconstruction->Closest_Match( data.partons.W, mypair.second, &write );
      reconstruction->Closest_Match( data.partons.B, mypair.first, & read);
      double del2_r = jetcuts->Return_DR( read.first, write.first);

      if(debug) cout<<"you finshed!"<<endl;



	  //delta R space      
      if(tmatch)
      	dr_vs_pseudomass->Fill( leptonpseudotop.m(), checkmatch->second );
      dr_vs_wlep_mass->Fill(lpw.m(), checkmatch->second);
      deltas->Fill(write.second, leptonpseudotop.m(), read.second );
      drpartons_vs_drpseudos->Fill(del1_r, del2_r);



      // keep track on the matching efficiency
      // this is a profile (TPofile) - fill with (x,y), and it will keep track
      // of the y-mean and plot it in bins of x
      matchingEff_vs_topPt->Fill(      leptonpseudotop.pt(),  tmatch);
      matchingEff_vs_topmass->Fill(    leptonpseudotop.m(),   tmatch);
      matchingEff_vs_topeta->Fill(     leptonpseudotop.eta(), tmatch);


      //Let's fill some histograms!
      if(tmatch)
      {
      	pseudo_top_match_mass->Fill(   leptonpseudotop.m()   );
      	pseudo_top_match_eta->Fill(    leptonpseudotop.eta() );
      	pseudo_top_match_pt ->Fill(    leptonpseudotop.pt()  );
      }
      else if(tbarmatch) 
      {
      	pseudo_tbar_match_mass->Fill(     leptonpseudotop.m());
      	pseudo_tbar_match_eta->Fill(      leptonpseudotop.m());
      	pseudo_tbar_match_pt->Fill(       leptonpseudotop.m());
      } 

      //information about truth parton
      if(tmatch || tbarmatch)
      {
      	top_mass->Fill(        data.partons.top[0].m() );
      	top_eta->Fill(         data.partons.top[0].eta() );
      	top_y->Fill(           data.partons.top[0].y() );
      	top_pt->Fill(          data.partons.top[0].pT() );
      }




    } //end of event loop





  //_____________________________________________________




  // Statistical summary
    if(debug) pythia.stat(); 

  //write to file and close it
    outFile->Write();
    outFile->Close();

    return 0;

} //end main function 



//check for b hadrons
bool isBhadron(int pdg) 
{
	if (pdg<0) pdg *=-1;
	if (pdg<500) return false;

	int q3 = (pdg/10)%10;
	if( ( pdg/1000%5==0 && pdg/100%5 == 0 && q3>0 && q3<5) || pdg/1000%5==0 ) return true ;
	return false;
}


void PrintPseudojets( const vector<fastjet::PseudoJet> &myjets) 
{
	for(size_t n(0); n < myjets.size(); ++n)
	{
		printf(" (pT,eta,phi,m) = (%6.1f GeV,%5.2f,%5.2f,%5.1f GeV) : \n",
		       myjets[n].pt(),
		       myjets[n].eta(),
		       myjets[n].phi(),
		       myjets[n].m());	     
	}  
}	

