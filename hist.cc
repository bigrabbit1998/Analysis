#include <iostream>
#include <fstream>
#include "root.h"
#include "Pythia.h"
#include "fjClustering.h"
#include "JetSeparation.h"
#include "MyTopEvent.h"


using namespace Pythia8;


//__________________________


#define _Do_RECON_true
#define _Do_HISTOGRAMS_true
#define _Do_MATCHING_true



//_____________________________



void PrintPseudojets( const vector<fastjet::PseudoJet> &);

bool isElectron(const Pythia8::Particle &ptcl) { return ptcl.idAbs()==11;} 

bool isTau(const Particle &ptcl) { return ptcl.idAbs() ==15;}

bool isMu(const Particle &ptcl) { return ptcl.idAbs()==13;}

bool isNu(const Particle &ptcl) { return ptcl.idAbs()==12 || ptcl.idAbs()==14 || ptcl.idAbs()==16; }

void fatal(TString msg) { printf("ERROR:\n\n Check your logic %s\n\n",msg.Data()); abort(); }


int main(int argc, char* argv[]) 
{
	
  if (argc != 4 || ( atol(argv[1]) < 1) )
    {
      printf("ERROR: %s did not recieve the right arguments \n",argv[0]);
      return 1;
    }
  
  bool Debug= false;
  bool excludeTau = true;

  int Nev=atol(argv[2]);  
  double  R = .1*double( atof( argv[1]) );

  // Create the ROOT application environment. 
  TApplication theApp("hist", &argc, argv);

  ofstream outfile;

  //Create file on which histograms will be saved.
  TString of_name = "ParticleLevel_clusteringR_."+string(argv[1])+"_tune_"+string(argv[3])+"_.root";
  TFile* outFile = new TFile("rootplots/"+of_name, "RECREATE");
  //TFile* outFile = new TFile(of_name, "RECREATE");

  Pythia pythia("",false);
  pythia.settings.resetMode("tune:pp");
  pythia.settings.mode("tune:pp", 3 ); 
  pythia.readFile("particlesettings.cmnd");
  pythia.init(2212,2212,7000);
  pythia.settings.flag("print:quiet",true);


  //_____________________________________________________________________________________


  vector<Pythia8::Particle> tops;
  vector<Pythia8::Particle> B;
  vector<Pythia8::Particle> W;

  vector<Pythia8::Particle> neutrinos;
  vector<Pythia8::Particle> muons;
  vector<Pythia8::Particle> electrons;

  vector<Pythia8::Particle> ETmiss;


  fjClustering *clustering = new fjClustering(fastjet::kt_algorithm, R, fastjet::E_scheme, fastjet::Best);
  JetMatching  *matching   = new JetMatching();
  MyTopEvent   *recons     = new MyTopEvent();


  //___________________________________________________________________--



  //not matched pseudos
  TH1D* pseudo_mass  = new TH1D("pseudo_mass"," Leptonic pseudotop mass ; mass [GeV]", 100,130,220.);
  TH1D* pseudo_eta   = new TH1D("pseudo_eta"," leptonic pseudotop eta; #eta", 100 ,-6, 6. );
  TH1D* pseudo_pt    = new TH1D("pseudo_pt"," leptonic pseudotop pT; #it{p}_{T} [GeV]", 100 ,0 ,300. );


  //TH1D* pseudo_tbar_match_mass  = new TH1D("pseudo_mass_tbar_matched"," tbarmatch leptonic pseudotop mass; mass [GeV]", 100,130,220.);
  //TH1D* pseudo_tbar_match_eta   = new TH1D("pseudo_eta_tbar_matched"," tbarmatch leptonic pseudotop eta; #eta", 100 ,-6,6.);
  //TH1D* pseudo_tbar_match_pt    = new TH1D("pseudo_pt_tbar_matched"," tbarmatch leptonic pseudotop pt;#it{p}_{T} [GeV]", 100 ,0,300.);
  
  //matched pseudos
  //TH1D* pseudo_top_match_mass  = new TH1D("pseudo_mass_matched"," tmatched leptonic pseudotop; mass [GeV]", 100 ,120,220.);
  //TH1D* pseudo_top_match_pt    = new TH1D("pseudo_pt_matched",  " tmatched leptonic pseudotop; #it{p}_{T} [GeV]", 100 ,0,300.);
  //TH1D* pseudo_top_match_eta   = new TH1D("pseudo_eta_matched","  tmatched leptonic pseudotop ; #eta ", 100 ,-6,6.);

    
  //hadronic and leptonic Ws
  TH1D* whad_mass = new TH1D("whadronic_mass","hadronic pseudo W; mass [GeV]", 100,30.,120.);
  TH1D* whad_pt   = new TH1D("whadronic_pt","  hadronic pseudo W; pT [GeV]", 100,0.,300.);
  TH1D* whad_eta  = new TH1D("whadronic_eta"," hadronic pseudo W; eta [GeV]", 100, -6.,6.);

  //leptonic W from definition 
  TH1D* wlep_eta   = new TH1D("wleptoniceta"," leptonic pseudo W; #eta ", 100,-6,6.);
  TH1D* wlep_pt    = new TH1D("wleptonicpt","   leptonic pseudo W; p_{T}", 100,0, 290.);
  TH1D* wlep_mass  = new TH1D("wleptonicmass","   leptonic pseudo W; mass [GeV]", 100,40, 90.);

  //TH1D* wlep_eta2   = new TH1D("2wleptoniceta"," leptonic pseudo W; #eta ", 100,-6,6.);
  //TH1D* wlep_pt2    = new TH1D("2wleptonicpt","   leptonic pseudo W; p_{T}", 100,0, 290.);
  //TH1D* wlep_mass2  = new TH1D("2wleptonicmass","   leptonic pseudo W; mass [GeV]", 100,40, 90.);

  TH1D* bjetmass   = new TH1D("bjetmass", " Bjet mass; Mass [GeV]", 100, 0, 10.);
  

  //delta R relationships
  TH2D* dr_vs_pseudomass           = new TH2D("core1"," delta R_vs mass of pseudo top; mass [GeV]; delta_R", 100 , 150, 190., 20, 0., .7);
  TH2D* delbbjets_vs_pseudotopmass = new TH2D("core3"," Delta_R b-bjets vs pseudotop mass ; Mass ", 100, 120, 200 , 100 , 0., 1.9 );
  TH2D* jethadw_vs_truehadw_mass   = new TH2D("had_hadt"," Mass hadronic: pseudo W vs true W ; PseudoJet; True W ", 100, 40., 90 , 100 , 0, 90. );


  //matching efficiency  
  TProfile*  matchingEff_vs_topPt    =  new TProfile("eff_vs_topPt"," Matching efficiency vs pT ; #it{p}_{T} [Gev]; efficiency",300,0,300);
  TProfile*  matchingEff_vs_topmass  =  new TProfile("eff_vs_topMass","Matching Efficiency vs mass_{top} ; mass [GeV]; efficieny ",100,130,210);
  TProfile*  matchingEff_vs_topeta   =  new TProfile("eff_vs_topeta","Matching efficiency vs #eta ; #eta ; efficiency",20,-5,5);


  //3D matching plots 
  //TH3D* del_R_relations_1  = new TH3D("deltas"," delR ratios; del_ratio ; del_partons; del_pseudos ", 100, 0., .6 , 100 , 0, .6 , 100, 0, .6);
  //TH3D* mass_relations_1   = new TH3D("Dmass1"," pseudo_top-leponticW-bjet_delR ; pseudotop ; leptonic W; bjet ", 100, 140., 200. , 100 ,50, 140. , 100, 0, 8.);



  //___________________________________________________________________________


  for (size_t iev(0); iev < Nev; ++iev)
    {

      if (!pythia.next()) continue;
      Event &event = pythia.event;

      tops.clear();
      B.clear();
      W.clear();
      neutrinos.clear();
      muons.clear();
      electrons.clear();
      ETmiss.clear();

      clustering->ClearJets();
      matching->Clear();
      recons->Clear();


      vector<std::pair< Particle, bool> > temporary,temporary2;
      std::pair<Particle,bool> temp2,temp;
      //true for hadronic, false for leptonic


      bool skipevent2(false), skipevent1(false), skip(false);
      int positionhadw(0), position1(0), position2(0);

      for(size_t i(0); i < event.size() ; ++i)
	{

	  const Particle & particle = event[i];
	  const Particle & d1 = event[particle.daughter1()];
	  const Particle & d2 = event[particle.daughter2()];
	  const Particle & m1 = event[particle.mother1()];
	  const Particle & m2 = event[particle.mother2()];

	  int M1 = particle.mother1();


	  if( particle.idAbs() == 6) 
	    {
	      if( d1.id() ==d2.id() ) continue;
	      if(d1.idAbs() == 6 || d2.idAbs() == 6 ) continue;

	      if(particle.idAbs() == 6) tops.push_back(particle);
	      if( d1.idAbs()!= 24 && d2.idAbs() != 5) 
		{
		  if(d1.idAbs() == 5 && d2.idAbs() == 24)
		    { B.push_back(d1); W.push_back(d2); }
		  else
		    { cout<<" top doesn't go to b-W"<<endl; break; }
		}
	      else 
		{ 
		  if ( d1.idAbs() ==24 && d2.idAbs() ==5 )
		    { B.push_back(d2); W.push_back(d1); }
		  else 
		    { cout<< " top dose't go to bW"<<endl; break;}
		}
	      if( tops.size() == 2 )
		{
		  for(int f(0); f < tops.size(); ++f)
		    {
		      const Particle & topnow = tops[f];
		      int D1 = tops[f].daughter1(), D2 = tops[f].daughter2();
		      for( int h(0); h < 30 ; ++h)
			{
			  const Particle &part1 = event[D1];
			  const Particle &part2 = event[D2];

			  if(part1.idAbs() ==24 ) { D1 = part1.daughter1(); D2 = part1.daughter2(); continue; }
			  if(part2.idAbs() ==24 ) { D1 = part2.daughter1(); D2 = part2.daughter2(); continue; }

			  if(part1.isLepton()) 
			    {

			      if( event[part1.mother1()].id() ==24)
				{ position1 = part1.mother1(); skipevent2 = true; }

			      if( event[part1.mother1()].id() == -24 )
				{ position2 = part1.mother1(); skipevent1 = true; }
			      temp2.first= event[part1.mother1()]; temp2.second =false; 
			      temp.first = topnow; temp.second = false;
			      temporary.push_back(temp2); 
			      temporary2.push_back(temp); break;

			    }

			  if(part1.isQuark() )
			    {
			      if( event[part1.mother1()].id() ==24) position1 = part1.mother1(); 
           
			      if( event[part1.mother1()].id() == -24 ) position2 = part1.mother1(); 
           
			      temp2.first =event[part1.mother1()] ; temp2.second = true ; 
			      temp.first = topnow; temp.second = true;
			      temporary.push_back(temp2);
			      temporary2.push_back(temp); 
			      break;
			    }

			}
		    }
		}

	      continue;
	    }


	  if( skipevent2 && skipevent1 ) break;

	  if(! particle.isFinal()) continue;


	  if(isTau( particle) && ( M1 == position1 || M1 == position2 ) ) 
	    { skip = true; break; }

	  if( isNu(particle) && ( M1 == position1 || M1 == position2) ) 
	    { neutrinos.push_back(particle); continue; }

	  if( isMu(particle) && ( M1 ==position2 || M1 == position1 ) )  
	    { muons.push_back(particle); continue; }

	  if( isElectron( particle) && ( M1 == position2 || M1 ==position1)  )  
	    { electrons.push_back (particle); continue; }

	  if( fabs( particle.eta() ) > 4.4 && !isNu(particle))
	    { ETmiss.push_back(particle); continue; } 
	  
	   
	  clustering->push_back(particle,-1); 
	}

      if( skipevent1 && skipevent2 ) { cout<<"dileptonic decay from t-tbar \n \n"<<endl;  continue;}
      if( tops.size() != 2) { cout <<"The number of tops is not 2! \n\n"<<endl; continue; }
      if( W.size() != tops.size() ) 
	{
	  if(Debug) 
	    { 
	      cout<<"size of Ws is: " << W.size() << endl;
	      cout <<"size of tops is: " << tops.size() << endl; 
	      event.list();
	    }
	  continue; 
	} 


      if( skip && excludeTau ) { cout<<"we found a tau! \n\n "<<endl; continue; }
      if(  electrons.size() ==0 && muons.size() == 0  ) {cout<<"all hadronic decay"<<endl;  continue; }
      if(  neutrinos.size() == 0 ) { cout<<"no neutrinos \n\n"<<endl; continue; }

      clustering->doClustering();

      if (Debug  && iev < 20)
	{
	  printf("\nEvent %d:\n",iev);
	  printf("  We found %d b-quarks\n",int(B.size())) ;
	  clustering->PrintJets(); 
	}


#ifdef _Do_MATCHING_true

      vector<fastjet::PseudoJet> btaggedjets, lightjets, all_jets; 
      all_jets = *( clustering->GetJets() ) ;

      //apply Pt/eta cuts on pseudojets
      //matching->Cuts( 30.5, 4.5, &all_jets);

      matching->Match_method_2(B, all_jets, & btaggedjets); 
      matching->RemoveSubset( btaggedjets, all_jets, &lightjets );


#endif 




#ifdef _Do_RECON_true


      if( !( matching->SelectedEvent(2,2,4) ) ) { continue; } 


      bool we_match(false), tbarbool(false), topbool(false);
      fastjet::PseudoJet leptonpseudotop, leptonicw, hadronicw, lastbtop, bestbop ;
      Particle nearestparton, usefulW, usefultop, nearestB;


      if(! temporary[0].second ) { usefulW = temporary[0].first;}
      else if( ! temporary[1].second) { usefulW = temporary[1].first;}	
      else { cout<<"something wrong with finding the hadronic w decay"<<endl; }

      if( ! temporary2[0].second) { usefultop = temporary2[0].first ;}
      else if( ! temporary2[1].second) { usefultop = temporary2[1].first ;}
      else { cout<<"something wrong with fiding the correct top lep" <<endl; }


      //this send in the leptonic top and leptonic w partons
      recons->Parton_leptonic( usefultop, usefulW);

      //Here we begin top reconstruction
      recons->Initialize_Reconstruction( btaggedjets, lightjets, neutrinos); 
      recons->Recon_Mass_Method_2( muons, electrons, &leptonpseudotop);
      recons->Closest_Match( tops, leptonpseudotop, &nearestparton); 

      if(nearestparton.id() ==  6 ) topbool = true;
      if(nearestparton.id() == -6 ) tbarbool = true;

      hadronicw = *(recons->Return_Recon_Particles( "HW"));
      leptonicw = *(recons->Return_Recon_Particles( "LW"));

      pair<fastjet::PseudoJet, fastjet::PseudoJet> B_Wpair = *( recons->Return_B_W() );


      recons->Closest_Match(B, B_Wpair.first, & nearestB);

      double del_partons_pseudotop = recons->Return_DR( nearestparton, leptonpseudotop); 
      double del_B; 
      del_B = recons->Return_DR( nearestB, B_Wpair.first);

      if( nearestparton.id() == usefultop.id() && ( fabs(leptonpseudotop.m() - usefultop.m()/usefultop.m()) < .15  )) we_match = true;



#endif




#ifdef _Do_HISTOGRAMS_true


      pseudo_mass->Fill( leptonpseudotop.m());
      pseudo_pt->Fill(    leptonpseudotop.pt());
      pseudo_eta->Fill( leptonpseudotop.eta());

      jethadw_vs_truehadw_mass->Fill(hadronicw.m(), usefulW.m() );  

      dr_vs_pseudomass->Fill( leptonpseudotop.m(), del_partons_pseudotop );

      bjetmass->Fill(B_Wpair.first.m());


      wlep_pt->Fill(      leptonicw.pt() );
      wlep_eta->Fill(     leptonicw.eta());
      wlep_mass ->Fill(   leptonicw.m() );

      //wlep_pt2->Fill(    B_Wpair.second.pt() );
      //wlep_eta2->Fill(    B_Wpair.second.eta() );
      //wlep_mass2->Fill(  B_Wpair.second.m() );

      whad_eta->Fill(  hadronicw.eta());
      whad_mass->Fill( hadronicw.m() );
      whad_pt->Fill(   hadronicw.pt()); 


      matchingEff_vs_topPt->Fill(      leptonpseudotop.pt(),  we_match);
      matchingEff_vs_topmass->Fill(    leptonpseudotop.m(),   we_match);
      matchingEff_vs_topeta->Fill(     leptonpseudotop.eta(), we_match);


#endif 


    }




  //meanfile.open("txt", ios::app | ios::out );
  //meanfile << std::setprecision(3) << mean << " " << std::setprecision(3) << selectedratio<<endl;
  //meanfile.close(); 
  // 


  //write to file and close it
  outFile->Write();
  outFile->Close();

  return 0;

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

//check for b hadrons
bool ISBHADRON(int pdg) 
{
  if (pdg<0) pdg *=-1;
  if (pdg<500) return false;

  int q3 = (pdg/10)%10;
  if( ( pdg/1000%5==0 && pdg/100%5 == 0 && q3>0 && q3<5) || pdg/1000%5==0 ) return true ;
  return false;
}

