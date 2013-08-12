//       Stdlib header file for input and output.
#include <iostream>

// Header file to access Pythia 8 program elements.
#include "Pythia.h"
#include "RECON.h"

//Jets
#include "JetSeparation.h"
#include "MyTopEvent.h"
#include "Strucpie.h"

using namespace Pythia8;

Reconstruct::Reconstruct()
{
	cout<<"constructor called"<<endl;
}


mystruct data; 


bool Reconstruct::RECON( int ReconMethod, double btagDR  )
{
/*
	jetcuts->Clear();
	vector<fastjet::PseudoJet> all_jets = (*clustering->GetJets());
	vector<fastjet::PseudoJet> btagjets, lightjets, trimmedjets, temp;

	bool skipoverlapremoval = true;
	if(skipoverlapremoval) goto step_2;

	if(data.muons.size()==0) goto step_1;
	jetcuts->SetParam(.4,2.4,3.0);
	jetcuts->OverlapRemoval(data.muons, all_jets,& trimmedjets );
	all_jets = trimmedjets;

	step_1:
	if(data.electrons.size() == 0) goto step_2;
	jetcuts->SetParam(.4, 2.4, 3.0);
	jetcuts->OverlapRemoval(data.electrons, all_jets, &trimmedjets);

	step_2: 
	jetcuts->SetParam(.4, 2.5, 10.0); 
	jetcuts->Match_method_1(data.B, all_jets, &btagjets);
	jetcuts->RemoveSubset(btagjets, all_jets, &lightjets);


	if(btagjets.size() < 2) cout<<"seems to be working"<<endl;

	//if (! jetcuts->SelectedEvent(2, 2, 4) ) { cout<<"bad event"<<endl; continue; } 

	fastjet::PseudoJet leptonpseudotop = reconstruction->Recon_Mass_Method_1(btagjets, lightjets,  
	                                                                         data.neutrinos, data.muons,data.electrons );
	
	pair<Particle, double> checkmatch;
	reconstruction->Closest_Match(data.tops, leptonpseudotop, &checkmatch);
	
    //retrieve the pseudojets for the b amd w jets
	std::pair<fastjet::PseudoJet, fastjet::PseudoJet> mypair =  *(reconstruction->Return_B_W());

    //delta r of w and b
	pair<Particle, double> wpair, bpair;

    //get the closest match w
	reconstruction->Closest_Match( data.W, mypair.second, &wpair );

    //get the closest match b
	reconstruction->Closest_Match( data.B, mypair.first, &bpair);

    //get the delta_r of two pseudojets, and two partons
	double del1_r = jetcuts->Return_DR( mypair.first, mypair.second);
	double del2_r = jetcuts->Return_DR( bpair.first, wpair.first);
*/	
}