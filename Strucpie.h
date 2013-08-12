#ifndef DATASTRUCTPOO__HH
#define DATASTRUCTPOO__HH 1

#include "Pythia.h"
#include "fjClustering.h"
#include "MyTopEvent.h"


struct mystruct
{
	vector<Pythia8::Particle> tops;
	vector<Pythia8::Particle> B;
	vector<Pythia8::Particle> W;

	vector<Pythia8::Particle> neutrinos;
	vector<Pythia8::Particle> muons;
	vector<Pythia8::Particle> electrons;

} *data;

fjClustering *clustering    = new fjClustering();
JetMatching *jetcuts        = new JetMatching();
MyTopEvent  *reconstruction = new MyTopEvent();

	//lepton pesudo top, lepW, hadW, B_W pair
	//vector< fastjet::PseudoJet> reconstructed_partons;
	//vector< fastjet::PseudoJet> reconstructed_2_partons;


#endif