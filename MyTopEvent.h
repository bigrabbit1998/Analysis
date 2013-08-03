#ifndef MYTOPEVENT__HH
#define MYTOPEVENT__HH 1

#include "fjClustering.h"
#include "MyEvent.h"
#include <utility> 
//#include "PartonTop.h"

//this class will recieve all information of top event 
//and reconstruct the top mass. 
class MyTopEvent/*:: public TNamed*/
{

public:

	typedef vector<fastjet::PseudoJet> mypseudojets;
	typedef vector<Pythia8::Particle> myparticlejets;

  //constructor
	MyTopEvent();


	//Use parton information to determine decay channel
	//bool AllHadronic() { return (m_t.hadronicDecay() && m_tbar.hadronicDecay() ); }
	//bool DiLeponic() { return (m_t.leptonicDecay()&&m_tbar.leptonicDecay());  }
	//bool SemiLeptonic() { return !AllHadroic() && !DiLeptonic(); }



	fastjet::PseudoJet Recon_Mass_Method_1(const mypseudojets &, const mypseudojets&, 
		vector<Pythia8::Particle> & nu, vector<Pythia8::Particle> & mu, vector<Pythia8::Particle> & els );

	fastjet::PseudoJet Recon_Mass_Method_2(const mypseudojets &, const mypseudojets&, 
		vector<Pythia8::Particle> & nu, vector<Pythia8::Particle> & mu, vector<Pythia8::Particle> & els );

	fastjet::PseudoJet LeptonicW(myparticlejets & , myparticlejets & , myparticlejets & );


	int BestCombination(const mypseudojets &, fastjet::PseudoJet & );
	double TopMatch( vector<Pythia8::Particle> & , fastjet::PseudoJet& );
	double Returnhadronicw();
	double Returnleptonicw();
	double Returnbestbtop();
	double Returnlastbtop();
	void Clear(); 

	vector<Pythia8::Vec4> ConvertToVec4(const mypseudojets &  );
	Pythia8::Vec4 Summation(myparticlejets &);
	Pythia8::Vec4 ConvertToVec4(const fastjet::PseudoJet& );



private:
	mypseudojets m_tops, m_bjets, m_lightjets, m_extra, m_tbars;
	myparticlejets m_muons, m_electrons, m_nutrinos;

	double m_mass, m_masswlep, m_masswhad, m_massbestb, m_lastbtop;

	void fatal(TString msg) { printf("ERROR:\n\n  %s\n\n",msg.Data()); abort(); }

};

#endif
