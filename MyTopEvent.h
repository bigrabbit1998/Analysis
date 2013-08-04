#ifndef MYTOPEVENT__HH
#define MYTOPEVENT__HH 1

#include "fjClustering.h"
#include "MyEvent.h"
#include <utility> 


//this class will recieve all information of top event 
//and reconstruct the top mass. 
class MyTopEvent/*:: public TNamed*/
{

public:

	typedef vector<fastjet::PseudoJet> mypseudojets;
	typedef vector<Pythia8::Particle> myparticlejets;

  //constructor
	MyTopEvent();


	fastjet::PseudoJet Recon_Mass_Method_1(const mypseudojets &, const mypseudojets&, 
	                                       mypseudojets & nu, myparticlejets & mu, myparticlejets & els );

	fastjet::PseudoJet Recon_Mass_Method_2(const mypseudojets &, const mypseudojets&, 
	                                       myparticlejets & nu, myparticlejets & mu, myparticlejets & els );

	fastjet::PseudoJet LeptonicW(myparticlejets & , myparticlejets & , myparticlejets & );


	int BestCombination(const mypseudojets &, fastjet::PseudoJet & );

	double TopsMatch( vector<Pythia8::Particle> & , fastjet::PseudoJet& );



	double Returnhadronicw();
	double Returnleptonicw();
	double Returnbestbtop();
	double Returnlastbtop();

	void Clear(); 

	//check size of bjets, all jets, light jets, with met, 
	bool SelectedEvent(const myparticlejets &, const mypseudojets&, double &);

	vector<Pythia8::Vec4> ConvertToVec4(const mypseudojets &  );
	Pythia8::Vec4 Summation(myparticlejets &);
	Pythia8::Vec4 ConvertToVec4(const fastjet::PseudoJet& );



private:
	mypseudojets m_tops, m_bjets, m_lightjets, m_extra, m_tbars;
	myparticlejets m_muons, m_electrons, m_nutrinos;

	//mass of particles/pseudojets
	double m_mass, m_masswlep, m_masswhad, m_massbestb, m_lastbtop;
	double m_topmass = 172.50;

	void fatal(TString msg) { printf("ERROR:\n\n  %s\n\n",msg.Data()); abort(); }

};

#endif
