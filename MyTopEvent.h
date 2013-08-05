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


	fastjet::PseudoJet Recon_Mass_Method_1( mypseudojets &, mypseudojets&, 
	                                        myparticlejets & nu, myparticlejets & mu, myparticlejets & els );

	fastjet::PseudoJet Recon_Mass_Method_2(const mypseudojets &, const mypseudojets&, 
	                                       const myparticlejets & nu, const myparticlejets & mu, const myparticlejets & els );


	int  BestPairs(mypseudojets &, mypseudojets &);

	double Return_DR( Pythia8::Particle & , fastjet::PseudoJet &);

	void TopsMatch_Closest( myparticlejets & , fastjet::PseudoJet & , std::pair<Pythia8::Particle, double> *);

	double Returnhadronicw();
	double Returnleptonicw();
	double Returnbestbtop();
	double Returnlastbtop();

	void Clear();  

	void BuildMET(myparticlejets &);
	

	vector<Pythia8::Vec4> ConvertToVec4(const mypseudojets &  );
	Pythia8::Vec4 Summation(const myparticlejets &);

	Pythia8::Vec4 ConvertToVec4(const fastjet::PseudoJet& );

	fastjet::PseudoJet LeptonicW(const myparticlejets &, const myparticlejets& , const myparticlejets&) ;
	fastjet::PseudoJet Return_W4vec();

private:
	mypseudojets m_tops, m_bjets, m_lightjets, m_extra, m_tbars;
	myparticlejets m_muons, m_electrons, m_nutrinos;

	//mass of particles/pseudojets
	double m_mass, m_masswlep, m_masswhad, m_massbestb, m_lastbtop;
	static const double m_topmass = 172.50, m_true_w = 80.5; 
	fastjet::PseudoJet m_leptonw;

	void fatal(TString msg) { printf("ERROR:\n\n  %s\n\n",msg.Data()); abort(); }

};

#endif
