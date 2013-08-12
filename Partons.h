#ifndef SEARCH__HH
#define SEARCH_HH 1

#include <iostream>

#include "Pythia.h"

#include "root.h"

//fastjet 
#include "fjClustering.h"

//Jets
#include "JetSeparation.h"
#include "MyTopEvent.h"

//input/output
#include <fstream>

class Search
{
public:
	Search();
	~Search();

	
	bool GENERATE();
	bool ISHADRON( int ); 
	void SetClustering( double );
	void Clear();

	bool isElectron(const Pythia8::Particle &ptcl) { return ptcl.idAbs()==11;} 

	bool isTau(const Pythia8::Particle &ptcl) { return ptcl.idAbs() ==15;}

	bool isMu(const Pythia8::Particle &ptcl) { return ptcl.idAbs()==13;}

	bool isNu(const Pythia8::Particle &ptcl) { return ptcl.idAbs()==12 || ptcl.idAbs()==14 || ptcl.idAbs()==16; }

	void fatal(TString msg) { printf("ERROR:\n\n  %s\n\n",msg.Data()); abort(); }

	double radius;

	/*first Copy( first copy) 
	{
		copy->neutrinos = data.neutrinos;
		copy->muons = data.muons;
		copy->electrons = data.electrons;

		copy->tops = data.tops;
		copy->B = data.B;
		copy->W = data.W;
	}*/


};

#endif 