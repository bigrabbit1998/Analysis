#include "Partons.h"
#include "Strucpie.h"

using namespace Pythia8;

mystruct mypie;
fjClustering * clustering = new fjClustering();

void Search::SetClustering( double rad)
{
	radius = rad;

}

Search::Search()
{
	//constructor 
}


void Search::Clear()
{
	mypie.tops.clear();
	mypie.W.clear();
	mypie.B.clear();

	mypie.electrons.clear();
	mypie.muons.clear();
	mypie.neutrinos.clear();
}


/*
void PrintPtcl(const Particle &ptcl, TString comment="")
{
	printf("  (pT,eta,phi,m) = (%6.1f GeV,%5.2f,%5.2f,%5.1f GeV) pdgID %4d : %s\n",
	       ptcl.pT(),ptcl.eta(),ptcl.phi(),ptcl.m(),ptcl.id(),comment.mypie());
} 
*/

bool Search::GENERATE()
{
	
	Pythia pythia("",false);
	//pythia.settings.mode("tune:pp", 6);
	pythia.settings.resetMode("tune:pp");
	pythia.readFile("particlesettings.cmnd");
	pythia.init(2212,2212,7000);
	pythia.settings.flag("print:quiet",true);

	if (!pythia.next()) return false;
	
	Event &event = pythia.event;
	cout<<"we made an event" <<endl;

	for(size_t i(0); i < event.size() ; ++i)
	{

		const Particle & particle = event[i];
		const Particle & d1 = event[particle.daughter1()];
		const Particle & d2 = event[particle.daughter2()];


		if( d1.id() ==d2.id() ) continue;
		if(d1.idAbs() == 6 && d2.idAbs()== 6 ) continue;
		if(particle.idAbs() == 6)
		{
			mypie.tops.push_back(particle);
			if( d1.idAbs()!= 24 ) { mypie.B.push_back(d1);mypie.W.push_back( d2 ); }
			else if( d1.idAbs() != 5 ) { mypie.B.push_back(d2); mypie.W.push_back( d1); }
			else return false;
		}

		//if(bquark.idAbs()!=5 || Wboson.idAbs()!=24) { printf(" We have top -> %d %d\n",Wboson.id(),bquark.id()); continue;}
		//if( mypie.tops.size() == 2) return false;
		//{
		//	pythia.event.list();
		//	fatal(Form("Event %d has too many good tops!?",)); 
		//}
	}

	bool skipevent2(false), skipevent1(false), skip(false);
	for( size_t  i= 0; i < event.size(); ++i)
	{
		const Particle & particle = event[i];
		const Particle & m1 = event[particle.mother1()];
		const Particle & m2 = event[particle.mother2()];
		

		if(!particle.isFinal()) continue;


		if(isTau( particle) ) 
			{ skip = true; break; }

		if( isNu(particle) && (m2.idAbs() == 24 || m1.idAbs() ==24) ) 
			{ mypie.neutrinos.push_back(particle); continue; }

		if( isMu(particle) && (m2.idAbs() == 24 || m1.idAbs() ==24) )  
		{ 
			mypie.muons.push_back(particle); continue; 
			if(m2.id() == -24) skipevent2 = true;
			if(m1.id() == 24) skipevent1=true;
		}
		if( isElectron( particle) && (m2.idAbs() == 24 || m1.idAbs() ==24) )  
		{ 
			mypie.electrons.push_back(particle); continue; 
			if(m2.id() == -24) skipevent2 = true;
			if(m1.id() == 24) skipevent1=true;
		}

		clustering->push_back(particle,-1);
		

	}

	if( skipevent1 && skipevent2 ) 
	{
		cout<<"dileptonic decay from t-tbar \n "<<endl;  
		return false; 
	}
	if( skip )
	{
		cout<<"we found a tau!\n "<<endl; 
		return false;; 
	} 

	if(  mypie.electrons.size() ==0 && mypie.muons.size() == 0  ) { cout<<"all hadronic decay"<<endl; return false; }
	if(  mypie.neutrinos.size() == 0 ) { cout<<"no neutrinos"<<endl; return false; }

	clustering->doClustering();

	//if (Debug  && iEv < 10)
	//{
	//	printf("\nEvent %d:\n",iEv);
	//	printf("  We found %d b-quarks\n",int(mypie.B.size()));
	//	for (size_t i=0;i < mypie.B.size();++i) PrintPtcl(mypie.B[i],Form("b-quark %d",i) );
	//		clustering->PrintJets(); 
	//}

	clustering->doClustering();

	return true;

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
