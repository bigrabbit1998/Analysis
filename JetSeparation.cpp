#include "fjClustering.h"
#include "JetSeparation.h"
using namespace std;
typedef vector<fastjet::PseudoJet> Pseudovector;
typedef vector<Pythia8::Particle> Particlevector;

bool debug(false);

//default constructor
JetMatching::JetMatching()
{
  cout<<"constructor called for JetMatching"<<endl;
}


//clear information of golbal variables
void JetMatching::Clear()
{
  m_size_of_bjets = 0;
  m_size_of_lightjets = 0;

}

//this function will return the delta r of two jets or two vectors
double JetMatching::Return_DR( const Pythia8::Particle & a, const Pythia8::Particle & b)
{
  double DR(0);
  DR = sqrt( pow((a.eta() - b.eta()), 2) + pow((a.phi() - b.phi()), 2 )  );
  return DR;
}

double JetMatching::Return_DR(const fastjet::PseudoJet & c,const fastjet::PseudoJet & d)
{
  double DR(0);
  DR = sqrt( pow( (c.eta() - d.eta()), 2) + pow( (c.phi() - d.phi()), 2 )  );
  return DR;
}

double JetMatching::Return_DR( const Pythia8::Particle & c,const fastjet::PseudoJet & d)
{
  double DR(0);
  DR = sqrt( pow( (c.eta() - d.eta()), 2) + pow( (c.phi() - d.phi()), 2 )  );
  return DR;
}


//this sets the parameters for the matching
void JetMatching::SetParam(double DeltaR , double etamax, double ptmin) 
{
  //values set
  m_DeltaR= DeltaR; 
  m_etamax = etamax;
  m_ptmin = ptmin;

  if (debug) printf(" deltaR = %4.2f, max eta = %4.2f, min pt = %4.2f !\n",m_DeltaR, m_etamax, m_ptmin);  
}



//match closet pseudo jet to a parton
void JetMatching::Match_method_2( const Particlevector & particles, const  Pseudovector &input_jets, Pseudovector *sendback )  
{
  
  for(size_t i(0); i < particles.size(); ++i)
  {
    fastjet::PseudoJet send; 
    Closest_Match(particles[i], input_jets, &send);

    sendback->push_back(send);
  }
  if(particles[0].idAbs() == 5) { m_size_of_bjets = sendback->size(); }
}


void JetMatching::Closest_Match(const Pythia8::Particle & parton, const Pseudovector &input_jets, fastjet::PseudoJet * sendback)
{
  double temp_del( 0 ), del(999);
  int position(0);
  for(size_t n(0); n < input_jets.size(); ++n)
  {
    double temp_del = Return_DR(parton, input_jets[n]);
    if( del > temp_del) 
      { position = n; del = temp_del;}
  }

  *sendback = input_jets[position];
}


void JetMatching::Closest_Match(const Particlevector & partons, const fastjet::PseudoJet & pseudotop, Pythia8::Particle * booger)
{
  double temp_del( 0 ), del(999);
  int position(0);

  for(size_t i(0) ; i < partons.size(); ++i)
  {
    const Pythia8::Particle &parton_i = partons[i];
    double temp_del= Return_DR(parton_i, pseudotop);

    if(del > temp_del ) { del = temp_del ;  position = i;  }
  }

  *booger = partons[position];

}


//does matching
void JetMatching::Match_method_1(const Particlevector &particles, const  Pseudovector &input_jets, Pseudovector *sendback ) 
{

  for( size_t ijet=0 ; ijet < input_jets.size() ; ++ijet) 
  {
    const fastjet::PseudoJet &nth_jet = input_jets[ijet];

    bool matched=false;
    for( size_t ith_p(0) ; ith_p < particles.size(); ++ith_p )
    {
     const Pythia8::Particle &pe = particles[ith_p];

     double tempd = Return_DR(pe, nth_jet);

     if( pe.pT() > m_ptmin && fabs(pe.eta()) < m_etamax && tempd < m_DeltaR ) matched=true; 
   }

   if(matched) sendback->push_back(nth_jet);  
 }

 if(particles[0].idAbs() == 5) { m_size_of_bjets = sendback->size(); }

}


//we need something to check if a jet is to be used. try different eta cuts
//this decides if the event is to be used
bool JetMatching::SelectedEvent( int size_b, int size_l, int size_all )
{

  if(m_size_of_bjets < size_b) return false;

  if( m_size_of_lightjets < size_l) return false;

  if(m_size_of_lightjets + m_size_of_bjets < size_all) return false;

  return true; 
}


//pt, eta cuts
void JetMatching::Cuts( double ptcuts, double etacuts, Pseudovector * sendback)
{

  Pseudovector temp;
  for (size_t ijet=0; ijet < sendback->size(); ++ijet) 
  {     
    if( (*sendback)[ijet].pt() > ptcuts && fabs((*sendback)[ijet].eta()) < etacuts ) temp.push_back( (*sendback)[ijet] );
  }

  sendback->clear();
  *sendback = temp;

};


bool JetMatching::ComparePt(fastjet::PseudoJet a, fastjet::PseudoJet b) 
{
  return a.pt() > b.pt();
}

bool JetMatching::ChiSquare(TLorentzVector parton, TLorentzVector jet )
{
  //compare angles, energies,pt, and delta R
}

void JetMatching::PrintMatches() 
{
  std::cout <<"hello"<<std::endl;
}


//this will be not be needed later
void JetMatching::RemoveSubset(const Pseudovector& subset,const Pseudovector& set, Pseudovector *newset) 
{
  newset->clear();

  for( size_t n(0); n < set.size(); ++n) 
  {
    const fastjet::PseudoJet  &s = set[n];

    for(size_t m(0); m < subset.size(); ++m)
    {
     const fastjet::PseudoJet &ss= subset[m];

     if(operator==(s,ss)) break;
     else if(m==(int(subset.size())-1))
       newset->push_back(s);
   }
 }

 m_size_of_lightjets = newset->size();
 
}


//send in clustered jets and particles that should be removed from jets
void JetMatching::OverlapRemoval(const Particlevector &input_particles, const Pseudovector &complete_jets, Pseudovector * sendback)
{
  sendback->clear();
  for (size_t ijet=0; ijet < complete_jets.size();++ijet) 
  {
    const fastjet:: PseudoJet &jet = complete_jets[ijet];

    bool isCloseToParticle = false;  

    for (size_t i=0; i < input_particles.size(); ++i) 
    { 

     const Pythia8::Particle &pe = input_particles[i];

     double dr = Return_DR( pe, jet);

     if ( pe.pT()>m_ptmin && fabs(pe.eta()) < m_etamax && dr < m_DeltaR ) isCloseToParticle = true; 
   }

   if (isCloseToParticle) continue; 

   sendback->push_back(jet);

 }

}



