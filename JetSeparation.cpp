#include "fjClustering.h"
#include "MyEvent.h"
#include "JetSeparation.h"
using namespace std;
 typedef vector<fastjet::PseudoJet> Pseudovector;
 typedef vector<Pythia8::Particle> Particlevector;

//default constructor
JetMatching::JetMatching()
{
  cout<<"constructor"<<endl;
}

//clear information of golbal variables
void JetMatching::Clear()
{
  m_size_of_bjets= 0;
  m_size_of_lightjets =0;

}



void JetMatching::SetParam(double DeltaR , double etamax, double ptmin) 
{
  //values set
  m_DeltaR= DeltaR; 
  m_etamax = etamax;
  m_ptmin = ptmin;

  printf(" deltaR = %4.2f, max eta = %4.2f, min pt = %4.2f",m_DeltaR, m_etamax, m_ptmin);  
}


//does matching
Pseudovector JetMatching::Match(const Particlevector &particles, const  Pseudovector &input_jets ) 
{

  //vector of matched jets and pairs
//  vector<std::pair<Pythia8::Particle, Pseudovector> > matchedpairs;
  Pseudovector matchedjets;

  for( size_t ijet=0 ; ijet < input_jets.size() ; ++ijet) 
  {
    const fastjet::PseudoJet &nth_jet = input_jets[ijet];

    bool match_j=false;
    for( size_t ith_p(0) ; ith_p < particles.size(); ++ith_p )
    {
      const Pythia8::Particle &pe = particles[ith_p];

      TLorentzVector p,j;
      p.SetPtEtaPhiM(pe.pT(),pe.eta(),pe.phi(),pe.m());
      j.SetPtEtaPhiM(nth_jet.pt(),nth_jet.eta(),nth_jet.phi(),nth_jet.m());
      
      if( p.Pt() > m_ptmin && fabs(p.Eta()) < m_etamax && p.DeltaR(j) < m_DeltaR )
        match_j=true;	

    }

    if(match_j)
     matchedjets.push_back(nth_jet);  
 }

 if(particles[0].idAbs() == 5)
  m_size_of_bjets = matchedjets.size();
else 
  m_size_of_lightjets = matchedjets.size();

return matchedjets;
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


//send in clustered jets and particles that should be removed from jets
Pseudovector JetMatching::OverlapRemoval(const Particlevector &input_particles, const Pseudovector &complete_jets)
{
  Pseudovector jets;

  //overlap removal
  for (size_t ijet=0; ijet < complete_jets.size();++ijet) 
  {
    const fastjet:: PseudoJet &jet = complete_jets[ijet];
    bool isCloseToParticle = false;  

    for (size_t i=0;i < input_particles.size(); ++i) 
    {

      const Pythia8::Particle &pe = input_particles[i];

      TLorentzVector p,j;
      p.SetPtEtaPhiM(pe.pT(),pe.eta(),pe.phi(),pe.m());
      j.SetPtEtaPhiM(jet.pt(),jet.eta(),jet.phi(),jet.m());

        //how many leptons should there be in the container?
      if ( p.Pt()>m_ptmin && fabs(p.Eta()) < m_etamax && p.DeltaR(j) < m_DeltaR )
       isCloseToParticle = true;
   }

   if (isCloseToParticle) continue; 
   jets.push_back(jet);


 }
 

 return jets;
}


//pt, eta cuts
Pseudovector Cuts ( const Pseudovector& input_jets, double ptcuts, double etacuts)
{
 Pseudovector cut_jets;

 for (size_t ijet=0; ijet < input_jets.size();++ijet) 
 {
  const fastjet:: PseudoJet &jet = input_jets[ijet];
  bool notremoved = false;    
  if( jet.pt() > ptcuts && fabs(jet.eta()) < etacuts )
    notremoved = true;

  if (notremoved) 
    cut_jets.push_back(jet);  
 }


return cut_jets;

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


// ( particles to match, new jets, all jets) so that jets and all jets are modified

//this will be not be needed later
vector<fastjet::PseudoJet> JetMatching::RemoveSubset(const vector<fastjet::PseudoJet>& subset,const vector<fastjet::PseudoJet> &set) 
{
  Pseudovector newset;

  for( size_t n(0); n < set.size(); ++n) 
  {
    const fastjet::PseudoJet  &s = set[n];

    for(size_t m(0); m < subset.size(); ++m)
    {
      const fastjet::PseudoJet &ss= subset[m];
      if(operator==(s,ss)) break;
      else if(m==(int(subset.size())-1))
       newset.push_back(s);
   }
 }

 return newset;
}
