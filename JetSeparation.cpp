#include "fjClustering.h"
#include "MyEvent.h"
#include "JetSeparation.h"
using namespace std;

//default constructor
JetMatching::JetMatching(){
  cout<<"constructor"<<endl;
}

void JetMatching::SetParam(double *DeltaR , double etamax, double ptmin) {
  //values set
  m_DeltaR= DeltaR[x][y]; 
  m_etamax = etamax;
  m_ptmin = ptmin;

  printf(" deltaR = %4.2f, max eta = %4.2f, min pt = %4.2f",m_DeltaR, m_etamax, m_ptmin);  
 }

 

//does matching
vector<fastjet::PseudoJet> JetMatching::Match(const vector<Pythia8::Particle> &particles, const  vector<fastjet::PseudoJet> &input_jets ) {

  vector<fastjet::PseudoJet> matchedjets;
  for( size_t ijet=0 ; ijet < input_jets.size() ; ++ijet) {
    const fastjet::PseudoJet &nth_jet = input_jets[ijet];

    bool match_j=false;
    for( size_t ith_p(0) ; ith_p < particles.size(); ++ith_p ){
      const  Pythia8::Particle &pe = particles[ith_p];

      TLorentzVector p,j;
      p.SetPtEtaPhiM(pe.pT(),pe.eta(),pe.phi(),pe.m());
      j.SetPtEtaPhiM(nth_jet.pt(),nth_jet.eta(),nth_jet.phi(),nth_jet.m());
      
      if( p.Pt() > m_ptmin && fabs(p.Eta()) < m_etamax && p.DeltaR(j) < m_DeltaR ){
       match_j=true;	
      }
    }
    
    if(match_j)
       matchedjets.push_back(nth_jet);  
  }

 
  return matchedjets;
}


//send in clustered jets and particles that should be removed from jets
vector<fastjet::PseudoJet> JetMatching::OverlapRemoval(const vector<Pythia8::Particle> &input_particles, const vector<fastjet::PseudoJet> &removal_jets){
  vector< fastjet::PseudoJet> jets;

  //overlap removal
  for (size_t ijet=0; ijet < removal_jets.size();++ijet) {
    const fastjet:: PseudoJet &jet = removal_jets[ijet];

    bool isCloseToParticle = false;    
    for (size_t i=0;i<input_particles.size();++i) {
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
vector<fastjet::PseudoJet> Cuts ( const vector<fastjet::PseudoJet> & input_jets, double ptcuts, double etacuts) {
 vector< fastjet::PseudoJet> cutjets;

 for (size_t ijet=0; ijet < input_jets.size();++ijet) {
  const fastjet:: PseudoJet &jet = input_jets[ijet];
  bool notremoved = false;    
  if( jet.pt() > ptcuts && fabs(jet.eta()) < etacuts )
    notremoved = true;

  if (notremoved) 
    cutjets.push_back(jet);
}
  
return cutjets;
};

bool JetMatching::ComparePt(fastjet::PseudoJet a, fastjet::PseudoJet b) {
  return a.pt() > b.pt();
}
 
bool JetMatching::ChiSquare(TLorentzVector parton, TLorentzVector jet ){
  //compare angles, energies,pt, and delta R
}

void JetMatching::PrintMatches() {
  std::cout <<"hello"<<std::endl;
}
  
vector<fastjet::PseudoJet> JetMatching::RemoveSubset(const vector<fastjet::PseudoJet>& subset,const vector<fastjet::PseudoJet> &set) {
  vector< fastjet::PseudoJet> newset;

  for( size_t n(0); n < set.size(); ++n) {
    const fastjet::PseudoJet  &s = set[n];

    for(size_t m(0); m < subset.size(); ++m){
      const fastjet::PseudoJet &ss= subset[m];
      
      if(operator==(s,ss)) break;
      else if(m==(int(subset.size())-1))
	      newset.push_back(s);
    }
  }
  return newset;
}
