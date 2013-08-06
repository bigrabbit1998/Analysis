#ifndef JETMATCHING__HH
#define JETMATCHING__HH 1

#include "TLorentzVector.h"
#include "fjClustering.h"
#include "MyEvent.h"

class JetMatching{

 public:
  
  typedef vector<fastjet::PseudoJet> Pseudovector;
  typedef vector<Pythia8::Particle> Particlevector;
  
  //constructor
  JetMatching();
  void Clear();
  

  //this removes Pseudovector that are close to input particles
  Pseudovector OverlapRemoval(const Particlevector&,const Pseudovector&);
  

  //check for event requirements
  bool SelectedEvent(int,int, int);


  //delta r
  double Return_DR( Pythia8::Particle &, Pythia8::Particle &);
  double Return_DR( fastjet::PseudoJet &, fastjet::PseudoJet &);

  // matching is done here
  void SetParam( double,double,double);  
  void PrintMatches();
  bool ComparePt(fastjet::PseudoJet, fastjet::PseudoJet );
  bool ChiSquare(TLorentzVector, TLorentzVector );
  

  //recieves ptcut, etacut, and jets to be cut
  Pseudovector Cuts(const Pseudovector&, double ptcut, double etacut );

  //matched a set of particle sot Pseudovector
  Pseudovector Match(const Particlevector &,const Pseudovector&);

  //this removes selected jets from 
  Pseudovector RemoveSubset(const Pseudovector &,const Pseudovector& );



 private:
  
  double m_DeltaR, m_ptmin, m_etamax;

  int m_size_of_bjets, m_size_of_lightjets;  
  
};

#endif
