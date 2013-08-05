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
  
  
  /*Jets(void);
  Jets(const Jets &src);
  void operator=(const Jets &src);
  */

  //this removes Pseudovector that are close to input particles
  Pseudovector OverlapRemoval(const Particlevector&,const Pseudovector&);
  

  bool SelectedEvent(int,int, int);

  // matching parameters
  void SetParam( double,double,double);

  void Clear();
  
  //matching is done here
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
