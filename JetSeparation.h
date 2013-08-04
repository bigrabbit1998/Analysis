#ifndef JETMATCHING__HH
#define JETMATCHING__HH 1

#include "TLorentzVector.h"
#include "fjClustering.h"
#include "MyEvent.h"

class JetMatching{

 public:
  
  typedef vector<fastjet::PseudoJet> Pseudojets;
  typedef vector<Pythia8::Particle> Particlejets;
  
  //constructor
  JetMatching();
  
  /*Jets(void);
  Jets(const Jets &src);
  void operator=(const Jets &src);
  */

  //this removes pseudojets that are close to input particles
  Pseudojets OverlapRemoval(const Particlejets&,const Pseudojets&);
  
  // matching parameters
  void SetParam( double,double,double);
  void Clear();
  
  //matching is done here
  void PrintMatches();
  bool ComparePt(fastjet::PseudoJet, fastjet::PseudoJet );
  bool ChiSquare(TLorentzVector, TLorentzVector );
  bool SelectedEvent(const Pseudojets &, const Pseudojets&, double &);
  

  //recieves ptcut, etacut, and jets to be cut
  Pseudojets Cuts(const Pseudojets&, double ptcut, double etacut );

  //matched a set of particle sot pseudojets
  Pseudojets Match(const Particlejets&,const Pseudojets&);

  //this removes selected jets from 
  Pseudojets RemoveSubset(const Pseudojets &,const Pseudojets& );


 private:
  
  double m_DeltaR, m_ptmin, m_etamax;

  int m_size_of_bjets, m_size_of_lightjets;
  
  
};

#endif
