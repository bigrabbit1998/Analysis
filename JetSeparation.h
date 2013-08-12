#ifndef JETMATCHING__HH 
#define JETMATCHING__HH 1

#include "TLorentzVector.h"
#include <iostream> 
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>

class JetMatching{

 public:
  
  typedef std::vector<fastjet::PseudoJet> Pseudovector;
  typedef std::vector<Pythia8::Particle> Particlevector;
  
  //constructor
  JetMatching();
  
  //delta r
  double Return_DR(const Pythia8::Particle &,  const Pythia8::Particle &);
  double Return_DR(const fastjet::PseudoJet &, const fastjet::PseudoJet &);
  double Return_DR(const Pythia8::Particle &,  const fastjet::PseudoJet &);

  // matching is done here
  void Clear();
  void Closest_Match( const Particlevector & , const fastjet::PseudoJet & , Pythia8::Particle * );
  void Closest_Match(const Pythia8::Particle &, const Pseudovector &, fastjet::PseudoJet * );
  void cuts( double ptcuts, double etacuts, Pseudovector * sendback);
  void Match_method_1(const Particlevector & , const Pseudovector &, Pseudovector * );
  void Match_method_2(const Particlevector & , const Pseudovector &, Pseudovector *);
  //void Match_method_3(const Particlevector & , const Pseudovector &, std::vector<std::pair<Pythia8::Particle, fastjet::PseudoJet> *);
  void OverlapRemoval(const Particlevector &,  const Pseudovector &, Pseudovector * );
  void PrintMatches();
  void RemoveSubset(const Pseudovector& ,const Pseudovector& , Pseudovector *);
  void SetParam( double,double,double);  
  

  bool ChiSquare(TLorentzVector, TLorentzVector ); 
  bool ComparePt(fastjet::PseudoJet, fastjet::PseudoJet );
  bool SelectedEvent(int,int, int);
      

 private:
  
  double m_DeltaR, m_ptmin, m_etamax;

  int m_size_of_bjets, m_size_of_lightjets;  
  
};

#endif
