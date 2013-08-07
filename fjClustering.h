#ifndef FJCLUSTERING__HH
#define FJCLUSTERING__HH 1

#include <string>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector>
#include <cstdio>
#include "Pythia.h"

class fjClustering
{

 public:
  //constructors
  fjClustering();
  fjClustering(fastjet::JetAlgorithm, double, fastjet::RecombinationScheme, fastjet::Strategy);
  ~fjClustering();

  void ClearJets();
  void PrintJets();
  void doClustering();
  // void push_back(double px, double py, double pz, double E);
  void push_back(const Pythia8::Particle &part, int );
  void push_back(const Pythia8::Particle &part);

  std::vector<fastjet::PseudoJet> *GetJets() { return &outputJets;};

 private:
  
  std::vector<fastjet::PseudoJet> inputJets;
  std::vector<fastjet::PseudoJet> outputJets;

  fastjet::JetAlgorithm fjAlgorithm;
  fastjet::JetDefinition fjJetDefinition;
  double rParameter;
  fastjet::Strategy fjStrategy;
  fastjet::RecombinationScheme fjRecombScheme;
};

#endif
