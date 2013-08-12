#include "fjClustering.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"

void fjClustering::ChangeRParam( double R )
{
  fjJetDefinition = fastjet::JetDefinition(fastjet::kt_algorithm, R, fastjet::E_scheme, fastjet::Best);
  std::cout<<fjJetDefinition.description()<<std::endl;

}

fjClustering::fjClustering()
{
  std::cout<<" fjClustering called with default parameters "<<std::endl;
  double R = 1.0;
  *this = fjClustering(fastjet::kt_algorithm, R, fastjet::E_scheme, fastjet::Best);
}

fjClustering::fjClustering(fastjet::JetAlgorithm ja, 
			   double rparam,
			   fastjet::RecombinationScheme rs,
			   fastjet::Strategy s)
{
  rParameter = rparam;
  fjAlgorithm = ja;
  fjStrategy = s;
  fjRecombScheme = rs;
  fjJetDefinition = fastjet::JetDefinition(fjAlgorithm, rParameter, fjRecombScheme, fjStrategy);
  std::cout<<fjJetDefinition.description()<<std::endl;
}

fjClustering::~fjClustering()
{
  inputJets.clear();
  outputJets.clear();
 }

//we now clear the vector of jets
void fjClustering::ClearJets()
{
  outputJets.clear();
  //outputJets.resize(1, 0);
  inputJets.clear();
  //inputJets.resize(1, 0);
 }

//this function will print the current jets that have been clustered
void fjClustering::PrintJets()
{
  std::cout<<"Number of  output jets is: "<<outputJets.size()<<std::endl;
  for (int  i = 0; i < outputJets.size(); i++)
  {
    std::printf(" (pT ,eta ,phi ,m) = ( %lf,%lf, %lf, %lf, %d) \n\n",
		outputJets[i].pt(),
		outputJets[i].eta(),
		outputJets[i].phi(),
		outputJets[i].m(),
		outputJets[i].user_index());
  }
 }

//This is a function to compare the pT of two jets
bool ComparePt(fastjet::PseudoJet a, fastjet::PseudoJet b) {
  return a.pt() > b.pt();
}

//this does the clustering
void fjClustering::doClustering()
{
  fastjet::ClusterSequence cluster_seq(inputJets, fjJetDefinition);
  double pTcut=8.0;
  outputJets = cluster_seq.inclusive_jets(pTcut);
  std::sort(outputJets.begin(), outputJets.end(), ComparePt);
}

void fjClustering::push_back(const Pythia8::Particle &part, int pid)
{
  //fastjet::PseudoJet *jet= new fastjet::PseudoJet(px,py,pz,E);
  //inputJets.push_back(jet);
  //fastjet::PseudoJet *myjet2 = new fastjet::PseudoJet(part.px(),part.py(),part.pz(),part.e());
  fastjet::PseudoJet myjet3(part.px(),part.py(),part.pz(),part.e());
  myjet3.set_user_index (pid); 
  inputJets.push_back(myjet3);
  //event.push_back(fastjet::PseudoJet(px,py,pz,E));
  //delete myjet2;
  //  inputJets.push_back(fastjet::PseudoJet(px,py,pz,E));
  //inputJets.push_back(fastjet::PseudoJet(px,py,pz,E));
}

void fjClustering::push_back(const Pythia8::Particle &part)
{
  int pid=-1;
  push_back(part,pid);
    
}
