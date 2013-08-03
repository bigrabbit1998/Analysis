#include "fjClustering.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"

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

void fjClustering::ClearJets()
{
  outputJets.clear();
  //outputJets.resize(1, 0);
  inputJets.clear();
  //inputJets.resize(1, 0);
 }

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

bool ComparePt(fastjet::PseudoJet a, fastjet::PseudoJet b) {
  return a.pt() > b.pt();
}

void fjClustering::doClustering()
{
  
  fastjet::ClusterSequence cluster_seq(inputJets, fjJetDefinition);
  double pTcut=8.0;
  outputJets = cluster_seq.inclusive_jets(pTcut);
  std::sort(outputJets.begin(), outputJets.end(), ComparePt);
  
}

/*void fjClustering::push_back(double px, double py, double pz, double E)
{
 //
  //fastjet::PseudoJet *jet= new fastjet::PseudoJet(px,py,pz,E);
  //inputJets.push_back(jet);
  inputJets.push_back(fastjet::PseudoJet(px,py,pz,E));

  }*/

void fjClustering::push_back(const Pythia8::Particle &part)
{
  //fastjet::PseudoJet *jet= new fastjet::PseudoJet(px,py,pz,E);
  //inputJets.push_back(jet);
  int pid=part.id();
  fastjet::PseudoJet *myjet2 = new fastjet::PseudoJet(part.px(),part.py(),part.pz(),part.e());
  fastjet::PseudoJet myjet3 = *myjet2; 
  myjet3.set_user_index (pid);  
  inputJets.push_back(myjet3);
  //event.push_back(fastjet::PseudoJet(px,py,pz,E));
  delete myjet2;
  //  inputJets.push_back(fastjet::PseudoJet(px,py,pz,E));
  //inputJets.push_back(fastjet::PseudoJet(px,py,pz,E));

}

