#ifndef MYTOPEVENT__HH
#define MYTOPEVENT__HH 1

#include "TLorentzVector.h"
#include <utility> 
#include <vector>
#include <cstdio>
#include "Pythia.h"
#include <string>
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"


class MyTopEvent/*:: public TNamed*/
{

public:

    typedef std::vector<fastjet::PseudoJet> pseudovector;
    typedef std::vector<Pythia8::Particle> particlevector;


    MyTopEvent();


    void BestPairs(const pseudovector &  , const fastjet::PseudoJet &, fastjet::PseudoJet * );
    double Return_DR( const Pythia8::Particle & , const fastjet::PseudoJet &);
    
    void BuildMET(const particlevector &);
    void Clear();  
    void Closest_Match( const particlevector & , const fastjet::PseudoJet & , Pythia8::Particle *);
    void Ininitalize_Reconstruction(const pseudovector &, const pseudovector& ,const particlevector & , const particlevector & ) ;
    void LeptonicW(const particlevector&, const particlevector& , const particlevector&, fastjet::PseudoJet * ) ;
    void Neutrino_Pz_Solutions(const fastjet::PseudoJet &, const fastjet::PseudoJet &);
    void Recon_Mass_Method_1(const pseudovector &,const pseudovector&, const particlevector & nu, 
                             const particlevector & mu, const particlevector & els, fastjet::PseudoJet * );

    void Recon_Mass_Method_2(const pseudovector &, const pseudovector&, const particlevector & nu, 
                             const particlevector & mu, const particlevector & els, fastjet::PseudoJet * );




    std::pair<fastjet::PseudoJet, fastjet::PseudoJet>  *Return_B_W() 
    { 
        return &m_pair;
    };
    fastjet::PseudoJet Summation(const particlevector &);
    fastjet::PseudoJet *Return_Recon_Particles( std::string type )
    {
      if ( type == "HW")       return &m_hadronicw; 
      if ( type == "LW")       return &m_leptonicw;
      if ( type == "BST"  )    return &m_bestbtop;
      if ( type == "LBT"  )    return &m_lastbtop;

    }


    static const double m_topmass = 172.50, m_true_w = 80.5; 

    std::pair<fastjet::PseudoJet, fastjet::PseudoJet> m_pair;

    fastjet::PseudoJet m_leptonicw, m_hadronicw, m_MET, m_bestbtop, m_lastbtop;

    pseudovector m_bjets, m_lightjets, m_neutrinosolutions;

    particlevector m_muons, m_electrons, m_neutrinos;



private:

    void fatal(TString msg) { printf("ERROR:\n\n  %s\n\n",msg.Data()); abort(); }


};

#endif
