#include "fjClustering.h"
#include "MyEvent.h"
#include "MyTopEvent.h"
#include "JetSeparation.h"
//#include "PartonTop.h"

//constructor 
MyTopEvent::MyTopEvent(){
  cout<<"Default constructor called"<<endl;
}


void MyTopEvent::Clear() {

  cout<<"clear"<<std::endl;

}


double MyTopEvent::TopMatch( vector< Pythia8::Particle> & in_particle, fastjet::PseudoJet &in_pseudotop){
  double DR(0);
  for(size_t n(0) ; n < in_particle.size() ; ++ n ) {
    DR = sqrt( (in_pseudotop.eta() - in_particle[n].eta() )*(in_pseudotop.eta() - in_particle[n].eta() ) + 
      (in_pseudotop.phi() - in_particle[n].phi() )*(in_pseudotop.phi() - in_particle[n].phi() ) );
  }

  return DR;
}
  

/*
bool MYTopEvent:: AllHadronic( ) {

  return (m_t.hadronicDecay() && m_tbar.hadronicDecay());
}

bool MyTopEvent::DiLeponic() { 
  return (m_t.leptonicDecay() && m_tbar.leptonicDecay()); 
}

bool MyTopEvent::SemiLeptonic() {
  return !AllHadroic() && !DiLeptonic(); 
}

*/


// leptonic pseudo top from leptonic 
fastjet::PseudoJet MyTopEvent::Recon_Mass_Method_1(const mypseudojets & Bjets, const mypseudojets &Lightjets,
 vector<Pythia8::Particle> & nu, vector<Pythia8::Particle> & mu, vector<Pythia8::Particle> & els )
{
  fastjet::PseudoJet errorjet(0,0,0,0);
  if(Lightjets.size() < 2 || Bjets.size() < 2) return errorjet ;

  //hadronic w by two highest pt non-bjets
  fastjet::PseudoJet hadronicw =operator+(Lightjets[1],Lightjets[0]);

  //find the best hadronic w and b jet pair
  int bestb = BestCombination(Bjets, hadronicw);

  vector<fastjet::PseudoJet> remainingbjet;
  for(size_t k(0); k < Bjets.size(); ++k)
    if(k != bestb)
      remainingbjet.push_back(Bjets[k]);


    fastjet::PseudoJet hadronictop= operator+(Bjets[bestb],hadronicw);

    //make call to create leptonic W
    fastjet::PseudoJet leptonicVV = LeptonicW(nu, mu, els);

    //create leptonic pseudow
    fastjet::PseudoJet lighttop = operator+(remainingbjet[0], leptonicVV);
    fastjet::PseudoJet bestb_top = operator+(Bjets[bestb], leptonicVV);

    //these contain the mass of different tops and w:s
    m_massbestb = bestb_top.m();
    m_masswlep = leptonicVV.m();
    m_masswhad =hadronicw.m();
    m_lastbtop = hadronictop.m();
    // return heavytop;
    return lighttop;
  }


  //method two for reconstructing the pseudotop
  fastjet::PseudoJet MyTopEvent::Recon_Mass_Method_2(const mypseudojets &, const mypseudojets&, 
   vector<Pythia8::Particle> & nu, vector<Pythia8::Particle> & mu, vector<Pythia8::Particle> & els )
  {


  }


  double MyTopEvent::Returnhadronicw(){ return m_masswhad; }
  double MyTopEvent::Returnleptonicw() {return m_masswlep;}
  double MyTopEvent::Returnbestbtop(){ return m_massbestb; }
  double MyTopEvent::Returnlastbtop() { return m_lastbtop;}

// Given one lepton and one neutrino
//  solve for pz of the neutrino by assuming e+nu make the W mass
//  return 4-vector of W for the solution with smaller nu |pz|
  fastjet::PseudoJet MyTopEvent::LeptonicW(vector<Pythia8::Particle> & nutrino, vector<Pythia8::Particle> & muon, vector< Pythia8::Particle> &electron) {

  // make sure input makes sense
    if (muon.size()+electron.size()!=1) 
      fatal(Form("There must only be one lepton! You are giving me %d electrons and %d muons.",
        electron.size(),muon.size()));


    Pythia8::Vec4 met, mu, el, lepton;
  // let's make the MET as a sum of all nu:s
    met = Summation(nutrino);
  // fastjet::PseudoJet nu_sum(0,0,0,0);
  // for (size_t i=0;i < nutrino.size();++i) 
  //      nu_sum = operator+(nu_sum, nutrino[i]);

  // let's make the MET as a sum of all nu:s
  // Pythia8::Vec4 nu_sum4vec;
  //for (size_t i=0;i< nutrino.size();++i) 
  //nu_sum4vec += Pythia8::MakeVec(nutrino[i]);

    mu = Summation(muon);
    el = Summation(electron);
    lepton = operator+(mu,el);
    double MET_phi = met.phi();
    double MET_et = met.eta();
    double M_W=80.3999;
    double M_W2 = M_W*M_W;

  vector<TLorentzVector> ret; //try an array of vector or arrays
  double nu_x = MET_et*cos(MET_phi),
  nu_y = MET_et*sin(MET_phi),
  M = .5*(M_W2-pow(lepton.mCalc(),2)),
  lTnuT = (nu_x*lepton.px() + nu_y*lepton.py()),
  A = pow(lepton.e(),2) - pow(lepton.pz(),2),
  B = 2 * lepton.pz() * (M + lTnuT),
  C = pow(MET_et,2) * pow(lepton.e(),2) - pow(M,2) - pow(lTnuT,2) - 2*M*lTnuT,
  discr = B*B - 4*A*C;
  if (discr > 0.) {
    double rtDiscr = sqrt(discr);
    double pz1 = (B+rtDiscr)/(2*A),
    pz2 = (B-rtDiscr)/(2*A),
    e1 = sqrt(nu_x*nu_x+nu_y*nu_y+pz1*pz1),
    e2 = sqrt(nu_x*nu_x+nu_y*nu_y+pz2*pz2);
    ret.push_back(TLorentzVector(nu_x,nu_y,pz1,e1));
    ret.push_back(TLorentzVector(nu_x,nu_y,pz2,e2));
  } else {
    double pz = B/(2*A),
    e = sqrt(nu_x*nu_x+nu_y*nu_y+pz*pz);
    ret.push_back(TLorentzVector(nu_x,nu_y,pz,e));
  }


  double maxpz(0);
  for(size_t i(0); i < ret.size(); ++i) {
    double 
    n_x = ret[i].Px(),
    n_y = ret[i].Py(),
    n_z = ret[i].Pz(),
    n_e = ret[i].E();

    if(fabs( n_z ) > maxpz) 
      maxpz = n_z;  
  }

  double E = sqrt(nu_x*nu_x+nu_y*nu_y+maxpz*maxpz);

  double 
  pwx = lepton.px() + nu_x,
  pwy = lepton.py() + nu_y,
  pwz = lepton.pz() + maxpz,
  pwE = lepton.e() + E;

  //  double pE = sqrt( pwx*pwx + pwy*pwy + pwz*pwz); 
  /* 
  double A = ( met.px()* met.px() + met.py()* met.py() );
  double C = total.pz();
  double D = total.e(); 
  double gamma = ( A + D*D  - pwx*pwx - pwy*pwy - mw*mw - C*C)/2.0;
  double pz = ( gamma* D - sqrt( A*(C*D) - A*D*D*D*D + gamma*gamma*D*D) ) / ( C*C - D*D) ;
  */ 

  //use larger value in the sum of the z components
  // double pwz = (fabs( total.pz() + pz) >  fabs( total.pz() - pz ) ) ? (total.pz() + pz) : (total.pz() - pz) ;
  
  
  //double pwE = sqrt( mw*mw + pwz*pwx + pwy*pwy + pwz*pwz ) ;

  fastjet::PseudoJet leptonic(pwx,pwy,pwz,pwE);
  
  return leptonic;

}

// convert psuedojet to Pythia Vec4
Pythia8::Vec4 MyTopEvent::ConvertToVec4(const fastjet::PseudoJet& pj ){
  double px= pj.px(), py=pj.py(), pz=pj.pz(), E=pj.e();
  return Pythia8::Vec4(px,py,pz,E);
}

//moves from psedojets to 4 vectors
vector<Pythia8::Vec4> MyTopEvent::ConvertToVec4(const mypseudojets & changetype ){
  vector<Pythia8::Vec4> returnjets;
  for(size_t n(0); n < changetype.size(); ++n)
    returnjets.push_back(ConvertToVec4(changetype[n]));

  return returnjets;
}

Pythia8::Vec4 MyTopEvent::Summation(myparticlejets &temporary){
 double px(0),py(0),pz(0),e(0);

  //sum over all leptons
 for(size_t v(0); v < temporary.size(); ++v){
  px+=temporary[v].px();
  py+=temporary[v].py();
  pz+=temporary[v].pz();
  e+= temporary[v].e();
}

Pythia8::Vec4 sendback(px,py,pz,e);
return sendback;
}



void MyTopEvent::EventRequirements(double ptmin, double etamax, double ETmissmin, )

  //we need something to check if a jet is to be used. try different eta cuts
  //this decides if the event is to be used
bool MyTopEvent::SelectedEvent( ){

  //get information from global variables and
  //check for requirements like eta, pt, and etc
  /*
  lepton.size() !=1 
  return false ;
  eta < some vale
  
  pt < ptmin
  
  */

}

void MyTopEvent::BuildMET(objets & obj){
  //check object type


}

//returns position of best top-w pair
int MyTopEvent::BestCombination(const mypseudojets & btaggedjets, /* vector<Pythia8::Vec4> &w*/ fastjet::PseudoJet &in ){
mypseudojets bestbjets;

double tpmass(99999);
double tempdif(0);
int bpos(0);
double m_topmass=172.5;

  //choose the first two elements of each vector
  //i.e. the two highest pt jets
for(size_t a(0) ; a < btaggedjets.size() ; ++a){
  const fastjet::PseudoJet &bf = btaggedjets[a];

  fastjet::PseudoJet final = operator+(bf,in);
  tempdif=fabs(m_topmass - final.m());

  if(tempdif < tpmass){
    tpmass=tempdif;
    bpos=a;
  } 
}

return bpos;
}



