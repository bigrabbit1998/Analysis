#include "MyTopEvent.h"
#include "fjClustering.h"
#include "MyEvent.h"
#include <utility> 
#include "JetSeparation.h"

typedef vector<fastjet::PseudoJet> mypseudojets;
typedef vector<Pythia8::Particle> myparticlejets;

//constructor 
MyTopEvent::MyTopEvent()
{
  cout<<"Default constructor called"<<endl;
}

//clear storage devices
void MyTopEvent::Clear() 
{ 
  cout<<"clear"<<std::endl;
  m_tvectors.clear();
  m_mass=0;
  m_masswlep=0;
  m_masswhad =0;
  m_massbestb = 0;
  m_lastbtop=0;
 
}


//this returns the delta r of a particle to a pseudojet
double MyTopEvent::Return_DR( Pythia8::Particle & a, fastjet::PseudoJet & b)
{
  double DR(0);
  DR = sqrt(  pow((a.eta() +- b.eta()),2) +  pow((a.phi() - b.phi() ),2)  );
  return DR;
}


//this returns the closest partons
void MyTopEvent::Closest_Match( myparticlejets & partons, fastjet::PseudoJet & pseudotop, 
                               std::pair<Pythia8::Particle, double> * booger)
{
  double temp_del( 99999 );
  int position(0);
  
  for(size_t i(0) ; i < partons.size(); ++i)
  {
    Pythia8::Particle &parton_i = partons[i];
    double del= Return_DR(parton_i, pseudotop);

    if(del < temp_del )
    {
     temp_del = del;
     position = i;
   }

 }

 booger->first  = partons[position];
 booger->second = temp_del;

}


// leptonic pseudo top from leptonic with its w and b jets 
fastjet::PseudoJet MyTopEvent::Recon_Mass_Method_1( mypseudojets & Bjets, mypseudojets &Lightjets,
                                                   myparticlejets & nu, myparticlejets & mu, myparticlejets & els )
{

  //why is this called!!!!
  if(Bjets.size() < 2 || Lightjets.size() < 2 ){ cout<< "something went terriblly wrong"; abort(); }  

  //hadronic w by two highest pt non-bjets
  fastjet::PseudoJet hadronicw =operator+(Lightjets[1],Lightjets[0]);
  
  vector<fastjet::PseudoJet> hadronw;
  hadronw.push_back(hadronicw);

  //find the best hadronic w and b jet pair
  int  bestcombos = BestPairs(Bjets, hadronw );

  mypseudojets remainingbjet;
  for(size_t i(0) ; i < Bjets.size(); ++i)
  {
    fastjet::PseudoJet &temp = Bjets[i];
    if(i != bestcombos )
     remainingbjet.push_back(temp);
 }

  //this forms the hadronic top 
 fastjet::PseudoJet hadronictop= operator+(Bjets[bestcombos],hadronicw);

  //make call to create leptonic W
 fastjet::PseudoJet leptonicVV = LeptonicW(nu, mu, els);
 m_leptonw = leptonicVV;

  //create leptonic pseudow
 fastjet::PseudoJet lighttop  = operator+(remainingbjet[0], leptonicVV);
 fastjet::PseudoJet bestb_top = operator+(Bjets[bestcombos], leptonicVV);


  //these contain the mass of different tops and w:s
 m_massbestb = bestb_top.m();
 m_masswlep = leptonicVV.m();
 m_masswhad =hadronicw.m();
 m_lastbtop = hadronictop.m();

 m_pair.first = remainingbjet[0];
 m_pair.second = leptonicVV ;


 return lighttop;
  //  return leptonicVV;
}




//---------------------------------------------------------------------




//method two for reconstructing the pseudotop
fastjet::PseudoJet MyTopEvent::Recon_Mass_Method_2(const mypseudojets & Bjets, const mypseudojets& Lightjets, 
                                                   const myparticlejets & ETmiss, const myparticlejets & mu, const myparticlejets & els )
{

  //why is this called!!!!
  fastjet::PseudoJet lightop;
  if(Bjets.size() < 2 || Lightjets.size() < 2 ){ cout<< "something went terriblly wrong"; abort(); }  
  
  //get the pz solutions
  fastjet::PseudoJet leptonicVV = LeptonicW(ETmiss, mu, els);
  m_leptonw = leptonicVV;

  vector<fastjet::PseudoJet> leptons;
  fastjet::PseudoJet lepton;
  if( mu.size() != 0 ) 
    { leptons.push_back(Summation(mu)); }
  else {  leptons.push_back(Summation(els)); }

  

  fastjet::PseudoJet one, two, three, final;
  double dif(99999), tempdif(0), mass(0);

  int position_b, position_nu, position_lep;
  //introduce uncertainty for chisquare
  //here we do the choosing of the pseudojets
  for( size_t i(0); i < m_tvectors.size(); ++i)
  {
    for(size_t m(0) ; m < Bjets.size(); ++m)
    {
     one = operator+(m_tvectors[i], Bjets[m]);
     for(size_t n(0); n < leptons.size(); ++n )
     {
       two = operator+(one, leptons[n]);
       mass = two.m();
       //eta = two.eta();
       //pt = two.pt();

       //deta = 
       tempdif = fabs(mass - m_topmass ) ;
       //dpt = fabs( )
       if(tempdif < dif ){dif = tempdif; final = two; position_lep = n; position_nu = i; position_b = m; }
     } 
   }

 }

 fastjet::PseudoJet lepw = operator+( leptons[position_lep], m_tvectors[position_nu] );
 fastjet::PseudoJet bjet = Bjets[position_b];
 std::pair< fastjet::PseudoJet, fastjet::PseudoJet> sendback (bjet, lepw);
 m_pair = sendback;

  //these contain the mass of different tops and w:s
  // m_massbestb = bestb_top.m();
  //  m_masswlep = leptonicVV.m();
  //m_masswhad =hadronicw.m();
  // m_lastbtop = hadronictop.m();

 return final;
}



double MyTopEvent::Returnhadronicw(){ return m_masswhad; }
double MyTopEvent::Returnleptonicw(){ return m_masswlep;}
double MyTopEvent::Returnbestbtop() { return m_massbestb; }
double MyTopEvent::Returnlastbtop() { return m_lastbtop;}
fastjet::PseudoJet MyTopEvent::Return_W4vec() { return m_leptonw;}





//-----------------------------------------------------------------




// Given one lepton and one neutrino
//  solve for pz of the neutrino by assuming e+nu make the W mass
//  return 4-vector of W for the solution with smaller nu |pz|
fastjet::PseudoJet MyTopEvent::LeptonicW(const myparticlejets & nutrino, const myparticlejets& muon, const myparticlejets& electron) 
{

  // make sure input makes sense
  if (muon.size()+electron.size()!=1) 
    fatal(Form("There must only be one lepton! You are giving me %d electrons and %d muons.",
          electron.size(),muon.size()));


  fastjet::PseudoJet met, mu, el, lepton;
  
  met = Summation(nutrino);

  mu = Summation(muon);
  el = Summation(electron);
  lepton = operator+(mu,el);

  double MET_phi = met.phi();
  double MET_et = met.eta();
  double M_W=80.3999;
  double M_W2 = M_W*M_W;

  vector<fastjet::PseudoJet> ret; //try an array of vector or arrays
  //  double nu_x = MET_et*cos(MET_phi),
  // nu_y = MET_et*sin(MET_phi);

  double nu_x = met.px(),
  nu_y = met. py();

  double M = .5*(M_W2-pow(lepton.m(),2)),
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
    fastjet::PseudoJet one(nu_x,nu_y,pz1,e1);
    fastjet::PseudoJet two(nu_x,nu_y, pz2, e2);
    ret.push_back( one);    
    ret.push_back(two);
  } else {
    double pz = B/(2*A),
    e = sqrt(nu_x*nu_x+nu_y*nu_y+pz*pz);
    fastjet::PseudoJet three(nu_x, nu_y, pz, e);
    ret.push_back(three);
  }

  m_tvectors = ret;
  double maxpz(0), tempdif(0), diff(99990), tempmass(0);
  //minimize the diff to the matched parton 
  for(size_t i(0); i < ret.size(); ++i) {
    double 
    n_x = ret[i].px(),
    n_y = ret[i].py(),
    n_z = ret[i].pz();

    double 
    pwx = lepton.px() + nu_x,
    pwy = lepton.py() + nu_y,
    pwz = lepton.pz() + n_z,   
    pwE = sqrt(M_W2 + pwx*pwx + pwy*pwy + pwz*pwz);

    diff = n_z;
    //tempmass = sqrt(pwE*pwE - pwx*pwx - pwy*pwy - pwz*pwz );
    
    //use the matched eta here
    //tempdif = fabs(M_W - tempmass); 
    
    //get the pz that minimizes the diff in eta and mass
    if(tempdif < diff)
    {
     diff= tempdif;
     maxpz = n_z;
   }

   
 }

  // double E = sqrt(nu_x*nu_x+nu_y*nu_y+maxpz*maxpz);

 double 
 pwx = lepton.px() + nu_x,
 pwy = lepton.py() + nu_y,
 pwz = lepton.pz() + maxpz,   
 pwE = sqrt(M_W2 + pwx*pwx + pwy*pwy + pwz*pwz);



  //for the leptonic pseudo w
 fastjet::PseudoJet leptonic(pwx,pwy,pwz,pwE);

 return leptonic;

}

// convert psuedojet to Pythia Vec4
Pythia8::Vec4 MyTopEvent::ConvertToVec4(const fastjet::PseudoJet& pj )
{
  double px= pj.px(), py=pj.py(), pz=pj.pz(), E=pj.e();
  return Pythia8::Vec4(px,py,pz,E);
}

//moves from psedojets to 4 vectors
vector<Pythia8::Vec4> MyTopEvent::ConvertToVec4(const mypseudojets & changetype )
{
  vector<Pythia8::Vec4> returnjets;
  for(size_t n(0); n < changetype.size(); ++n)
    returnjets.push_back(ConvertToVec4(changetype[n]));

  return returnjets;
}

fastjet::PseudoJet MyTopEvent::Summation(const myparticlejets &temporary)
{
  double px(0),py(0),pz(0),e(0);

  //sum over all leptons
  for(size_t v(0); v < temporary.size(); ++v)
  {
    px+=temporary[v].px();
    py+=temporary[v].py();
    pz+=temporary[v].pz();
    e+= temporary[v].e();
  }

  fastjet::PseudoJet sendback(px,py,pz,e);
  return sendback;
}



//here we take the particle considered to contribute to MET
void MyTopEvent::BuildMET(std::vector<Pythia8::Particle> &input_particles)
{

  fastjet::PseudoJet ETmiss = Summation( input_particles );

  //m_MET = ETmiss; 
};


//returns position of best top-w pair
int MyTopEvent::BestPairs( mypseudojets & btaggedjets,  mypseudojets& ws)
{
  //vector<std::pair<fastjet::PseudoJet, fastjet::PseudoJet> > *bestcombinations = new vector<std::pair<fastjet::PseudoJet, fastjet::PseudoJet> >;
  vector<std::pair<fastjet::PseudoJet, fastjet::PseudoJet> > bestpair;

  double dif(99999);
  double tempdif(0);
  int besti(0);
  std::pair<fastjet::PseudoJet, fastjet::PseudoJet> p;
  for(size_t a(0) ; a < btaggedjets.size() ; ++a)
  {
    const fastjet::PseudoJet &bj = btaggedjets[a];
    for (int i = 0; i < ws.size(); ++i)
    {
     fastjet::PseudoJet final = operator+(bj,ws[i]);
     tempdif=fabs(m_topmass - final.m());

     if(tempdif < dif) 
     { 
       dif=tempdif;
       besti=a;
     }
	  //p.first = btaggedjets[i];
	  //p.second =  besti;
	  //bestpair.push_back(p);

   }
 }


 return besti;
}





