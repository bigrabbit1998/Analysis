#include "MyTopEvent.h"


typedef vector<fastjet::PseudoJet> pseudovector;
typedef vector<Pythia8::Particle> particlevector;

//constructor 
MyTopEvent::MyTopEvent()
{
  cout<<"Default constructor called"<<endl;
}

//clear storage devices
void MyTopEvent::Clear() 
{ 
  //clear pseudojets
  m_leptonw.reset(0,0,0,0);
  m_leptonicw.reset(0,0,0,0);
  m_hadronicw.reset(0,0,0,0);
  m_MET.reset(0,0,0,0);
  m_bestbtop.reset(0,0,0,0);
  m_lastbtop.reset(0,0,0,0);

  //clear vectors
  m_tvectors.clear();
  m_bjets.clear();
  m_lightjets.clear();

  //clear pairs
  m_pair.first.reset(0,0,0,0);
  m_pair.second.reset(0,0,0,0);
}


//this returns the delta r of a particle to a pseudojet
double MyTopEvent::Return_DR(const Pythia8::Particle & a,const fastjet::PseudoJet & b)
{
  double DR(0);
  DR = sqrt(  pow((a.eta() - b.eta()),2) +  pow((a.phi() - b.phi() ),2)  );
  return DR;
}


//this returns the closest partons
void MyTopEvent::Closest_Match(const particlevector & partons, const fastjet::PseudoJet & pseudotop, 
                               Pythia8::Particle * sendback)
{
  double del( 99 );
  int position(0);
  
  for(size_t i(0) ; i < partons.size(); ++i)
  {
    const Pythia8::Particle &parton_i = partons[i];
    double temp_del= Return_DR(parton_i, pseudotop);

    if(del > temp_del ) { del = temp_del;  position = i;  }

  }

  *sendback = partons[position];

}


// leptonic pseudo top from leptonic with its w and b jets 
void MyTopEvent::Recon_Mass_Method_1( const pseudovector & Bjets,const pseudovector &Lightjets, const particlevector & nu, 
                                     const particlevector & mu,const particlevector & els, fastjet::PseudoJet *sendback )
{

  //why is this called!!!!
  if(Bjets.size() < 2 || Lightjets.size() < 2 ){ cout<< "something went terriblly wrong"; abort(); }  

  //hadronic w by two highest pt non-bjets
  fastjet::PseudoJet hadronicw =operator+(Lightjets[1],Lightjets[0]);
  fastjet::PseudoJet remainingbjet; 

  //find the best hadronic w and b jet pair
  BestPairs( Bjets, hadronicw, &remainingbjet);

  //make call to create leptonic W
  fastjet::PseudoJet leptonicVV;
  LeptonicW(nu, mu, els, &leptonicVV);

  //create leptonic pseudow
  fastjet::PseudoJet lighttop  = operator+(remainingbjet, leptonicVV);
  
  m_leptonicw = leptonicVV;
  m_hadronicw = hadronicw;

  m_pair.first = remainingbjet;
  m_pair.second = leptonicVV ;

  *sendback = lighttop;
}



//method two for reconstructing the pseudotop
void MyTopEvent::Recon_Mass_Method_2(const pseudovector & Bjets, const pseudovector& Lightjets, const particlevector & ETmiss, 
                                     const particlevector & mu, const particlevector & els, fastjet::PseudoJet * sendback )
{

  //why is this called!!!!
  fastjet::PseudoJet lightop;
  if(Bjets.size() < 2 || Lightjets.size() < 2 ){ cout<< "something went terriblly wrong"; abort(); }  
  
  //make the hadronic pseudow from talk prescriptio
  m_hadronicw = operator+(Lightjets[0], Lightjets[1]);
  
  LeptonicW(ETmiss, mu, els, & m_leptonw);

  vector<fastjet::PseudoJet> leptons;
  fastjet::PseudoJet one, two, three, final;
  fastjet::PseudoJet lepton;

  
  if( mu.size() != 0 ) 
    { leptons.push_back(Summation(mu)); }
  else 
    {  leptons.push_back(Summation(els)); }
  

  double dif(999), tempdif(0), mass(0);
  int position_b, position_nu, position_lep;

  for( size_t i(0); i < m_tvectors.size(); ++i)
  {
    for(size_t m(0) ; m < Bjets.size(); ++m)
    {
     one = operator+(m_tvectors[i], Bjets[m]);
     for(size_t n(0); n < leptons.size(); ++n )
     {
       two = operator+(one, leptons[n]);
       mass = two.m();
       tempdif = fabs(mass - m_topmass ) ;
       if(tempdif < dif ) {dif = tempdif; final = two; position_lep = n; position_nu = i; position_b = m; }
     } 
   }

 }

 fastjet::PseudoJet lepw = operator+( leptons[position_lep],  m_tvectors[position_nu] );
 fastjet::PseudoJet bjet = Bjets[position_b];
 std::pair< fastjet::PseudoJet, fastjet::PseudoJet> send (bjet, lepw);
 m_pair = send;

 
 *sendback = final;
}



void MyTopEvent::Ininitalize_Reconstruction(const pseudovector & Bjets, const pseudovector& Lightjets,
                                            const particlevector & mu, const particlevector & els ) 
{

  //send in everything 

}



//this forms the leptonic w using one lepton and one nutrino from top decay
/*void MyTopEvent::LeptonicW(const particlevector &nutrino, const particlevector& muon, const particlevector& electron, 
                           fastjet::PseudoJet * sendback) 
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

  vector<fastjet::PseudoJet> ret; 
  double 
  nu_x = MET_et*cos(MET_phi),
  nu_y = MET_et*sin(MET_phi);

  //double 
  //nu_x = met.px(),
  //nu_y = met. py();

  double 
  M = .5*(M_W2-pow(lepton.m(),2)),
  lTnuT = (nu_x*lepton.px() + nu_y*lepton.py()),
  A = pow(lepton.e(),2) - pow(lepton.pz(),2),
  B = 2 * lepton.pz() * (M + lTnuT),
  C = pow(MET_et,2) * pow(lepton.e(),2) - pow(M,2) - pow(lTnuT,2) - 2*M*lTnuT,
  discr = B*B - 4*A*C;

  if (discr > 0.) 
  {
    double rtDiscr = sqrt(discr);
    double pz1 = (B+rtDiscr)/(2*A),
    pz2 = (B-rtDiscr)/(2*A),
    e1 = sqrt(nu_x*nu_x+nu_y*nu_y+pz1*pz1),
    e2 = sqrt(nu_x*nu_x+nu_y*nu_y+pz2*pz2);
    fastjet::PseudoJet one(nu_x,nu_y,pz1,e1);
    fastjet::PseudoJet two(nu_x,nu_y, pz2, e2);
    ret.push_back( one);    
    ret.push_back(two);
  } else 
  {
    double pz = B/(2*A),
    e = sqrt(nu_x*nu_x+nu_y*nu_y+pz*pz);
    fastjet::PseudoJet three(nu_x, nu_y, pz, e);
    ret.push_back(three);
  }

  m_tvectors = ret;
  double maxpz(0), tempdif(0), diff(0), tempmass(0);
  for(size_t i(0); i < ret.size(); ++i)
  {
    double n_z = ret[i].pz();

    tempdif = fabs(n_z);
    if(tempdif > diff)
    {     
      diff= tempdif;
      maxpz = n_z;
    }  
  }

  double 
  pwx = lepton.px() + nu_x,
  pwy = lepton.py() + nu_y,
  pwz = lepton.pz() + maxpz,   
  pwE = sqrt(M_W2 + pwx*pwx + pwy*pwy + pwz*pwz);


  //for the leptonic pseudo w
  sendback->reset(pwx,pwy,pwz,pwE);

}
*/


void MyTopEvent::Neutrino_Pz_Solutions( const fastjet::PseudoJet & lepton, const fastjet::PseudoJet &MET)
{

  vector<fastjet::PseudoJet> ret; 


  double MET_phi = m_MET.phi();
  double MET_et = m_MET.eta();
  double M_W=80.3999;
  double M_W2 = M_W*M_W;

  double 
  nu_x = MET_et*cos(MET_phi),
  nu_y = MET_et*sin(MET_phi);

  //double 
  //nu_x = met.px(),
  //nu_y = met. py();

  double 
  M = .5*(M_W2-pow(lepton.m(),2)),
  lTnuT = (nu_x*lepton.px() + nu_y*lepton.py()),
  A = pow(lepton.e(),2) - pow(lepton.pz(),2),
  B = 2 * lepton.pz() * (M + lTnuT),
  C = pow(MET_et,2) * pow(lepton.e(),2) - pow(M,2) - pow(lTnuT,2) - 2*M*lTnuT,
  discr = B*B - 4*A*C;

  if (discr > 0.) 
  {
    double rtDiscr = sqrt(discr);
    double pz1 = (B+rtDiscr)/(2*A),
    pz2 = (B-rtDiscr)/(2*A),
    e1 = sqrt(nu_x*nu_x+nu_y*nu_y+pz1*pz1),
    e2 = sqrt(nu_x*nu_x+nu_y*nu_y+pz2*pz2);
    fastjet::PseudoJet one(nu_x,nu_y,pz1,e1);
    fastjet::PseudoJet two(nu_x,nu_y, pz2, e2);
    ret.push_back( one);    
    ret.push_back(two);
  } else 
  {
    double pz = B/(2*A),
    e = sqrt(nu_x*nu_x+nu_y*nu_y+pz*pz);
    fastjet::PseudoJet three(nu_x, nu_y, pz, e);
    ret.push_back(three);
  }


  m_tvectors = ret;
  
}




fastjet::PseudoJet MyTopEvent::Summation(const particlevector &temporary)
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
/*void MyTopEvent::BuildMET(const std::vector<Pythia8::Particle> & input_particles)
{

  fastjet::PseudoJet ETmiss = Summation( input_particles );

  m_MET = ETmiss; 

};*/


//returns position of best top-w pair
  void MyTopEvent::BestPairs( const pseudovector & btaggedjets, const fastjet::PseudoJet & ws, fastjet::PseudoJet * sendback)
  {

    double dif(99999);
    double tempdif(0);
    int besti(0);
    
    for(size_t a(0) ; a < btaggedjets.size() ; ++a)
    {
      
      fastjet::PseudoJet final = operator+(btaggedjets[a],ws);
      tempdif=fabs(m_topmass - final.m());

      if(tempdif < dif)  { dif=tempdif; besti=a; }
    }

    if( besti == 0) *sendback = btaggedjets[1];
    else if( besti == 1 ) *sendback = btaggedjets[0];
    else *sendback = btaggedjets[0];


  } 





void MyTopEvent::LeptonicW(const particlevector &nutrino, const particlevector& muon, const particlevector& electron, 
                           fastjet::PseudoJet * sendback) 
{

  // make sure input makes sense
  if (muon.size()+electron.size()!=1) 
    fatal(Form("There must only be one lepton! You are giving me %d electrons and %d muons.",
          electron.size(),muon.size()));

  fastjet::PseudoJet met,mu,els, lepton;
  
  met = Summation(nutrino);
  mu = Summation(muon);
  els = Summation(electron);

  if( muon.size() !=0 ) lepton = mu;
  else if( electron.size() !=0) lepton = els;

  Neutrino_Pz_Solutions(lepton, met);

  double maxpz(0), tempdif(0), diff(0), tempmass(0), M_W2 = ( m_true_w * m_true_w ) ;
  
  for(size_t i(0); i < m_tvectors.size(); ++i)
  {
    double n_z = m_tvectors[i].pz();

    tempdif = fabs(n_z);
    if(tempdif > diff)
    {     
      diff= tempdif;
      maxpz = n_z;
    }  
  }

  double 
  pwx = lepton.px() + met.px(),
  pwy = lepton.py() + met.py(),
  pwz = lepton.pz() + maxpz,   
  pwE = sqrt(M_W2 + pwx*pwx + pwy*pwy + pwz*pwz);


  //for the leptonic pseudo w
  sendback->reset(pwx,pwy,pwz,pwE);

}


