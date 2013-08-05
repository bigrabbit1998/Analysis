  #include "fjClustering.h"
  #include "MyEvent.h"
  #include "MyTopEvent.h"
  #include "JetSeparation.h"
  //#include "PartonTop.h"

  typedef vector<fastjet::PseudoJet> mypseudojets;
  typedef vector<Pythia8::Particle> myparticlejets;

  //constructor 
MyTopEvent::MyTopEvent()
{
  cout<<"Default constructor called"<<endl;
}


void MyTopEvent::Clear() 
{ 
  cout<<"clear"<<std::endl;
}


double MyTopEvent::Return_DR( Pythia8::Particle & in_particle, fastjet::PseudoJet &in_pseudotop)
{
  double DR(0);
  
  DR = sqrt( (in_pseudotop.eta() - in_particle.eta() )*(in_pseudotop.eta() - in_particle.eta() ) + 
            (in_pseudotop.phi() - in_particle.phi() )*(in_pseudotop.phi() - in_particle.phi() ) );
  

  return DR;
}


  void MyTopEvent::TopsMatch_Closest( myparticlejets & partons, fastjet::PseudoJet & pseudotop, 
                                                                   std::pair<Pythia8::Particle, double> * booger)
{
  double temp_del( 99999);
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


  // leptonic pseudo top from leptonic 
fastjet::PseudoJet MyTopEvent::Recon_Mass_Method_1( mypseudojets & Bjets, mypseudojets &Lightjets,
                                                   myparticlejets & nu, myparticlejets & mu, myparticlejets & els )
{
  fastjet::PseudoJet errorjet(0,0,0,0);
  if(Lightjets.size() < 2 || Bjets.size() < 2) return errorjet ;

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


  fastjet::PseudoJet hadronictop= operator+(Bjets[bestcombos],hadronicw);

    //make call to create leptonic W
  fastjet::PseudoJet leptonicVV = LeptonicW(nu, mu, els);
  m_leptonw = leptonicVV;

  //

    //create leptonic pseudow
  fastjet::PseudoJet lighttop  = operator+(remainingbjet[0], leptonicVV);
  fastjet::PseudoJet bestb_top = operator+(Bjets[bestcombos], leptonicVV);

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
                                                    const myparticlejets & nu, const myparticlejets & mu, const myparticlejets & els )
  {
    cout<<"something"<<endl;
  }


  double MyTopEvent::Returnhadronicw(){ return m_masswhad; }
  double MyTopEvent::Returnleptonicw(){ return m_masswlep;}
  double MyTopEvent::Returnbestbtop() { return m_massbestb; }
  double MyTopEvent::Returnlastbtop() { return m_lastbtop;}
  fastjet::PseudoJet MyTopEvent::Return_W4vec() { return m_leptonw;}


  // Given one lepton and one neutrino
  //  solve for pz of the neutrino by assuming e+nu make the W mass
  //  return 4-vector of W for the solution with smaller nu |pz|
  fastjet::PseudoJet MyTopEvent::LeptonicW(const myparticlejets & nutrino, const myparticlejets& muon, const myparticlejets& electron) 
  {

   // make sure input makes sense
    if (muon.size()+electron.size()!=1) 
      fatal(Form("There must only be one lepton! You are giving me %d electrons and %d muons.",
        electron.size(),muon.size()));


    Pythia8::Vec4 met, mu, el, lepton;
  
    met = Summation(nutrino);
 
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

Pythia8::Vec4 MyTopEvent::Summation(const myparticlejets &temporary)
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

Pythia8::Vec4 sendback(px,py,pz,e);
return sendback;

}




void MyTopEvent::BuildMET(std::vector<Pythia8::Particle> &met)
{
  //check object type
  cout<<"something"<<endl;
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



