#include "PartonTop.h"
#include "MyEvent.h"
#include "fjClustering.h"



//main constructor
PartonTop::PartonTop(Pythia8::Particle parton) {
  //constructor called
}


bool PartonTop::IsTop() {  return  parton.id() ==6;}

bool PartonTop::IsAntiTop() {  return m_parton.id()==-6;}

bool PartonTop::IsGluon(){ return m_parton.id()==21;}
