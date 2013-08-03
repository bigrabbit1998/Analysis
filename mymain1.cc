// Stdlib header file for input and output.
#include <iostream>

// Header file to access Pythia 8 program elements.
#include "Pythia.h"

// ROOT, for histogramming.
#include "TH1.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h""

bool isBHadron(int id) {
  // This snippet is meant to capture all B hadrons
  // as given in the PDG.
  if (id<0) id*=-1;
  if (id<500) return false;
  return (fmod(id/100,5.)==0.0 || id/1000==5);
}

using namespace Pythia8; 
int main() {
  // Initialize root
  TROOT root("Manuel's ROOT Session","PYTHIA Histograms");

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  
  // Uncomment line below to turns on all HardQCD processses 
  // These are 111-116 and  121-124
  //pythia.readString("HardQCD:all = on"); 
  
  // Turn on only bbar production:
  // g g    -> b bbar (subprocess 123)
  // q qbar -> b bbar (subprocess 124)
  pythia.readString("HardQCD:gg2bbbar = on");
  pythia.readString("HardQCD:qqbar2bbbar = on");
  
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  
  // Random number Generator Should be Set Here if needed (before pythia.init())
  // On seeds:
  // seed = -1 : default (any number < 0 will revert to the default).  seed = 19780503
  // seed = 0 : calls Stdlib time(0) to provide a seed based on the unix time
  // seed = 1 through 900 000 000: various numbers that can be used as seeds
  
  //pythia.readString("Random.setSeed = on");// doesn't work needs fixing
  //pythia.readString("Random.seed = 3000000");
 
  pythia.init( 2212, 2212, 200.);
  Hist mult("charged multiplicity", 100, -0.5, 799.5);
  
  TH1D* multHist = new TH1D("multHist","Multiplicity",100,-0.5,99.5);
  TH1D* bquarkPt = new TH1D("bquarkPt","bquarkPt",100,0,50);
  TH1D* bbarquarkPt = new TH1D("bbarquarkPt","bbar quark Pt",100,0,50);
  TH1D* B0mesonPt = new TH1D("BOmesonPt","B0mesonPt",100,0,50);
  TH1D* B0barmesonPt = new TH1D("BObarmesonPt","B0bar meson Pt",100,0,50);
  TH1D* electronFrombPt = new TH1D("electronFrombPt","electrons from b",100,0,30);
  TH1D* positronFrombPt = new TH1D("positronFrombPt","positrons from b",100,0,30);
  TH1D* epluseminusMinv = new TH1D("epluseminusMinv","e+ e- Inv. Mass",100,0,30);

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 10000; ++iEvent) {
    if (!pythia.next()) continue;
    if (iEvent < 1) {pythia.info.list(); pythia.event.list();} 
    // Find number of all final charged particles and fill histogram.
    // Find the b (id = 5) and bbar (id = -5), find their daughters, 
    // if daughters include electron (id = 11) and positron (id=-11), calculate their
    // invariant mass
    // Status flags:
    //   21 incoming particles of hardest subprocess
    //   23 outgoing particles of hardest subprocess
    //   81-89 primary hadrons produced by hadronization process (B mesons, e.g.)
    //   91-99 particles produced in decay process or by B-E effects (e.g. the electrons)

    int nCharged = 0;
    int indexBQuark(0), indexBbarQuark(0);
    for (int i = 0; i < pythia.event.size(); ++i) {
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) {
        ++nCharged; 
      }
      Particle& theParticle = pythia.event[i];
    
      if (theParticle.id() == 5 ) {
    indexBQuark = i;
    //cout << "Mother 1, Mother 2 = " << theParticle.mother1() << ", " << theParticle.mother2() << endl;
      }
      if (theParticle.id() == -5) {
    indexBbarQuark = i;
    //cout << "Mother 1, Mother 2 = " << theParticle.mother1() << ", " << theParticle.mother2() << endl;
      }
    } // particle loop

    cout << "Found b quark at index " << indexBQuark << endl;
    cout << "Found bbar quark at index " << indexBbarQuark << endl;
    bquarkPt->Fill(pythia.event[indexBQuark].pT());
    bbarquarkPt->Fill(pythia.event[indexBbarQuark].pT());
    mult.fill( nCharged );
    multHist->Fill(nCharged);
    //cout << "Event " << iEvent << ", Nch= " << nCharged << endl;
    
    
    //Find hadronization products of b and bbar.
    int bQuarkDaughter1 = pythia.event[indexBQuark].daughter1();
    int bQuarkDaughter2 = pythia.event[indexBQuark].daughter2();
    int bbarQuarkDaughter1 = pythia.event[indexBbarQuark].daughter1();
    int bbarQuarkDaughter2 = pythia.event[indexBbarQuark].daughter2();
    
    // Obtain the two hadrons from the fragmentation process
    // Use the PDG id's for this.  All B mesons id's are of the form xx5xx, and
    // all B baryons are of the form 5xxx.
    // So we obtain the id, (make it positive if needed) and then test
    // to see if it is a meson with fmod(currId/100,5)==0.0
    // to see if it is a baryon with currId/1000==5
    int HadronFromBQuark(0), HadronFromBbarQuark(0);
    if (bQuarkDaughter1<bQuarkDaughter2) {
      cout << "Daughters of b Quark" << endl;
      for (int j=bQuarkDaughter1; j<=bQuarkDaughter2; ++j) {
    if (isBHadron(pythia.event[j].id())) {
      cout << "Fragmentation: b -> " << pythia.event[j].name() << endl;
      cout << "                 id " << pythia.event[j].id() << " at index " << j << endl; 
      HadronFromBQuark = j;
    }
      }
    }
    if (bbarQuarkDaughter1<bbarQuarkDaughter2) {
      cout << "Daughters of bbar Quark" << endl;
      for (int k=bbarQuarkDaughter1; k<=bbarQuarkDaughter2; ++k) {
    if (isBHadron(pythia.event[k].id())) {
      cout << "Fragmentation : bbar -> " << pythia.event[k].name()  << endl;
      cout << "                     id " << pythia.event[k].id() << " at index " << k << endl; 
      HadronFromBbarQuark = k;
    }
      }
    }
    // Search the daughters of the hadrons until electrons and positrons are found
    // if there are any from a semileptonic decay of a beauty hadron
    // Start with the b quark
    int Daughter(HadronFromBQuark), electronIndex(0), positronIndex(0);
    while (Daughter!=0) {
      cout << "Checking " << pythia.event[Daughter].name() << " for e+/e- daughters" << endl;
      if (pythia.event[Daughter].id()==-511) {
    // This is a Bbar0, enter its pT
    cout << "Filling Bbar0 pT" << endl;
    B0barmesonPt->Fill(pythia.event[Daughter].pT());
      }
      if (pythia.event[Daughter].id()==511) {
    // This is a B0, enter its pT
    cout << "Filling Bbar0 pT" << endl;
    B0mesonPt->Fill(pythia.event[Daughter].pT());
      }
      int nextDaughter1 = pythia.event[Daughter].daughter1();
      int nextDaughter2 = pythia.event[Daughter].daughter2();
      // search for electron or positron
      for (int iDaughter = nextDaughter1; iDaughter<=nextDaughter2; ++iDaughter) {
    if (pythia.event[iDaughter].id()==11) {
      cout << "Found electron" << endl;
      cout << pythia.event[iDaughter].name() << endl;
      electronIndex=iDaughter;
      electronFrombPt->Fill(pythia.event[electronIndex].pT());
      break;
    }
    if (pythia.event[iDaughter].id()==-11) {
      cout << "Found positron" << endl;
      cout << pythia.event[iDaughter].name() << endl;
      positronIndex=iDaughter;
      positronFrombPt->Fill(pythia.event[positronIndex].pT());
      break;
    }
      }// loop over daughters to check for e+e-
      
      // If we get here, that means there were no electrons nor positrons.
      // Set the Daughter index to zero now.
      Daughter = 0;
      // If any of the daughters is still a beauty-hadron, we can try again
      // and reset the Daughter index, but only if one of the daughters contains a
      // b quark.
      for (int jDaughter = nextDaughter1; jDaughter<=nextDaughter2; ++jDaughter) {
    if (isBHadron(pythia.event[jDaughter].id())) {
      //One of the daughters is a beauty hadron.
      Daughter = jDaughter;
    }
      }// loop over daughters to check for another b hadron
    }// end of search for electrons in all the daughters of the b quark
    
    // Now search among the daughters of the bbar quark
    Daughter=HadronFromBbarQuark;
    while (Daughter!=0) {
      cout << "Checking " << pythia.event[Daughter].name() << " for e+/e- daughters" << endl;
      if (pythia.event[Daughter].id()==-511) {
    // This is a Bbar0, enter its pT
    cout << "Filling Bbar0 pT" << endl;
    B0barmesonPt->Fill(pythia.event[Daughter].pT());
      }
      if (pythia.event[Daughter].id()==511) {
    // This is a B0, enter its pT
    cout << "Filling B0 pT" << endl;
    B0mesonPt->Fill(pythia.event[Daughter].pT());
      }
      int nextDaughter1 = pythia.event[Daughter].daughter1();
      int nextDaughter2 = pythia.event[Daughter].daughter2();
      // search for electron or positron
      for (int iDaughter = nextDaughter1; iDaughter<=nextDaughter2; ++iDaughter) {
    //cout << "daughter is a " << pythia.event[iDaughter].name() << endl;
    if (pythia.event[iDaughter].id()==11) {
      cout << "Found electron" << endl;
      cout << pythia.event[iDaughter].name() << endl;
      electronIndex=iDaughter;
      electronFrombPt->Fill(pythia.event[electronIndex].pT());
      break;
    }
    if (pythia.event[iDaughter].id()==-11) {
      cout << "Found positron" << endl;
      cout << pythia.event[iDaughter].name() << endl;
      positronIndex=iDaughter;
      positronFrombPt->Fill(pythia.event[positronIndex].pT());
      break;
    }
      }// loop over daughters to check for e+e-
      
      // If we get here, that means there were no electrons nor positrons.
      // Set the Daughter index to zero now.
      Daughter = 0;
      // If any of the daughters is still a beauty-hadron, we can try again
      // and reset the Daughter index, but only if one of the daughters contains a
      // b quark.
      for (int jDaughter = nextDaughter1; jDaughter<=nextDaughter2; ++jDaughter) {
    if (isBHadron(pythia.event[jDaughter].id())) {
      //One of the daughters is a beauty hadron.
      Daughter = jDaughter;
    }
      }// loop over daughters to check for another b hadron
    }//end of search for electron among daughters of bbar quark
    
    if (electronIndex!=0 && positronIndex!=0) {
      cout << "Found an e+e- pair from bbar" << endl;
      cout << "Ele 4-mom = " << pythia.event[electronIndex].p() << endl;
      cout << "Pos 4-mom = " << pythia.event[positronIndex].p() << endl;
      Vec4 epluseminus(pythia.event[electronIndex].p()+pythia.event[positronIndex].p());
      epluseminusMinv->Fill(epluseminus.mCalc());
    }
    else {
      cout << "No e+e- pair in event" << endl;
    }
    
  // End of event loop. Statistics. Histogram. Done.
  }// event loop
  pythia.statistics();
  //cout << mult << endl;

  //Write Output ROOT hisotgram into ROOT file
  TFile* outFile = new TFile("pythiaOutputHistos1M.root","RECREATE");
  multHist->Write();
  bquarkPt->Write();
  bbarquarkPt->Write();
  B0mesonPt->Write();
  B0barmesonPt->Write();
  electronFrombPt->Write();
  positronFrombPt->Write();
  epluseminusMinv->Write();
  outFile->Close();

  return 0;
}
