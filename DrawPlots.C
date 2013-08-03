//This is a macro to get histograms and print them to a pdf 
//or other format like png, jpeg, etc

void fatal(TString msg) { printf("Bad: %s\n\n",msg.Data()); abort(); }
TH1 *GetHisto(TFile *f, TString hn) {
  TH1* h = (TH1*)f->Get(hn);
  if (h==NULL) fatal("Cannot get histo "+hn+" in file "+f->GetName());
  return h;
}

int runs = 1;

std::string run="0";


void DrawPlots(TString fn="rootplots/ttbar_partonLevel_histos_"+run+".root") {
  TFile *f = TFile::Open(fn);
  if (f==NULL) fatal("Cannot open file "+fn);
  
  TString pdf("Pdfs/top_plot_"+run+".pdf");
  c = new TCanvas();

  c->Print(pdf+"[");

  //divide canvas into quadrants
    c->Divide(2,2);

    TLatex tex; tex.SetNDC();
    tex.DrawLatex(0.35,0.49,"Leptonic pseudotop");
    c->cd(1); GetHisto(f,"pseudo_mass_matched")->Draw("e1");
    c->Update();

    c->cd(2); GetHisto(f,"pseudo_mass_unmatched")->Draw("e1");
    c->Update();

    c->cd(3); GetHisto(f,"pseudo_pt_matched")->Draw("e1 ");
    c->Update();

    c->cd(4); GetHisto(f,"pseudo_pt_unmatched") ->Draw("e1 ");
    c->Print(pdf);



    c->cd(1); GetHisto(f,"pseudo_eta_matched")->Draw("e1 ");
    c->Update();

    c->cd(2); GetHisto(f,"pseudo_eta_unmatched")->Draw("e1 ");

    c->Update();

    c->cd(3); GetHisto(f,"whadronic")->Draw("e1 y+ ");
    c->Update();

    c->cd(4); GetHisto(f,"wleptonic")->Draw("e1");
    c->Print(pdf);


    c->cd(1); GetHisto(f,"m_top")->Draw("e1");

    c->cd(2); GetHisto(f,"eta_top")->Draw("e1 c ");

    c->cd(3); GetHisto(f,"pt_top") ->Draw("e1 ");

    c->cd(4); GetHisto(f,"y_top")->Draw("e1 ");
    c->Print(pdf);



    c->cd(1); GetHisto(f,"eff_vs_topPt")->Draw("e1 ");

    c->cd(2); GetHisto(f,"neff_vs_topPt")->Draw("e1  ");

    c->cd(3); GetHisto(f,"eff_vs_topMass") ->Draw("e1 ");

    c->cd(4); GetHisto(f,"eff_vs_topeta") ->Draw("e1");
    c->Print(pdf);


    c->cd(1); GetHisto(f,"pseudo_btop_mass")->Draw("e1 ");

    c->cd(2); GetHisto(f,"pseudo_HW_top")->Draw("e1 ");

    c->cd(3); GetHisto(f,"pseudb")->Draw("e1");

    c->cd(4); GetHisto(f,"eff_vs_bmass")->Draw("e1");
    c->Print(pdf);

    c->Print(pdf+"]");


  }
