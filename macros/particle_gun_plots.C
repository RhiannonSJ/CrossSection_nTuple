#include "TMath.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include <vector>
#include <iostream>
#include <string>
#include "TFile.h"

void particle_gun_plots() {

  // Read in the sample files
  TFile *f_proton = new TFile( "/sbnd/app/users/rsjones/LArSoft_v06_63_00/LArSoft-v06_63_00/srcs/recoparameters/recoparameters/output_files/proton_pid_information.root" );
  TFile *f_muon   = new TFile( "/sbnd/app/users/rsjones/LArSoft_v06_63_00/LArSoft-v06_63_00/srcs/recoparameters/recoparameters/output_files/muon_pid_information.root" );

  TTree *t_proton_pid = (TTree*)f_proton->Get("fNt_pid");
  TTree *t_muon_pid   = (TTree*)f_muon->Get("fNt_pid");

  std::vector< float > pr_chi2_pr;
  std::vector< float > mu_chi2_pr;
  std::vector< float > pr_chi2_mu;
  std::vector< float > mu_chi2_mu;
  std::vector< float > pr_chi2_pi;
  std::vector< float > mu_chi2_pi;
  std::vector< float > pr_chi2_ka;
  std::vector< float > mu_chi2_ka;
  std::vector< float > pr_chi2_pr_scaled;
  std::vector< float > mu_chi2_pr_scaled;
  std::vector< float > pr_chi2_mu_scaled;
  std::vector< float > mu_chi2_mu_scaled;
  std::vector< float > pr_chi2_pi_scaled;
  std::vector< float > mu_chi2_pi_scaled;
  std::vector< float > pr_chi2_ka_scaled;
  std::vector< float > mu_chi2_ka_scaled;
  std::vector< float > pr_pida;
  std::vector< float > mu_pida;

  int n_entries_pr = t_proton_pid->GetEntries();
  int n_entries_mu = t_muon_pid->GetEntries();

  for( unsigned int i = 0; i < n_entries_pr; ++i ){
  
    t_proton_pid->GetEntry(i);

    pr_chi2_pr_scaled.push_back( t_proton_pid->GetLeaf("Chi2Proton")->GetValue() / double( t_proton_pid->GetLeaf("Ndf")->GetValue() ) );
    pr_chi2_mu_scaled.push_back( t_proton_pid->GetLeaf("Chi2Muon")->GetValue()   / double( t_proton_pid->GetLeaf("Ndf")->GetValue() ) );
    pr_chi2_pi_scaled.push_back( t_proton_pid->GetLeaf("Chi2Pion")->GetValue()   / double( t_proton_pid->GetLeaf("Ndf")->GetValue() ) );
    pr_chi2_ka_scaled.push_back( t_proton_pid->GetLeaf("Chi2Kaon")->GetValue()   / double( t_proton_pid->GetLeaf("Ndf")->GetValue() ) );
    pr_chi2_pr.push_back(        t_proton_pid->GetLeaf("Chi2Proton")->GetValue() );
    pr_chi2_mu.push_back(        t_proton_pid->GetLeaf("Chi2Muon")->GetValue() );
    pr_chi2_pi.push_back(        t_proton_pid->GetLeaf("Chi2Pion")->GetValue() );
    pr_chi2_ka.push_back(        t_proton_pid->GetLeaf("Chi2Kaon")->GetValue() );
    pr_pida.push_back(           t_proton_pid->GetLeaf("PIDA")->GetValue() );

  }
   
  for( unsigned int i = 0; i < n_entries_mu; ++i ){
  
    t_muon_pid->GetEntry(i);

    mu_chi2_pr_scaled.push_back( t_muon_pid->GetLeaf("Chi2Proton")->GetValue() / double( t_muon_pid->GetLeaf("Ndf")->GetValue() ) );
    mu_chi2_mu_scaled.push_back( t_muon_pid->GetLeaf("Chi2Muon")->GetValue()   / double( t_muon_pid->GetLeaf("Ndf")->GetValue() ) );
    mu_chi2_pi_scaled.push_back( t_muon_pid->GetLeaf("Chi2Pion")->GetValue()   / double( t_muon_pid->GetLeaf("Ndf")->GetValue() ) );
    mu_chi2_ka_scaled.push_back( t_muon_pid->GetLeaf("Chi2Kaon")->GetValue()   / double( t_muon_pid->GetLeaf("Ndf")->GetValue() ) );
    mu_chi2_pr.push_back(        t_muon_pid->GetLeaf("Chi2Proton")->GetValue() );
    mu_chi2_mu.push_back(        t_muon_pid->GetLeaf("Chi2Muon")->GetValue() );
    mu_chi2_pi.push_back(        t_muon_pid->GetLeaf("Chi2Pion")->GetValue() );
    mu_chi2_ka.push_back(        t_muon_pid->GetLeaf("Chi2Kaon")->GetValue() );
    mu_pida.push_back(           t_muon_pid->GetLeaf("PIDA")->GetValue() );

  }
  
  TH1D *h_pr_chi2_pr_scaled = new TH1D( "h_pr_chi2_pr_scaled", "Proton particle gun #chi^{2}/NDF proton", 50, 0, 10 );
  TH1D *h_mu_chi2_pr_scaled = new TH1D( "h_mu_chi2_pr_scaled", "Muon particle gun #chi^{2}/NDF proton",   50, 0, 4 );

  TH1D *h_pr_chi2_mu_scaled = new TH1D( "h_pr_chi2_mu_scaled", "Proton particle gun #chi^{2}/NDF muon",   50, 0, 10 );
  TH1D *h_mu_chi2_mu_scaled = new TH1D( "h_mu_chi2_mu_scaled", "Muon particle gun #chi^{2}/NDF muon",     50, 0, 4 );
  
  TH1D *h_pr_chi2_pi_scaled = new TH1D( "h_pr_chi2_pi_scaled", "Proton particle gun #chi^{2}/NDF pion",   50, 0, 10 );
  TH1D *h_mu_chi2_pi_scaled = new TH1D( "h_mu_chi2_pi_scaled", "Muon particle gun #chi^{2}/NDF pion",     50, 0, 4 );
  
  TH1D *h_pr_chi2_ka_scaled = new TH1D( "h_pr_chi2_ka_scaled", "Proton particle gun #chi^{2}/NDF kaon",   50, 0, 10 );
  TH1D *h_mu_chi2_ka_scaled = new TH1D( "h_mu_chi2_ka_scaled", "Muon particle gun #chi^{2}/NDF kaon",     50, 0, 4 );
  
  TH1D *h_pr_chi2_pr = new TH1D( "h_pr_chi2_pr", "Proton particle gun #chi^{2} proton", 60, 0, 60 );
  TH1D *h_mu_chi2_pr = new TH1D( "h_mu_chi2_pr", "Muon particle gun #chi^{2} proton",   60, 0, 60 );

  TH1D *h_pr_chi2_mu = new TH1D( "h_pr_chi2_mu", "Proton particle gun #chi^{2} muon",   60, 0, 60 );
  TH1D *h_mu_chi2_mu = new TH1D( "h_mu_chi2_mu", "Muon particle gun #chi^{2} muon",     60, 0, 60 );
  
  TH1D *h_pr_chi2_pi = new TH1D( "h_pr_chi2_pi", "Proton particle gun #chi^{2} pion",   60, 0, 60 );
  TH1D *h_mu_chi2_pi = new TH1D( "h_mu_chi2_pi", "Muon particle gun #chi^{2} pion",     60, 0, 60 );
  
  TH1D *h_pr_chi2_ka = new TH1D( "h_pr_chi2_ka", "Proton particle gun #chi^{2} kaon",   60, 0, 60 );
  TH1D *h_mu_chi2_ka = new TH1D( "h_mu_chi2_ka", "Muon particle gun #chi^{2} kaon",     60, 0, 60 );
  
  TH1D *h_pr_pida    = new TH1D( "h_pr_pida",    "Proton particle gun PIDA",            80, 0, 20 );
  TH1D *h_mu_pida    = new TH1D( "h_mu_pida",    "Muon particle gun PIDA",              80, 0, 20 );
  
  // Counter for number of events with more than 0 vertices below a certain distance
  for( int i = 0; i < n_entries_pr; ++i ){

    h_pr_chi2_pr_scaled->Fill( pr_chi2_pr_scaled[i] );
    h_pr_chi2_mu_scaled->Fill( pr_chi2_mu_scaled[i] );
    h_pr_chi2_pi_scaled->Fill( pr_chi2_pi_scaled[i] );
    h_pr_chi2_ka_scaled->Fill( pr_chi2_ka_scaled[i] );
    h_pr_chi2_pr->Fill( pr_chi2_pr[i] );
    h_pr_chi2_mu->Fill( pr_chi2_mu[i] );
    h_pr_chi2_pi->Fill( pr_chi2_pi[i] );
    h_pr_chi2_ka->Fill( pr_chi2_ka[i] );
    h_pr_pida->Fill(    pr_pida[i] );
  
  }
  for( int i = 0; i < n_entries_mu; ++i ){

    h_mu_chi2_pr_scaled->Fill( mu_chi2_pr_scaled[i] );
    h_mu_chi2_mu_scaled->Fill( mu_chi2_mu_scaled[i] );
    h_mu_chi2_pi_scaled->Fill( mu_chi2_pi_scaled[i] );
    h_mu_chi2_ka_scaled->Fill( mu_chi2_ka_scaled[i] );
    h_mu_chi2_pr->Fill( mu_chi2_pr[i] );
    h_mu_chi2_mu->Fill( mu_chi2_mu[i] );
    h_mu_chi2_pi->Fill( mu_chi2_pi[i] );
    h_mu_chi2_ka->Fill( mu_chi2_ka[i] );
    h_mu_pida->Fill(    mu_pida[i] );

  }

  gStyle->SetPalette(57);
  gStyle->SetNumberContours(250);

  h_pr_chi2_pr_scaled->SetStats(kFALSE);
  h_pr_chi2_pr_scaled->SetLineColor(2);
  h_pr_chi2_pr_scaled->GetXaxis()->SetTitle( "Proton #chi^{2}/NDF" );
  h_pr_chi2_pr_scaled->GetXaxis()->SetTitleOffset(1.2);
  
  h_pr_chi2_mu_scaled->SetStats(kFALSE);
  h_pr_chi2_mu_scaled->SetLineColor(3);
  h_pr_chi2_mu_scaled->GetXaxis()->SetTitle( "Proton #chi^{2}/NDF" );
  h_pr_chi2_mu_scaled->GetXaxis()->SetTitleOffset(1.2);

  h_pr_chi2_pi_scaled->SetStats(kFALSE);
  h_pr_chi2_pi_scaled->SetLineColor(4);
  h_pr_chi2_pi_scaled->GetXaxis()->SetTitle( "Proton #chi^{2}/NDF" );
  h_pr_chi2_pi_scaled->GetXaxis()->SetTitleOffset(1.2);

  h_pr_chi2_ka_scaled->SetStats(kFALSE);
  h_pr_chi2_ka_scaled->SetLineColor(6);
  h_pr_chi2_ka_scaled->GetXaxis()->SetTitle( "Proton #chi^{2}/NDF" );
  h_pr_chi2_ka_scaled->GetXaxis()->SetTitleOffset(1.2);

  h_mu_chi2_pr_scaled->SetStats(kFALSE);
  h_mu_chi2_pr_scaled->SetLineColor(2);
  h_mu_chi2_pr_scaled->GetXaxis()->SetTitle( "Muon #chi^{2}/NDF" );
  h_mu_chi2_pr_scaled->GetXaxis()->SetTitleOffset(1.2);

  h_mu_chi2_mu_scaled->SetStats(kFALSE);
  h_mu_chi2_mu_scaled->SetLineColor(3);
  h_mu_chi2_mu_scaled->GetXaxis()->SetTitle( "Muon #chi^{2}/NDF" );
  h_mu_chi2_mu_scaled->GetXaxis()->SetTitleOffset(1.2);
  
  h_mu_chi2_pi_scaled->SetStats(kFALSE);
  h_mu_chi2_pi_scaled->SetLineColor(4);
  h_mu_chi2_pi_scaled->GetXaxis()->SetTitle( "Muon #chi^{2}/NDF" );
  h_mu_chi2_pi_scaled->GetXaxis()->SetTitleOffset(1.2);
  
  h_mu_chi2_ka_scaled->SetStats(kFALSE);
  h_mu_chi2_ka_scaled->SetLineColor(6);
  h_mu_chi2_ka_scaled->GetXaxis()->SetTitle( "Muon #chi^{2}/NDF" );
  h_mu_chi2_ka_scaled->GetXaxis()->SetTitleOffset(1.2);
  
  h_pr_chi2_pr->SetStats(kFALSE);
  h_pr_chi2_pr->SetLineColor(2);
  h_pr_chi2_pr->GetXaxis()->SetTitle( "Proton #chi^{2}" );
  h_pr_chi2_pr->GetXaxis()->SetTitleOffset(1.2);
  
  h_pr_chi2_mu->SetStats(kFALSE);
  h_pr_chi2_mu->SetLineColor(3);
  h_pr_chi2_mu->GetXaxis()->SetTitle( "Proton #chi^{2}" );
  h_pr_chi2_mu->GetXaxis()->SetTitleOffset(1.2);

  h_pr_chi2_pi->SetStats(kFALSE);
  h_pr_chi2_pi->SetLineColor(4);
  h_pr_chi2_pi->GetXaxis()->SetTitle( "Proton #chi^{2}" );
  h_pr_chi2_pi->GetXaxis()->SetTitleOffset(1.2);

  h_pr_chi2_ka->SetStats(kFALSE);
  h_pr_chi2_ka->SetLineColor(6);
  h_pr_chi2_ka->GetXaxis()->SetTitle( "Proton #chi^{2}" );
  h_pr_chi2_ka->GetXaxis()->SetTitleOffset(1.2);

  h_mu_chi2_pr->SetStats(kFALSE);
  h_mu_chi2_pr->SetLineColor(2);
  h_mu_chi2_pr->GetXaxis()->SetTitle( "Muon #chi^{2}" );
  h_mu_chi2_pr->GetXaxis()->SetTitleOffset(1.2);

  h_mu_chi2_mu->SetStats(kFALSE);
  h_mu_chi2_mu->SetLineColor(3);
  h_mu_chi2_mu->GetXaxis()->SetTitle( "Muon #chi^{2}" );
  h_mu_chi2_mu->GetXaxis()->SetTitleOffset(1.2);
  
  h_mu_chi2_pi->SetStats(kFALSE);
  h_mu_chi2_pi->SetLineColor(4);
  h_mu_chi2_pi->GetXaxis()->SetTitle( "Muon #chi^{2}" );
  h_mu_chi2_pi->GetXaxis()->SetTitleOffset(1.2);
  
  h_mu_chi2_ka->SetStats(kFALSE);
  h_mu_chi2_ka->SetLineColor(6);
  h_mu_chi2_ka->GetXaxis()->SetTitle( "Muon #chi^{2}" );
  h_mu_chi2_ka->GetXaxis()->SetTitleOffset(1.2);
  
  h_pr_pida->SetStats(kFALSE);
  h_pr_pida->SetLineColor(2);
  h_pr_pida->GetXaxis()->SetTitle( "PIDA" );
  h_pr_pida->GetXaxis()->SetTitleOffset(1.2);

  h_mu_pida->SetStats(kFALSE);
  h_mu_pida->SetLineColor(4);
  h_mu_pida->GetXaxis()->SetTitle( "PIDA" );
  h_mu_pida->GetXaxis()->SetTitleOffset(1.2);

  TLegend *l_chi2_pr  = new TLegend(0.7, 0.6, 0.85, 0.85);
  TLegend *l_chi2_mu  = new TLegend(0.7, 0.6, 0.85, 0.85);
  TLegend *l_pida     = new TLegend(0.7, 0.6, 0.85, 0.85);

  l_chi2_pr->AddEntry( h_pr_chi2_pr, "Proton", "l");
  l_chi2_pr->AddEntry( h_pr_chi2_mu, "Muon",   "l");
  l_chi2_pr->AddEntry( h_pr_chi2_pi, "Pion",   "l");
  l_chi2_pr->AddEntry( h_pr_chi2_ka, "Kaon",   "l");

  l_chi2_mu->AddEntry( h_mu_chi2_pr, "Proton", "l");
  l_chi2_mu->AddEntry( h_mu_chi2_mu, "Muon",   "l");
  l_chi2_mu->AddEntry( h_mu_chi2_pi, "Pion",   "l");
  l_chi2_mu->AddEntry( h_mu_chi2_ka, "Kaon",   "l");

  l_pida->AddEntry( h_pr_pida,       "Proton", "l");
  l_pida->AddEntry( h_mu_pida,       "Muon", "l");

  TCanvas *c  = new TCanvas( "c", "Canvas", 800, 600 );

  //c->SetLogy();
  //
  // Draw!!
  //
  h_pr_chi2_pr_scaled->Draw();
  h_pr_chi2_mu_scaled->Draw("same");
  h_pr_chi2_pi_scaled->Draw("same");
  h_pr_chi2_ka_scaled->Draw("same");
  l_chi2_pr->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_63_00/LArSoft-v06_63_00/srcs/recoparameters/recoparameters/pg_plots/scaled_proton_chi2.root");
  c->Clear();

  h_mu_chi2_mu_scaled->Draw();
  h_mu_chi2_pr_scaled->Draw("same");
  h_mu_chi2_pi_scaled->Draw("same");
  h_mu_chi2_ka_scaled->Draw("same");
  l_chi2_mu->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_63_00/LArSoft-v06_63_00/srcs/recoparameters/recoparameters/pg_plots/scaled_muon_chi2.root");
  c->Clear();

  h_pr_chi2_pr->Draw();
  h_pr_chi2_mu->Draw("same");
  h_pr_chi2_pi->Draw("same");
  h_pr_chi2_ka->Draw("same");
  l_chi2_pr->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_63_00/LArSoft-v06_63_00/srcs/recoparameters/recoparameters/pg_plots/proton_chi2.root");
  c->Clear();

  h_mu_chi2_mu->Draw();
  h_mu_chi2_pr->Draw("same");
  h_mu_chi2_pi->Draw("same");
  h_mu_chi2_ka->Draw("same");
  l_chi2_mu->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_63_00/LArSoft-v06_63_00/srcs/recoparameters/recoparameters/pg_plots/muon_chi2.root");
  c->Clear();

  h_pr_pida->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_63_00/LArSoft-v06_63_00/srcs/recoparameters/recoparameters/pg_plots/proton_pida.root");
  c->Clear();
  
  h_mu_pida->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_63_00/LArSoft-v06_63_00/srcs/recoparameters/recoparameters/pg_plots/muon_pida.root");
  c->Clear();
  
  h_mu_pida->SetTitle("Particle gun PIDA distributions");
  h_pr_pida->SetTitle("Particle gun PIDA distributions");
  h_mu_pida->Draw();
  h_pr_pida->Draw("same");
  l_pida->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_63_00/LArSoft-v06_63_00/srcs/recoparameters/recoparameters/pg_plots/pida.root");
  c->Clear();

}
