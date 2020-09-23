/***********************************************
* Software developement for WASA-at-COSY
* (c) 2005-2020 The WASA-at-COSY Collaboration
* Aleksander K.                 2019-11
* This software is distributed under the terms
  of the GNU General Public Licence v3.0
*
* Modified 2020-05
***********************************************/

//Macro to calculate the efficiency for the pd -> (3He-eta)_bound -> dp pi0 reaction

#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TPaveLabel.h>
#include <TFrame.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TPaveText.h>
#include <TInterpreter.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPaletteAxis.h>
#include <TLegend.h>
#include <TLine.h>
#include <cassert>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TArrow.h>
#include <TObjArray.h>
#include <vector>
//#include <TFractionFitter.h>
#include <TMinuit.h>
#include <Riostream.h>

void efficiency() {

    TFile *myFile[20];

    myFile[0] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0.root","READ");    
    myFile[1] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0_G20_Bs10.root","READ");
    myFile[2] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0_G20_Bs20.root","READ");
    myFile[3] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0_G20_Bs30.root","READ");
    myFile[4] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0_G20_Bs40.root","READ");
    myFile[5] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0_G30_Bs20.root","READ");
    myFile[6] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0_G30_Bs30.root","READ");
    myFile[7] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0_G30_Bs40.root","READ");
    myFile[8] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0_G40_Bs20.root","READ");
    myFile[9] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0_G40_Bs30.root","READ");
    myFile[10] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0_G40_Bs40.root","READ");
    myFile[11] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0_G50_Bs10.root","READ");
    myFile[12] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0_G50_Bs20.root","READ");
    myFile[13] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0_G50_Bs30.root","READ");
    myFile[14] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0_G50_Bs40.root","READ");
    myFile[15] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0_G150_Bs20.root","READ");    
    myFile[16] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0_ProtDist.root","READ");
    myFile[17] = new TFile("input/MC-newcuts-AddGammaCut-pd-pdpi0.root","READ");

    TH1D *hGenerated[20];
    TH1D *hAccepted[20];
    TH1D *hGeomAccept[20];

    TH1D *hAcceptance[3];

    hAcceptance[0] = new TH1D("Acceptance_N","",40,-70,30);
    hAcceptance[1] = new TH1D("Acceptance_P","",40,-70,30);

    TH1D *hEfficiency[20];

    TH1D *hGeneratedSystErr[1];
    TH1D *hAcceptedSystErr[11];
    TH1D *hEfficiencySystErr[11];

    for (Int_t i = 0; i < 16; i++) {

        myFile[i]->cd("Histograms");

        hGenerated[i] = (TH1D*)gDirectory->Get("WMC/hGenerated_Q");
        hAccepted[i] = (TH1D*)gDirectory->Get("DATA_lev2_cut4/hQ_lev2_cut4");

    }

    myFile[0]->cd("Histograms");
    hGeomAccept[0] = (TH1D*)gDirectory->Get("WMC/hAccepted_Q");

    myFile[16]->cd("Histograms");
    hGenerated[16] = (TH1D*)gDirectory->Get("WMC/hGenerated_Q");
    hGeomAccept[16] = (TH1D*)gDirectory->Get("WMC/hAccepted_Q");

    myFile[17]->cd("Histograms");
    hGenerated[17] = (TH1D*)gDirectory->Get("WMC/hGenerated_Q");
    hAccepted[17] = (TH1D*)gDirectory->Get("DATA_lev2_cut4/hQ_lev2_cut4");

    Double_t Ngen = 0.;
    Double_t Nacc = 0.;
    Double_t statErrEff = 0.;

    for (Int_t j = 1; j < 41; j++) {

        Ngen = hGenerated[0]->GetBinContent(j);
        Nacc = hGeomAccept[0]->GetBinContent(j);
        statErrEff = TMath::Sqrt((Nacc/(Ngen*Ngen)) + (Nacc*Nacc)/(Ngen*Ngen*Ngen));

        hAcceptance[0]->SetBinContent(j,(Nacc/Ngen)*100);
        hAcceptance[0]->SetBinError(j,statErrEff*100);

    }

    for (Int_t j = 1; j < 41; j++) {

        Ngen = hGenerated[16]->GetBinContent(j);
        Nacc = hGeomAccept[16]->GetBinContent(j);
        statErrEff = TMath::Sqrt((Nacc/(Ngen*Ngen)) + (Nacc*Nacc)/(Ngen*Ngen*Ngen));

        hAcceptance[1]->SetBinContent(j,(Nacc/Ngen)*100);
        hAcceptance[1]->SetBinError(j,statErrEff*100);

    }

    for (Int_t k = 0; k < 16; k++) {

        hEfficiency[k] = new TH1D("Efficiency","",40,-70,30);

        for (Int_t l = 1; l < 41; l++) {

            Ngen = hGenerated[k]->GetBinContent(l);
            Nacc = hAccepted[k]->GetBinContent(l);

            statErrEff = TMath::Sqrt((Nacc/(Ngen*Ngen)) + (Nacc*Nacc)/(Ngen*Ngen*Ngen));

            hEfficiency[k]->SetBinContent(l,(Nacc/Ngen)*100);
            hEfficiency[k]->SetBinError(l,statErrEff*100);

        }

    }

    hEfficiency[17] = new TH1D("Efficiency","",40,-70,30);

    for (Int_t l = 1; l < 41; l++) {

        Ngen = hGenerated[17]->GetBinContent(l);
        Nacc = hAccepted[17]->GetBinContent(l);

        statErrEff = TMath::Sqrt((Nacc/(Ngen*Ngen)) + (Nacc*Nacc)/(Ngen*Ngen*Ngen));

        hEfficiency[17]->SetBinContent(l,(Nacc/Ngen)*100);
        hEfficiency[17]->SetBinError(l,statErrEff*100);

    }

    //Systematics
    Double_t NgenSystErr = 0.;
    Double_t NaccSystErr = 0.;
    Double_t statErrEffSystErr = 0.;

    myFile[0]->cd("Histograms");

    hGeneratedSystErr[0] = (TH1D*)gDirectory->Get("WMC/hGenerated_Q");
    hAcceptedSystErr[0] = (TH1D*)gDirectory->Get("DATA_lev2_cut4/hQ_lev2_cut4");

    hAcceptedSystErr[1] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_1_cut4");
    hAcceptedSystErr[2] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_2_cut4");
    hAcceptedSystErr[3] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_cut4_1_1");
    hAcceptedSystErr[4] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_cut4_1_2");
    hAcceptedSystErr[5] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_cut4_2_1");
    hAcceptedSystErr[6] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_cut4_2_2");
    hAcceptedSystErr[7] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_cut4_3_1");
    hAcceptedSystErr[8] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_cut4_3_2");
    hAcceptedSystErr[9] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_cut4_4_1");
    hAcceptedSystErr[10] = (TH1D*)gDirectory->Get("Systematcs/hQ_lev2_cut4_4_2");

    for (Int_t p = 0; p < 11; p++) {

        hEfficiencySystErr[p] = new TH1D(Form("hEfficiencySystErr_%d",p),"",40,-70,30);

        for (Int_t l = 1; l < 41; l++) {

            NgenSystErr = hGeneratedSystErr[0]->GetBinContent(l);
            NaccSystErr = hAcceptedSystErr[p]->GetBinContent(l);

            statErrEffSystErr = TMath::Sqrt((NaccSystErr/(NgenSystErr*NgenSystErr)) + (NaccSystErr*NaccSystErr)/(NgenSystErr*NgenSystErr*NgenSystErr));

            hEfficiencySystErr[p]->SetBinContent(l,(NaccSystErr/NgenSystErr));
            hEfficiencySystErr[p]->SetBinError(l,statErrEffSystErr);

        }

    }

    ////
    TFile* newFile = new TFile("output/files/Efficiency.root","RECREATE");

    newFile->cd();

    hGenerated[0]->Write("hGenerated");
    hGenerated[1]->Write("hGenerated_G20_Bs10");
    hGenerated[2]->Write("hGenerated_G20_Bs20");
    hGenerated[3]->Write("hGenerated_G20_Bs30");
    hGenerated[4]->Write("hGenerated_G20_Bs40");
    hGenerated[5]->Write("hGenerated_G30_Bs20");
    hGenerated[6]->Write("hGenerated_G30_Bs30");
    hGenerated[7]->Write("hGenerated_G30_Bs40");
    hGenerated[8]->Write("hGenerated_G40_Bs20");
    hGenerated[9]->Write("hGenerated_G40_Bs30");
    hGenerated[10]->Write("hGenerated_G40_Bs40");
    hGenerated[11]->Write("hGenerated_G50_Bs10");
    hGenerated[12]->Write("hGenerated_G50_Bs20");
    hGenerated[13]->Write("hGenerated_G50_Bs30");
    hGenerated[14]->Write("hGenerated_G50_Bs40");
    hGenerated[15]->Write("hGenerated_G150_Bs20");
    hGenerated[16]->Write("hGenerated_proton");
    hGenerated[17]->Write("hGenerated_bkgnd");

    hAccepted[0]->Write("hAccepted");
    hAccepted[1]->Write("hAccepted_G20_Bs10");
    hAccepted[2]->Write("hAccepted_G20_Bs20");
    hAccepted[3]->Write("hAccepted_G20_Bs30");
    hAccepted[4]->Write("hAccepted_G20_Bs40");
    hAccepted[5]->Write("hAccepted_G30_Bs20");
    hAccepted[6]->Write("hAccepted_G30_Bs30");
    hAccepted[7]->Write("hAccepted_G30_Bs40");
    hAccepted[8]->Write("hAccepted_G40_Bs20");
    hAccepted[9]->Write("hAccepted_G40_Bs30");
    hAccepted[10]->Write("hAccepted_G40_Bs40");
    hAccepted[11]->Write("hAccepted_G50_Bs10");
    hAccepted[12]->Write("hAccepted_G50_Bs20");
    hAccepted[13]->Write("hAccepted_G50_Bs30");
    hAccepted[14]->Write("hAccepted_G50_Bs40");
    hAccepted[15]->Write("hAccepted_G150_Bs20");
    hAccepted[17]->Write("hAccepted_bkgnd");

    hAcceptance[0]->Write("Acceptance_N*");
    hAcceptance[1]->Write("Acceptance_proton");

    hEfficiency[0]->Write("hEfficiency");
    hEfficiency[1]->Write("hEfficiency_G20_Bs10");
    hEfficiency[2]->Write("hEfficiency_G20_Bs20");
    hEfficiency[3]->Write("hEfficiency_G20_Bs30");
    hEfficiency[4]->Write("hEfficiency_G20_Bs40");
    hEfficiency[5]->Write("hEfficiency_G30_Bs20");
    hEfficiency[6]->Write("hEfficiency_G30_Bs30");
    hEfficiency[7]->Write("hEfficiency_G30_Bs40");
    hEfficiency[8]->Write("hEfficiency_G40_Bs20");
    hEfficiency[9]->Write("hEfficiency_G40_Bs30");
    hEfficiency[10]->Write("hEfficiency_G40_Bs40");
    hEfficiency[11]->Write("hEfficiency_G50_Bs10");
    hEfficiency[12]->Write("hEfficiency_G50_Bs20");
    hEfficiency[13]->Write("hEfficiency_G50_Bs30");
    hEfficiency[14]->Write("hEfficiency_G50_Bs40");
    hEfficiency[15]->Write("hEfficiency_G150_Bs20");
    hEfficiency[17]->Write("hEfficiency_bkgnd");

    hEfficiencySystErr[0]->Write("hEfficiency_lev2_cut4");
    hEfficiencySystErr[1]->Write("hEfficiency_lev2_1_cut4");
    hEfficiencySystErr[2]->Write("hEfficiency_lev2_2_cut4");
    hEfficiencySystErr[3]->Write("hEfficiency_lev2_cut4_1_1");
    hEfficiencySystErr[4]->Write("hEfficiency_lev2_cut4_1_2");
    hEfficiencySystErr[5]->Write("hEfficiency_lev2_cut4_2_1");
    hEfficiencySystErr[6]->Write("hEfficiency_lev2_cut4_2_2");
    hEfficiencySystErr[7]->Write("hEfficiency_lev2_cut4_3_1");
    hEfficiencySystErr[8]->Write("hEfficiency_lev2_cut4_3_2");
    hEfficiencySystErr[9]->Write("hEfficiency_lev2_cut4_4_1");
    hEfficiencySystErr[10]->Write("hEfficiency_lev2_cut4_4_2");

    newFile->Close();

    ////
    gStyle->SetOptStat(kFALSE);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadRightMargin(0.10);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadBottomMargin(0.15);

    TCanvas* myCanvas01 = new TCanvas;

    //hAcceptance[0]->SetTitle("Acceptance");    
    hAcceptance[0]->GetXaxis()->SetTitle("excess energy [MeV]");
    hAcceptance[0]->GetXaxis()->SetTitleSize(0.06);
    hAcceptance[0]->GetXaxis()->SetTitleOffset(0.95);
    hAcceptance[0]->GetXaxis()->SetLabelSize(0.05);
    //hAcceptance[0]->GetXaxis()->SetRangeUser(-70.,30.);
    hAcceptance[0]->GetYaxis()->SetTitle("acceptance [%]");
    hAcceptance[0]->GetYaxis()->SetTitleSize(0.06);
    hAcceptance[0]->GetYaxis()->SetTitleOffset(0.95);
    hAcceptance[0]->GetYaxis()->SetLabelSize(0.05);
    hAcceptance[0]->GetYaxis()->SetRangeUser(0.,100.);

    hAcceptance[0]->SetLineWidth(1);
    hAcceptance[0]->SetLineColor(kCyan-3);
    hAcceptance[0]->SetMarkerColor(kCyan-3);
    hAcceptance[0]->SetMarkerSize(0.85);
    hAcceptance[0]->SetMarkerStyle(23);
    //hAcceptance[0]->Scale(100);
    hAcceptance[0]->Draw("pE1");

    hAcceptance[1]->SetLineWidth(1);
    hAcceptance[1]->SetLineColor(kOrange+1);
    hAcceptance[1]->SetMarkerColor(kOrange+1);
    hAcceptance[1]->SetMarkerSize(0.75);
    hAcceptance[1]->SetMarkerStyle(20);
    hAcceptance[1]->Draw("same pE1");

    hEfficiency[0]->SetLineWidth(1);
    hEfficiency[0]->SetLineColor(kOrange+1);
    hEfficiency[0]->SetMarkerColor(kOrange+1);
    hEfficiency[0]->SetMarkerSize(0.75);
    hEfficiency[0]->SetMarkerStyle(20);
    //hEfficiency[0]->Scale(100);
    //hEfficiency[0]->Draw("same  pE1");

    hEfficiency[17]->SetLineWidth(1);
    hEfficiency[17]->SetLineColor(kMagenta+2);
    hEfficiency[17]->SetMarkerColor(kMagenta+2);
    hEfficiency[17]->SetMarkerSize(0.75);
    hEfficiency[17]->SetMarkerStyle(21);
    //hEfficiency[17]->Scale(100);
    //hEfficiency[17]->Draw("same  pE1");

    TLegend *myLegend01 = new TLegend(0.500, 0.750, 0.885, 0.885);
    myLegend01->SetFillStyle(1); myLegend01->SetFillColor(0); myLegend01->SetLineColor(0); myLegend01->SetTextSize(0.04);
    //myLegend01->AddEntry( hAcceptance[1], "geometrical acceptance (p inside ^{3}He)", "lp");
    //myLegend01->AddEntry( hAcceptance[0], "geometrical acceptance (N* in N*-d)", "lp");
    myLegend01->AddEntry( hAcceptance[0], "geometrical acceptance", "lp");
    myLegend01->AddEntry( hEfficiency[0], "efficiency", "lp");
    //myLegend01->AddEntry( hEfficiency[17], "bkgnd efficiency: pd #rightarrow dp#pi^{0}", "lp");
    myLegend01->Draw("same");

    myCanvas01->Print("output/plots/hAcceptance.png","png");
    myCanvas01->Print("output/plots/hAcceptance.eps","eps");

    //
    TCanvas* myCanvas01a = new TCanvas;

    hAcceptance[0]->GetXaxis()->SetTitle("\\hbox{energia dostępna, MeV}");
    hAcceptance[0]->GetYaxis()->SetTitle("\\hbox{akceptancja, %}");
    //hAcceptance[0]->GetYaxis()->SetTitle("\\hbox{wydajność, %}");

    hAcceptance[0]->Draw("pE1");
    hAcceptance[1]->Draw("same pE1");
    //hEfficiency[0]->Draw("same pE");
    //hEfficiency[17]->Draw("same pE");

    TLegend *myLegend01a = new TLegend(0.450, 0.720, 0.885, 0.885);
    //TLegend *myLegend01a = new TLegend(0.457, 0.607, 0.885, 0.885);
    myLegend01a->SetFillStyle(1001); myLegend01a->SetFillColor(19); myLegend01a->SetLineColor(1); myLegend01a->SetTextSize(0.04); myLegend01a->SetBorderSize(5);
    myLegend01a->AddEntry((TObject*)0, "Akceptancja geometryczna:", "");
    myLegend01a->AddEntry( hAcceptance[1], "proton w ^{3}He", "elp");
    myLegend01a->AddEntry( hAcceptance[0], "\\hbox{N* w układzie N*-d}", "elp");
    //myLegend01a->AddEntry( hAcceptance[0], "akceptancja geometryczna", "elp");
    //myLegend01a->AddEntry((TObject*)0, "\\hbox{Wydajność rekonstrukcji:}", "");
    //myLegend01a->AddEntry( hEfficiency[0], "pd #rightarrow (^{3}He#eta)_{            } #rightarrow dp#pi^{0}", "elp");
    //myLegend01a->AddEntry( hEfficiency[17], "pd #rightarrow dp#pi^{0}", "elp");
    myLegend01a->Draw();

    TPaveText *bnd = new TPaveText(9.5, 72.5, 9.5, 72.5,"boubd");
    bnd->SetTextSize(0.025);
    bnd->SetFillColor(0);
    bnd->SetTextColor(1);
    bnd->SetTextAlign(22);
    bnd->AddText("\\hbox{związany}");
    //bnd->Draw("same");

    myCanvas01a->Print("output/plots/hAcceptance_pl.png","png");
    myCanvas01a->Print("output/plots/hAcceptance_pl.eps","eps");

    //
    TCanvas* myCanvas02 = new TCanvas;

    //hEfficiency[1]->SetTitle("Efficiency");
    hEfficiency[1]->GetXaxis()->SetTitle("excess energy, MeV");
    hEfficiency[1]->GetXaxis()->SetTitleSize(0.06);
    hEfficiency[1]->GetXaxis()->SetTitleOffset(0.95);
    hEfficiency[1]->GetXaxis()->SetLabelSize(0.05);
    //hEfficiency[1]->GetXaxis()->SetRangeUser(-70.,30.);
    hEfficiency[1]->GetYaxis()->SetTitle("efficiency, %");
    hEfficiency[1]->GetYaxis()->SetTitleSize(0.06);
    hEfficiency[1]->GetYaxis()->SetTitleOffset(0.85);
    hEfficiency[1]->GetYaxis()->SetLabelSize(0.05);
    hEfficiency[1]->GetYaxis()->SetRangeUser(0.,25.);

    hEfficiency[1]->SetLineWidth(1);
    hEfficiency[1]->SetLineColor(51);
    hEfficiency[1]->SetMarkerColor(51);
    hEfficiency[1]->SetMarkerSize(0.62);
    hEfficiency[1]->SetMarkerStyle(20);
    hEfficiency[1]->Draw("p E1");

    hEfficiency[5]->SetLineWidth(1);
    hEfficiency[5]->SetLineColor(61);
    hEfficiency[5]->SetMarkerColor(61);
    hEfficiency[5]->SetMarkerSize(0.62);
    hEfficiency[5]->SetMarkerStyle(21);
    hEfficiency[5]->Draw("same p E1");

    hEfficiency[9]->SetLineWidth(1);
    hEfficiency[9]->SetLineColor(68);
    hEfficiency[9]->SetMarkerColor(68);
    hEfficiency[9]->SetMarkerSize(0.62);
    hEfficiency[9]->SetMarkerStyle(22);
    hEfficiency[9]->Draw("same p E1");

    hEfficiency[14]->SetLineWidth(1);
    hEfficiency[14]->SetLineColor(94);
    hEfficiency[14]->SetMarkerColor(94);
    hEfficiency[14]->SetMarkerSize(0.62);
    hEfficiency[14]->SetMarkerStyle(24);
    hEfficiency[14]->Draw("same p E1");

    hEfficiency[15]->SetLineWidth(1);
    hEfficiency[15]->SetLineColor(98);
    hEfficiency[15]->SetMarkerColor(98);
    hEfficiency[15]->SetMarkerSize(0.62);
    hEfficiency[15]->SetMarkerStyle(25);
    hEfficiency[15]->Draw("same p E1");

    TLegend *myLegend02 = new TLegend(0.460, 0.570, 0.885, 0.885);
    myLegend02->SetFillStyle(1); myLegend02->SetFillColor(0); myLegend02->SetLineColor(0); myLegend02->SetTextSize(0.04);
    myLegend02->AddEntry(hEfficiency[1], "B_{s} = 10 MeV, #Gamma = 20 MeV", "elp");
    myLegend02->AddEntry(hEfficiency[5], "B_{s} = 20 MeV, #Gamma = 30 MeV", "elp");
    myLegend02->AddEntry(hEfficiency[9], "B_{s} = 30 MeV, #Gamma = 40 MeV", "elp");
    myLegend02->AddEntry(hEfficiency[14], "B_{s} = 40 MeV, #Gamma = 50 MeV", "elp");
    myLegend02->AddEntry(hEfficiency[15], "B_{s} = 20 MeV, #Gamma = 150 MeV", "elp");
    myLegend02->Draw("same");

    myCanvas02->Print("output/plots/hEfficiency.png","png");
    myCanvas02->Print("output/plots/hEfficiency.eps","eps");

    //
    TCanvas* myCanvas02a = new TCanvas;

    hEfficiency[1]->GetXaxis()->SetTitle("\\hbox{energia dostępna, MeV}");
    hEfficiency[1]->GetYaxis()->SetTitle("\\hbox{wydajność, %}");
    hEfficiency[1]->Draw("p E1");

    hEfficiency[5]->Draw("same p E1");
    hEfficiency[9]->Draw("same p E1");
    hEfficiency[14]->Draw("same p E1");
    hEfficiency[15]->Draw("same p E1");

    TLegend *myLegend02a = new TLegend(0.460, 0.590, 0.885, 0.885);
    myLegend02a->SetFillStyle(1001); myLegend02a->SetFillColor(19); myLegend02a->SetLineColor(1); myLegend02a->SetTextSize(0.04); myLegend02a->SetBorderSize(5);
    myLegend02a->AddEntry(hEfficiency[1], "B_{s} = 10 MeV, #Gamma = 20 MeV", "elp");
    myLegend02a->AddEntry(hEfficiency[5], "B_{s} = 20 MeV, #Gamma = 30 MeV", "elp");
    myLegend02a->AddEntry(hEfficiency[9], "B_{s} = 30 MeV, #Gamma = 40 MeV", "elp");
    myLegend02a->AddEntry(hEfficiency[14], "B_{s} = 40 MeV, #Gamma = 50 MeV", "elp");
    myLegend02a->AddEntry(hEfficiency[15], "B_{s} = 20 MeV, #Gamma = 150 MeV", "elp");
    myLegend02a->Draw();

    myCanvas02a->Print("output/plots/hEfficiency_pl.png","png");
    myCanvas02a->Print("output/plots/hEfficiency_pl.eps","eps");

    //
    TCanvas* myCanvas03 = new TCanvas;

    //hEfficiencySystErr[0]->SetTitle("Efficiency");
    hEfficiencySystErr[0]->GetXaxis()->SetTitle("excess energy, MeV");
    hEfficiencySystErr[0]->GetXaxis()->SetTitleSize(0.06);
    hEfficiencySystErr[0]->GetXaxis()->SetTitleOffset(0.95);
    hEfficiencySystErr[0]->GetXaxis()->SetLabelSize(0.05);
    //hEfficiencySystErr[0]->GetXaxis()->SetRangeUser(-70.,30.);
    hEfficiencySystErr[0]->GetYaxis()->SetTitle("efficiency");
    hEfficiencySystErr[0]->GetYaxis()->SetTitleSize(0.06);
    hEfficiencySystErr[0]->GetYaxis()->SetTitleOffset(1.);
    hEfficiencySystErr[0]->GetYaxis()->SetLabelSize(0.05);
    hEfficiencySystErr[0]->GetYaxis()->SetRangeUser(0.,0.25);

    hEfficiencySystErr[0]->SetLineWidth(1);
    hEfficiencySystErr[0]->SetLineColor(1);
    hEfficiencySystErr[0]->Draw("E1");

    hEfficiencySystErr[1]->SetLineWidth(1);
    hEfficiencySystErr[1]->SetLineColor(51);
    hEfficiencySystErr[1]->Draw("same E1");

    hEfficiencySystErr[2]->SetLineWidth(1);
    hEfficiencySystErr[2]->SetLineColor(51);
    hEfficiencySystErr[2]->Draw("same E1");

    hEfficiencySystErr[3]->SetLineWidth(1);
    hEfficiencySystErr[3]->SetLineColor(61);
    hEfficiencySystErr[3]->Draw("same E1");

    hEfficiencySystErr[4]->SetLineWidth(1);
    hEfficiencySystErr[4]->SetLineColor(61);
    hEfficiencySystErr[4]->Draw("same E1");

    hEfficiencySystErr[5]->SetLineWidth(1);
    hEfficiencySystErr[5]->SetLineColor(68);
    hEfficiencySystErr[5]->Draw("same E1");

    hEfficiencySystErr[6]->SetLineWidth(1);
    hEfficiencySystErr[6]->SetLineColor(68);
    hEfficiencySystErr[6]->Draw("same E1");

    hEfficiencySystErr[7]->SetLineWidth(1);
    hEfficiencySystErr[7]->SetLineColor(94);
    hEfficiencySystErr[7]->Draw("same E1");

    hEfficiencySystErr[8]->SetLineWidth(1);
    hEfficiencySystErr[8]->SetLineColor(94);
    hEfficiencySystErr[8]->Draw("same E1");

    hEfficiencySystErr[9]->SetLineWidth(1);
    hEfficiencySystErr[9]->SetLineColor(98);
    hEfficiencySystErr[9]->Draw("same E1");

    hEfficiencySystErr[10]->SetLineWidth(1);
    hEfficiencySystErr[10]->SetLineColor(98);
    hEfficiencySystErr[10]->Draw("same E1");

    TLegend *myLegend03 = new TLegend(0.460, 0.570, 0.885, 0.885);
    myLegend03->SetFillStyle(1); myLegend03->SetFillColor(0); myLegend03->SetLineColor(0); myLegend03->SetTextSize(0.04);
    myLegend03->AddEntry(hEfficiencySystErr[0], "basic", "elp");
    myLegend03->AddEntry(hEfficiencySystErr[1], "#DeltaE(PSB)-#DeltaE(SEC)", "elp");
    myLegend03->AddEntry(hEfficiencySystErr[3], "invariant mass", "elp");
    myLegend03->AddEntry(hEfficiencySystErr[5], "#theta_{#pi^{0}-p}", "elp");
    myLegend03->AddEntry(hEfficiencySystErr[7], "missing mass", "elp");
    myLegend03->AddEntry(hEfficiencySystErr[9], "deuteron momentum", "elp");
    myLegend03->Draw("same");

    myCanvas03->Print("output/plots/hEfficiencySystErr.png","png");
    myCanvas03->Print("output/plots/hEfficiencySystErr.eps","eps");

    //
    TCanvas* myCanvas03a = new TCanvas;

    hEfficiencySystErr[0]->GetXaxis()->SetTitle("\\hbox{energia dostępna, MeV}");
    hEfficiencySystErr[0]->GetYaxis()->SetTitle("\\hbox{wydajność,%}");
    hEfficiencySystErr[0]->Draw("E1");

    hEfficiencySystErr[1]->Draw("same E1");
    hEfficiencySystErr[2]->Draw("same E1");
    hEfficiencySystErr[3]->Draw("same E1");
    hEfficiencySystErr[4]->Draw("same E1");
    hEfficiencySystErr[5]->Draw("same E1");
    hEfficiencySystErr[6]->Draw("same E1");
    hEfficiencySystErr[7]->Draw("same E1");
    hEfficiencySystErr[8]->Draw("same E1");
    hEfficiencySystErr[9]->Draw("same E1");
    hEfficiencySystErr[10]->Draw("same E1");

    TLegend *myLegend03a = new TLegend(0.530, 0.500, 0.885, 0.885);
    myLegend03a->SetFillStyle(1001); myLegend03a->SetFillColor(19); myLegend03a->SetLineColor(1); myLegend03a->SetTextSize(0.04); myLegend03a->SetBorderSize(5);
    myLegend03a->AddEntry(hEfficiencySystErr[0], "bazowa", "elp");
    myLegend03a->AddEntry(hEfficiencySystErr[1], "#DeltaE(PSB)-#DeltaE(SEC)", "elp");
    myLegend03a->AddEntry(hEfficiencySystErr[3], "masa niezmiennicza", "elp");
    myLegend03a->AddEntry(hEfficiencySystErr[5], "#theta_{#pi^{0}-p}", "elp");
    myLegend03a->AddEntry(hEfficiencySystErr[7], "\\hbox{masa brakująca}", "elp");
    myLegend03a->AddEntry(hEfficiencySystErr[9], "\\hbox{pęd deuteronu}", "elp");
    myLegend03a->Draw();

    myCanvas03a->Print("output/plots/hEfficiencySystErr_pl.png","png");
    myCanvas03a->Print("output/plots/hEfficiencySystErr_pl.eps","eps");

}
