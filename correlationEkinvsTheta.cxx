#include<TH1F.h>
#include<TH2F.h>
#include<TH3F.h>
#include<TVector3.h>
#include<TLorentzVector.h>
#include<TF1.h>
#include<TFile.h>
#include<TTree.h>
#include<TMath.h>
#include<TCanvas.h>
#include<TClonesArray.h>
#include<TPaveLabel.h>
#include<TFrame.h>
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
#include <vector>

void correlationEkinvsTheta() {

    TFile* myFile = new TFile("input/MC-acceptance.root","READ");

    TH2F* hist_EkinvsTheta[3];

    myFile->cd("Histograms");

    hist_EkinvsTheta[0] = (TH2F*)gDirectory->Get("WMC/hEkin_vs_Theta_d_lab_MC");
    hist_EkinvsTheta[1] = (TH2F*)gDirectory->Get("WMC/hEkin_vs_Theta_p_lab_MC");
    hist_EkinvsTheta[2] = (TH2F*)gDirectory->Get("WMC/hEkin_vs_Theta_g1_lab_MC");

    ////
    gStyle->SetOptStat(kFALSE);
    gStyle->SetPalette(1,0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.16);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPalette(55);

    gStyle->SetTitleFont(62,"XYZ");
    gStyle->SetLabelFont(62,"XYZ");
    gStyle->SetTextFont(62);

    //
    TCanvas* MyCanvas00 = new TCanvas;

    hist_EkinvsTheta[0]->GetXaxis()->SetTitle("kinetic energy, GeV");
    hist_EkinvsTheta[0]->GetYaxis()->SetTitle("#theta,#circ");
    hist_EkinvsTheta[0]->GetXaxis()->SetTitleSize(0.06);
    hist_EkinvsTheta[0]->GetXaxis()->SetTitleOffset(1.0);
    hist_EkinvsTheta[0]->GetXaxis()->SetLabelSize(0.05);
    hist_EkinvsTheta[0]->GetYaxis()->SetTitleSize(0.06);
    hist_EkinvsTheta[0]->GetYaxis()->SetTitleOffset(1.05);
    hist_EkinvsTheta[0]->GetYaxis()->SetLabelSize(0.05);
    hist_EkinvsTheta[0]->GetXaxis()->SetRangeUser(0.,0.5);
    hist_EkinvsTheta[0]->GetYaxis()->SetRangeUser(0.,20.);
    hist_EkinvsTheta[0]->GetZaxis()->SetLabelSize(0.05);
    hist_EkinvsTheta[0]->GetXaxis()->SetNdivisions(5,5,0, kTRUE);
    //hist_EkinvsTheta[0]->SetMaximum(8000.);
    hist_EkinvsTheta[0]->Draw("colz");

    TLine* line000 = new TLine(0.,3.,0.5,3.);
    line000->SetLineColor(2);
    line000->SetLineWidth(2);
    line000->SetLineStyle(1);
    line000->Draw("same");

    TLine* line001 = new TLine(0.,18.,0.5,18.);
    line001->SetLineColor(2);
    line001->SetLineWidth(2);
    line001->SetLineStyle(1);
    line001->Draw("same");

    TPaveText *detec00 = new TPaveText(0.47,10.,0.47,10.,"FD");
    detec00->SetTextSize(0.05);
    detec00->SetFillColor(0);
    detec00->SetTextColor(2);
    detec00->SetTextAlign(22);
    detec00->AddText("FD");
    detec00->Draw("same");

    MyCanvas00->Print("output/plots/hEkinvsTheta_deuteron.eps","eps");
    MyCanvas00->Print("output/plots/hEkinvsTheta_deuteron.png","png");

    //
    TCanvas* MyCanvas00a = new TCanvas("MyCanvas00a","",465,500);
    hist_EkinvsTheta[0]->GetXaxis()->SetTitle("energia kinetyczna, GeV");
    hist_EkinvsTheta[0]->Draw("colz");
    line000->Draw("same");
    line001->Draw("same");
    detec00->Draw("same");
    MyCanvas00a->Print("output/plots/hEkinvsTheta_deuteron_pl.eps","eps");
    MyCanvas00a->Print("output/plots/hEkinvsTheta_deuteron_pl.png","png");

    //
    TCanvas* MyCanvas01 = new TCanvas;

    hist_EkinvsTheta[1]->GetXaxis()->SetTitle("kinetic energy, GeV");
    hist_EkinvsTheta[1]->GetYaxis()->SetTitle("#theta,#circ");
    hist_EkinvsTheta[1]->GetXaxis()->SetTitleSize(0.06);
    hist_EkinvsTheta[1]->GetXaxis()->SetTitleOffset(1.0);
    hist_EkinvsTheta[1]->GetXaxis()->SetLabelSize(0.05);
    hist_EkinvsTheta[1]->GetYaxis()->SetTitleSize(0.06);
    hist_EkinvsTheta[1]->GetYaxis()->SetTitleOffset(1.05);
    hist_EkinvsTheta[1]->GetYaxis()->SetLabelSize(0.05);
    hist_EkinvsTheta[1]->GetXaxis()->SetRangeUser(0.,0.70);
    //hist_EkinvsTheta[1]->GetYaxis()->SetRangeUser(0.,180.);
    hist_EkinvsTheta[1]->GetZaxis()->SetLabelSize(0.05);
    //hist_EkinvsTheta[1]->SetMaximum(8000.);
    hist_EkinvsTheta[1]->Draw("colz");

    TLine* line010 = new TLine(0.,3.,0.7,3.);
    line010->SetLineColor(2);
    line010->SetLineWidth(2);
    line010->SetLineStyle(1);
    line010->Draw("same");

    TLine* line011 = new TLine(0.,18.,0.7,18.);
    line011->SetLineColor(2);
    line011->SetLineWidth(2);
    line011->SetLineStyle(1);
    line011->Draw("same");

    TLine* line012 = new TLine(0.,20.,0.7,20.);
    line012->SetLineColor(2);
    line012->SetLineWidth(2);
    line012->SetLineStyle(1);
    line012->Draw("same");

    TLine* line013 = new TLine(0.,169.,0.7,169.);
    line013->SetLineColor(2);
    line013->SetLineWidth(2);
    line013->SetLineStyle(1);
    line013->Draw("same");

    TPaveText *detec010 = new TPaveText(0.65,95.,0.65,95.,"CD");
    detec010->SetTextSize(0.05);
    detec010->SetFillColor(0);
    detec010->SetTextColor(4);
    detec010->SetTextAlign(22);
    detec010->AddText("CD");
    detec010->Draw("same");

    TPaveText *detec011 = new TPaveText(0.65,11.,0.65,11.,"FD");
    detec011->SetTextSize(0.05);
    detec011->SetFillColor(0);
    detec011->SetTextColor(2);
    detec011->SetTextAlign(22);
    detec011->AddText("FD");
    detec011->Draw("same");

    MyCanvas01->Print("output/plots/hEkinvsTheta_proton.eps","eps");
    MyCanvas01->Print("output/plots/hEkinvsTheta_proton.png","png");

    //
    TCanvas* MyCanvas01a = new TCanvas("MyCanvas01a","",465,500);
    hist_EkinvsTheta[1]->GetXaxis()->SetTitle("energia kinetyczna, GeV");
    hist_EkinvsTheta[1]->Draw("colz");
    line010->Draw("same");
    line011->Draw("same");
    line012->Draw("same");
    line013->Draw("same");
    detec010->Draw("same");
    detec011->Draw("same");
    MyCanvas01a->Print("output/plots/hEkinvsTheta_proton_pl.eps","eps");
    MyCanvas01a->Print("output/plots/hEkinvsTheta_proton_pl.png","png");

    //
    TCanvas* MyCanvas02 = new TCanvas;
    hist_EkinvsTheta[2]->GetXaxis()->SetTitle("kinetic energy, GeV");
    hist_EkinvsTheta[2]->GetYaxis()->SetTitle("#theta,#circ");
    hist_EkinvsTheta[2]->GetXaxis()->SetTitleSize(0.06);
    hist_EkinvsTheta[2]->GetXaxis()->SetTitleOffset(1.0);
    hist_EkinvsTheta[2]->GetXaxis()->SetLabelSize(0.05);
    hist_EkinvsTheta[2]->GetYaxis()->SetTitleSize(0.06);
    hist_EkinvsTheta[2]->GetYaxis()->SetTitleOffset(1.05);
    hist_EkinvsTheta[2]->GetYaxis()->SetLabelSize(0.05);
    hist_EkinvsTheta[2]->GetXaxis()->SetRangeUser(0.,0.8);
    //hist_EkinvsTheta[2]->GetYaxis()->SetRangeUser(0.,180.);
    hist_EkinvsTheta[2]->GetZaxis()->SetLabelSize(0.05);
    //hist_EkinvsTheta[2]->SetMaximum(8000.);
    hist_EkinvsTheta[2]->Draw("colz");

    TLine* line020 = new TLine(0.,20.,0.8,20.);
    line020->SetLineColor(2);
    line020->SetLineWidth(2);
    line020->SetLineStyle(1);
    line020->Draw("same");

    TLine* line021 = new TLine(0.,169.,0.8,169.);
    line021->SetLineColor(2);
    line021->SetLineWidth(2);
    line021->SetLineStyle(1);
    line021->Draw("same");

    TPaveText *detec02 = new TPaveText(0.75,95.,0.75,95.,"CD");
    detec02->SetTextSize(0.05);
    detec02->SetFillColor(0);
    detec02->SetTextColor(4);
    detec02->SetTextAlign(22);
    detec02->AddText("CD");
    detec02->Draw("same");

    MyCanvas02->Print("output/plots/hEkinvsTheta_gamma.eps","eps");
    MyCanvas02->Print("output/plots/hEkinvsTheta_gamma.png","png");

    //
    TCanvas* MyCanvas02a = new TCanvas("MyCanvas02a","",465,500);
    hist_EkinvsTheta[2]->GetXaxis()->SetTitle("energia kinetyczna, GeV");
    hist_EkinvsTheta[2]->Draw("colz");
    line020->Draw("same");
    line021->Draw("same");
    detec02->Draw("same");
    MyCanvas02a->Print("output/plots/hEkinvsTheta_gamma_pl.eps","eps");
    MyCanvas02a->Print("output/plots/hEkinvsTheta_gamma_pl.png","png");

}

