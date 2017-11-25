#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TList.h"
#include "TLegendEntry.h"
#include "string.h"
#include <iostream>
#include "TTreeReader.h"
#include "TChain.h"
#include "TProfile.h"
#include <cmath>
#include "TProfile2D.h"
#include "TH2D.h"
#include "TArrayD.h"
#include "TPaveStats.h"
#include "TPad.h"

TH1F* defineTH1F(const char *file="",const char *branch=""){
  TFile *f = new TFile(file);
  TH1F *h = (TH1F*)f->Get(branch);
  return h;
}
TH2F* defineTH2F(const char *file="",const char *branch=""){
  TFile *f = new TFile(file);
  TH2F *h = (TH2F*)f->Get(branch);
  return h;
}
TPaveStats* makeStatBox(TH1F* h,double X1,double Y1,double X2,double Y2,int color){
  gROOT->SetBatch(kTRUE);
  TCanvas *c = new TCanvas();
  h->Draw();
  gPad->Update();
  TPaveStats *tps = (TPaveStats*) h->FindObject("stats");
  tps->SetTextColor(color);
  tps->SetLineColor(color);
  tps->SetX1NDC(X1);
  tps->SetX2NDC(X2);
  tps->SetY1NDC(Y1-(Y2-Y1));
  tps->SetY2NDC(Y1);
  gROOT->SetBatch(kFALSE);
  return tps;
}

void plot_RecHit(){

  ////////////////////////////////////////////////////////////////////
  //                 Define everything that is needed               //
  ////////////////////////////////////////////////////////////////////

  //Define/Get the Histos to be plotted

  //Single Pion
  TH1F* h1 = defineTH1F("206_plots/203_oldPhase/Output_Histo_RecHit.root","chi2_energyHEQIE111D");
  TH1F* h2 = defineTH1F("206_plots/NewData_206shifted_Y11shifted/Output_Histo_RecHit.root","chi2_energyHEQIE111D");
  TH1F* h3 = defineTH1F("206_plots/NewData_Y11shifted_205/Output_Histo_RecHit.root","chi2_energyHEQIE111D");
  TH1F* h4 = defineTH1F("206_plots/203_oldPhase/Output_Histo_RecHit.root","chi2_energyHEQIE111D_noCut");
  TH1F* h5 = defineTH1F("206_plots/NewData_206shifted_Y11shifted/Output_Histo_RecHit.root","chi2_energyHEQIE111D_noCut");
  TH1F* h6 = defineTH1F("206_plots/NewData_Y11shifted_205/Output_Histo_RecHit.root","chi2_energyHEQIE111D_noCut");
  //TTbar 2017
  TH1F* h7  = defineTH1F("206_plots/testPR/beforePR/TTbar_2017/Output_Histo_RecHit.root","chi2_energyHEQIE111D_noCut");
  TH1F* h8  = defineTH1F("206_plots/testPR/afterPR/TTbar_2017/Output_Histo_RecHit.root","chi2_energyHEQIE111D_noCut");
  TH1F* h9  = defineTH1F("206_plots/testPR/beforePR/TTbar_2017/Output_Histo_RecHit.root","chi2_energyHEQIE111D");
  TH1F* h10 = defineTH1F("206_plots/testPR/afterPR/TTbar_2017/Output_Histo_RecHit.root","chi2_energyHEQIE111D");
  //TTbar 2023
  TH1F* h11 = defineTH1F("206_plots/testPR/beforePR/TTbar_2023/Output_Histo_RecHit.root","chi2_energyHB1D_noCut");
  TH1F* h12 = defineTH1F("206_plots/testPR/afterPR/TTbar_2023/Output_Histo_RecHit.root","chi2_energyHB1D_noCut");
  //TTbar 2018
  TH1F* h13 = defineTH1F("206_plots/testPR/beforePR/TTbar_2018/Output_Histo_RecHit.root","chi2_energyHE_All1D_noCut");
  TH1F* h14 = defineTH1F("206_plots/testPR/afterPR/TTbar_2018/Output_Histo_RecHit.root","chi2_energyHE_All1D_noCut");
  TH1F* h15 = defineTH1F("206_plots/testPR/beforePR/TTbar_2018/Output_Histo_RecHit.root","chi2_energyHE_All1D");
  TH1F* h16 = defineTH1F("206_plots/testPR/afterPR/TTbar_2018/Output_Histo_RecHit.root","chi2_energyHE_All1D");

  TH1F* h17  = defineTH1F("206_plots/testPR/beforePR/TTbar_2018/Output_Histo_RecHit.root","chi2_energyHEQIE111D_noCut");
  TH1F* h18  = defineTH1F("206_plots/testPR/afterPR/TTbar_2018/Output_Histo_RecHit.root","chi2_energyHEQIE111D_noCut");
  TH1F* h19  = defineTH1F("206_plots/testPR/beforePR/TTbar_2018/Output_Histo_RecHit.root","chi2_energyHEQIE111D");
  TH1F* h20  = defineTH1F("206_plots/testPR/afterPR/TTbar_2018/Output_Histo_RecHit.root","chi2_energyHEQIE111D");
  h1->SetName("Sim 203, Reco 203");
  h2->SetName("Sim 206, Reco 206");
  h3->SetName("Sim 206, Reco 205");
  h4->SetName("Sim 203, Reco 203");
  h5->SetName("Sim 206, Reco 206");
  h6->SetName("Sim 206, Reco 205");
  h7->SetName("Sim 203, Reco 203");
  h8->SetName("Sim 206, Reco 206");
  h11->SetName("Sim 203, Reco 203");
  h12->SetName("Sim 206, Reco 206");
  h13->SetName("Sim 203, Reco 203");
  h14->SetName("Sim 206, Reco 206");
  h15->SetName("Sim 203, Reco 203");
  h16->SetName("Sim 206, Reco 206");
  h17->SetName("Sim 203, Reco 203");
  h18->SetName("Sim 206, Reco 206");
  h19->SetName("Sim 203, Reco 203");
  h20->SetName("Sim 206, Reco 206");
  
  double x1,y1,x2,y2;
  x1=0.78; y1=0.94; x2=0.98; y2=1.09;
  TPaveStats *tps1 = makeStatBox(h1,x1,y1,x2,y2,kBlack);
  double X1 = tps1->GetX1NDC(); double Y1 = tps1->GetY1NDC(); double X2 = tps1->GetX2NDC(); double Y2 = tps1->GetY2NDC();
  TPaveStats *tps2 = makeStatBox(h2,X1,Y1,X2,Y2,kRed);
  double X3 = tps2->GetX1NDC(); double Y3 = tps2->GetY1NDC(); double X4 = tps2->GetX2NDC(); double Y4 = tps2->GetY2NDC();
  TPaveStats *tps3 = makeStatBox(h3,X3,Y3,X4,Y4,kBlue);
  TPaveStats *tps4 = makeStatBox(h4,x1,y1,x2,y2,kBlack);
  TPaveStats *tps5 = makeStatBox(h5,X1,Y1,X2,Y2,kRed);
  TPaveStats *tps6 = makeStatBox(h6,X3,Y3,X4,Y4,kBlue);
  TPaveStats *tps7 = makeStatBox(h7,x1,y1,x2,y2,kBlack);
  TPaveStats *tps8 = makeStatBox(h8,X1,Y1,X2,Y2,kRed);
  TPaveStats *tps11 = makeStatBox(h11,x1,y1,x2,y2,kBlack);
  TPaveStats *tps12 = makeStatBox(h12,X1,Y1,X2,Y2,kRed);
  TPaveStats *tps13 = makeStatBox(h13,x1,y1,x2,y2,kBlack);
  TPaveStats *tps14 = makeStatBox(h14,X1,Y1,X2,Y2,kRed);
  TPaveStats *tps15 = makeStatBox(h15,x1,y1,x2,y2,kBlack);
  TPaveStats *tps16 = makeStatBox(h16,X1,Y1,X2,Y2,kRed);
  TPaveStats *tps17 = makeStatBox(h17,x1,y1,x2,y2,kBlack);
  TPaveStats *tps18 = makeStatBox(h18,X1,Y1,X2,Y2,kRed);
  TPaveStats *tps19 = makeStatBox(h19,x1,y1,x2,y2,kBlack);
  TPaveStats *tps20 = makeStatBox(h20,X1,Y1,X2,Y2,kRed);

  //2017 Single Pion Gun (With Cuts)
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.13);

  TLegend* catLeg1 = new TLegend(0.17,0.7,0.4,0.9);
  h1->Draw("hist");
  h1->SetTitle("2017 Single Pion 500GeV: HEP17 ");
  h1->GetXaxis()->SetTitle("Log10(chi2)");
  h1->GetYaxis()->SetTitle("RecHits > 5 GeV");
  h1->SetTitleSize(0.001);
  h1->SetTitleSize(0.05,"Y");
  h1->SetTitleSize(0.05,"X");
  h1->SetTitleOffset(1.0,"Y");
  h1->SetTitleOffset(1.0,"X");
  h1->SetLabelSize(0.03,"Y");
  h1->SetLabelSize(0.03,"X");
  h1->SetLineColor(kBlack);
  h2->SetLineColor(kRed);
  h2->Draw("hist same");
  h3->SetLineColor(kBlue);
  h3->Draw("hist same");
  tps1->Draw("same");
  tps2->Draw("same");
  tps3->Draw("same");
  catLeg1->SetTextSize(0.03);
  catLeg1->AddEntry(h1,"Sim 203, Reco 203","l");
  catLeg1->AddEntry(h2,"Sim 206, Reco 206","l");
  catLeg1->AddEntry(h3,"Sim 206, Reco 205","l");  
  //catLeg1->Draw();

  //2017 Single Pion Gun (Without Cuts)
  TCanvas *c2 = new TCanvas("c2","c2",1200,800);
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.13);

  TLegend* catLeg2 = new TLegend(0.17,0.7,0.4,0.9);
  h4->Draw("hist");
  h4->SetTitle("2017 Single Pion 500GeV: HEP17");
  h4->GetXaxis()->SetTitle("Log10(chi2)");
  h4->GetYaxis()->SetTitle("RecHits");
  h4->SetTitleSize(0.001);
  h4->SetTitleSize(0.05,"Y");
  h4->SetTitleSize(0.05,"X");
  h4->SetTitleOffset(1.0,"Y");
  h4->SetTitleOffset(1.0,"X");
  h4->SetLabelSize(0.03,"Y");
  h4->SetLabelSize(0.03,"X");
  h4->SetLineColor(kBlack);
  h5->SetLineColor(kRed);
  h5->Draw("hist same");
  h6->SetLineColor(kBlue);
  h6->Draw("hist same");
  tps4->Draw("same");
  tps5->Draw("same");
  tps6->Draw("same");
  catLeg2->SetTextSize(0.03);
  catLeg2->AddEntry(h4,"Sim 203, Reco 203","l");
  catLeg2->AddEntry(h5,"Sim 206, Reco 206","l");
  catLeg2->AddEntry(h6,"Sim 206, Reco 205","l");  
  //catLeg2->Draw();

  //2017 TTBar (Without Cuts)
  TCanvas *c3 = new TCanvas("c3","c3",1200,800);
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.13);

  TLegend* catLeg3 = new TLegend(0.17,0.7,0.4,0.9);
  h7->Draw("hist");
  h7->SetTitle("2017 TTBar 13TeV: HEP17");
  h7->GetXaxis()->SetTitle("Log10(chi2)");
  h7->GetYaxis()->SetTitle("RecHits");
  h7->SetTitleSize(0.001);
  h7->SetTitleSize(0.05,"Y");
  h7->SetTitleSize(0.05,"X");
  h7->SetTitleOffset(1.0,"Y");
  h7->SetTitleOffset(1.0,"X");
  h7->SetLabelSize(0.03,"Y");
  h7->SetLabelSize(0.03,"X");
  h7->SetLineColor(kBlack);
  h8->SetLineColor(kRed);
  h8->Draw("hist same");
  tps7->Draw("same");
  tps8->Draw("same");
  catLeg3->SetTextSize(0.03);
  catLeg3->AddEntry(h7,"Sim 203, Reco 203","l");
  catLeg3->AddEntry(h8,"Sim 206, Reco 206","l");
  //catLeg3->Draw();

  //2023 TTBar (Without Cuts)
  TCanvas *c4 = new TCanvas("c4","c4",1200,800);
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.13);

  TLegend* catLeg4 = new TLegend(0.17,0.7,0.4,0.9);
  h11->Draw("hist");
  h11->SetTitle("2023 TTBar 14TeV: HB");
  h11->GetXaxis()->SetTitle("Log10(chi2)");
  h11->GetYaxis()->SetTitle("RecHits");
  h11->SetTitleSize(0.001);
  h11->SetTitleSize(0.05,"Y");
  h11->SetTitleSize(0.05,"X");
  h11->SetTitleOffset(1.0,"Y");
  h11->SetTitleOffset(1.0,"X");
  h11->SetLabelSize(0.03,"Y");
  h11->SetLabelSize(0.03,"X");
  h11->SetLineColor(kBlack);
  h12->SetLineColor(kRed);
  h12->Draw("hist same");
  tps11->Draw("same");
  tps12->Draw("same");
  catLeg4->SetTextSize(0.03);
  catLeg4->AddEntry(h11,"Sim 203, Reco 203","l");
  catLeg4->AddEntry(h12,"Sim 206, Reco 206","l");
  //catLeg4->Draw();
  
  //2018 TTBar (Without Cuts)
  TCanvas *c5 = new TCanvas("c5","c5",1200,800);
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.13);

  TLegend* catLeg5 = new TLegend(0.17,0.7,0.4,0.9);
  h13->Draw("hist");
  h13->SetTitle("2018 TTBar 13TeV: HE");
  h13->GetXaxis()->SetTitle("Log10(chi2)");
  h13->GetYaxis()->SetTitle("RecHits");
  h13->SetTitleSize(0.001);
  h13->SetTitleSize(0.05,"Y");
  h13->SetTitleSize(0.05,"X");
  h13->SetTitleOffset(1.0,"Y");
  h13->SetTitleOffset(1.0,"X");
  h13->SetLabelSize(0.03,"Y");
  h13->SetLabelSize(0.03,"X");
  h13->SetLineColor(kBlack);
  h14->SetLineColor(kRed);
  h14->Draw("hist same");
  tps13->Draw("same");
  tps14->Draw("same");
  catLeg5->SetTextSize(0.03);
  catLeg5->AddEntry(h13,"Sim 203, Reco 203","l");
  catLeg5->AddEntry(h14,"Sim 206, Reco 206","l");
  //catLeg5->Draw();

  //2018 TTBar (With Cuts)
  TCanvas *c6 = new TCanvas("c6","c6",1200,800);
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.13);

  TLegend* catLeg6 = new TLegend(0.17,0.7,0.4,0.9);
  h16->Draw("hist");
  h16->SetTitle("2018 TTBar 13TeV: HE");
  h16->GetXaxis()->SetTitle("Log10(chi2)");
  h16->GetYaxis()->SetTitle("RecHits > 5 GeV");
  h16->SetTitleSize(0.001);
  h16->SetTitleSize(0.05,"Y");
  h16->SetTitleSize(0.05,"X");
  h16->SetTitleOffset(1.0,"Y");
  h16->SetTitleOffset(1.0,"X");
  h16->SetLabelSize(0.03,"Y");
  h16->SetLabelSize(0.03,"X");
  h15->SetLineColor(kBlack);
  h16->SetLineColor(kRed);
  h15->Draw("hist same");
  tps15->Draw("same");
  tps16->Draw("same");
  catLeg6->SetTextSize(0.03);
  catLeg6->AddEntry(h15,"Sim 203, Reco 203","l");
  catLeg6->AddEntry(h16,"Sim 206, Reco 206","l");
  //catLeg6->Draw();

  //2018 TTBar HEP17 (Without Cuts)
  TCanvas *c7 = new TCanvas("c7","c7",1200,800);
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.13);

  TLegend* catLeg7 = new TLegend(0.17,0.7,0.4,0.9);
  h17->Draw("hist");
  h17->SetTitle("2018 TTBar 13TeV: HEP17");
  h17->GetXaxis()->SetTitle("Log10(chi2)");
  h17->GetYaxis()->SetTitle("RecHits");
  h17->SetTitleSize(0.001);
  h17->SetTitleSize(0.05,"Y");
  h17->SetTitleSize(0.05,"X");
  h17->SetTitleOffset(1.0,"Y");
  h17->SetTitleOffset(1.0,"X");
  h17->SetLabelSize(0.03,"Y");
  h17->SetLabelSize(0.03,"X");
  h17->SetLineColor(kBlack);
  h18->SetLineColor(kRed);
  h18->Draw("hist same");
  tps17->Draw("same");
  tps18->Draw("same");
  catLeg7->SetTextSize(0.03);
  catLeg7->AddEntry(h17,"Sim 203, Reco 203","l");
  catLeg7->AddEntry(h18,"Sim 206, Reco 206","l");
  //catLeg7->Draw();
  
  //2018 TTBar HEP17 (With Cuts)
  TCanvas *c8 = new TCanvas("c8","c8",1200,800);
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.13);

  TLegend* catLeg8 = new TLegend(0.17,0.7,0.4,0.9);
  h20->Draw("hist");
  h20->SetTitle("2018 TTBar 13TeV: HEP17");
  h20->GetXaxis()->SetTitle("Log10(chi2)");
  h20->GetYaxis()->SetTitle("RecHits > 5 GeV");
  h20->SetTitleSize(0.001);
  h20->SetTitleSize(0.05,"Y");
  h20->SetTitleSize(0.05,"X");
  h20->SetTitleOffset(1.0,"Y");
  h20->SetTitleOffset(1.0,"X");
  h20->SetLabelSize(0.03,"Y");
  h20->SetLabelSize(0.03,"X");
  h19->SetLineColor(kBlack);
  h20->SetLineColor(kRed);
  h19->Draw("hist same");
  tps19->Draw("same");
  tps20->Draw("same");
  catLeg8->SetTextSize(0.03);
  catLeg8->AddEntry(h19,"Sim 203, Reco 203","l");
  catLeg8->AddEntry(h20,"Sim 206, Reco 206","l");
  //catLeg8->Draw();

  //Save Plot
  c1->SaveAs("All_Chi2_1D.pdf");
  c2->SaveAs("All_Chi2_1D_noCut.pdf");
  c3->SaveAs("All_Chi2_1D_noCut_TTbar2017.pdf");
  c4->SaveAs("All_Chi2_1D_noCut_TTbar2023.pdf");
  c5->SaveAs("All_Chi2_1D_noCut_TTbar2018.pdf");
  c6->SaveAs("All_Chi2_1D_TTbar2018.pdf");
  c7->SaveAs("All_Chi2_1D_noCut_TTbar2018_HEP17.pdf");
  c8->SaveAs("All_Chi2_1D_TTbar2018_HEP17.pdf");

}//End of Macro
