#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TDirectory.h>
//pulse shape vectors
#include "pulses.h"

//Y11 pulse shape vectors
#include "dataMip.h"

#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TF1.h"

int shiftMC  = 74; //shift MC Pulse
int phase    = 0;//-81;//50; //Change phase of 203new/old and 205
int shift    = 0;  //shift components
int shiftpulses = 1; //Shift just 203pulses and 205pulses
double nor = 1.0;  //Make 203s bigger in plot

static const int nBinsSiPM_ = 251;
std::vector<float> siPMShape2017_(nBinsSiPM_,0.0);

//this function can take function pointers *or* functors!
template <class F1, class F2>
static std::vector<double> convolve(unsigned nbin, F1 f1, F2 f2){
  std::vector<double> result(2*nbin-1,0.);
  for(unsigned i = 0; i < 2*nbin-1; ++i){
    for(unsigned j = 0; j < std::min(i+1,nbin); ++j){
      double tmp = f1(j)*f2(i-j);
      if(std::isnan(tmp) or std::isinf(tmp)) continue;
      result[i] += tmp;
      //if(i==250)std::cout<<j<<" SiPM = "<<f1(j)<<"  Y11 = "<< f2(j)<<endl;
    }
  }
  return result;
}

TGraph* makeGraph(int num, std::vector<std::vector<double>>& data){
  Double_t x[num], y[num];
  for(unsigned i = 0; i < num; ++i){
    x[i] = data[i][0];
    y[i] = data[i][1];	
  }
  TGraph* gtemp = new TGraph(num,x,y);
  gtemp->SetMarkerStyle(kFullSquare);
  return gtemp;
}

TH1F* makeLandau(double sigma,std::string name){
  TH1F* gtemp = new TH1F(name.c_str(),"",252,0,251);
  for(int i = 1; i <= 252; i++){
    gtemp->SetBinContent(i,TMath::Landau(i,89+shiftpulses+phase,sigma));
  }
  gtemp->Scale(1.0/gtemp->Integral(0,251));
  return gtemp;
}

double frac = 0.11;
double sh = 7.2;

double fitY11results(double t){
  //Fit From Deconvolved Data
  double A,n,t0;
  A=0.104204; n=0.44064; t0=10.0186;
  if(t>sh) return A*(1-exp(-(t-sh)/n))*exp(-(t-sh)/t0);
  else return 0.0;
}

double corTerm(double t){
  double norm,mpv,sigma;
  //norm=0.0806123; //Unshifted norm
  //norm=0.0809721; //shift of 6.9ns
  //norm=0.0809775; //shift of 7ns
  norm=0.0809882; //shift of 7.2ns
  //norm=0.0810045; //shift of 7.5ns
  //norm=0.0810316; //shift of 8ns
  mpv=0; sigma=20;
  if(t>sh) return norm*TMath::Landau((t-sh),mpv,sigma);  
  else return 0.0;
}

double Y11TimePDF(double t) {
  double val = (1-frac)*fitY11results(t) + frac*corTerm(t);
  if(val >= 0) return val;
  else return 0.0;
}

double onePulse(double t, double A, double sigma, double theta, double m) {
  return (t<theta) ? 0 : A*TMath::LogNormal(t,sigma,theta,m);
}

double analyticPulseShapeSiPMHE(double t) {
  // taken from fit to laser measurement taken by Iouri M. in Spring 2016.
  double A1(5.204/6.94419), sigma1_shape(0.5387), theta1_loc(-0.3976), m1_scale(4.428);
  double A2(1.855/6.94419), sigma2_shape(0.8132), theta2_loc(7.025),   m2_scale(12.29);
  return
    onePulse(t,A1,sigma1_shape,theta1_loc,m1_scale) +
    onePulse(t,A2,sigma2_shape,theta2_loc,m2_scale);
}

TH1F* computeSiPMShape2017(string name){
  //numerical convolution of SiPM pulse + WLS fiber shape
  std::vector<double> nt = convolve(nBinsSiPM_,analyticPulseShapeSiPMHE,Y11TimePDF);

  double norm = 0.;
  for (unsigned int j = 1; j <= nBinsSiPM_; ++j) {
    norm += (nt[j]>0) ? nt[j] : 0.;
  }
  for (unsigned int j = 1; j <= nBinsSiPM_; ++j) {
    nt[j] /= norm;
    siPMShape2017_[j]=nt[j];
    //std::cout<<siPMShape2017_[j]<<","<<endl;
  }

  TH1F* htemp = new TH1F(name.c_str(),"",nBinsSiPM_+1,0.,nBinsSiPM_);
  for(unsigned int i = 1; i<=nBinsSiPM_+1; i++){
    htemp->SetBinContent(i+shiftMC+phase,siPMShape2017_[i]);
  }
  htemp->Scale(nor/htemp->Integral(0,htemp->GetNbinsX()+1));
  return htemp;
}

TH1F* makeY11Histo(string name){
  TH1F* htemp = new TH1F(name.c_str(),"",nBinsSiPM_+1,0.,nBinsSiPM_);
  for(unsigned int i = 1; i<=nBinsSiPM_+1; i++){
    htemp->SetBinContent(i+9,Y11TimePDF(i));
  }
  //std::cout<<"Norm of Y11: "<<htemp->Integral(0,htemp->GetNbinsX()+1)<<std::endl;
  //htemp->Scale(1.0/htemp->Integral(0,htemp->GetNbinsX()+1));
  htemp->Scale(1.0/htemp->GetBinContent(htemp->GetMaximumBin()));
  return htemp;
}

TH1F* makePulseHisto(vector<double>& pulse, int bstart, string name){
	//make pulse template histo (from convolution)
	TH1F* htempT = new TH1F(name.c_str(),"",252,0,251);
	for(unsigned b = 0; b < pulse.size(); ++b){
	  htempT->SetBinContent(b+bstart,pulse[b]);
	}
	htempT->Scale(nor/htempT->Integral(0,htempT->GetNbinsX()+1));
	return htempT;
}

void fineTuningY11Shape(){
//int main(){
  TCanvas *c0 = new TCanvas("PulseShape_203","",1000,800);  
  TLegend* catLeg0 = new TLegend(0.8,0.72,0.99,0.87);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.16);
  gPad->SetLogy();

  TH1F* h203_MC = NULL;
  h203_MC = computeSiPMShape2017("203_MC");
  h203_MC->SetMinimum(0.00000001);
  h203_MC->SetMaximum(0.05);
  //h203_MC->GetXaxis()->SetRange(0,100);
  h203_MC->GetXaxis()->SetTitle("time (ns)");
  h203_MC->GetYaxis()->SetTitle("A.U.");
  h203_MC->SetTitle("Pulse Shape Comparison");
  h203_MC->SetName("Pulse Shape  Comparison");
  h203_MC->SetTitleSize(0.002);
  h203_MC->SetTitleSize(0.05,"X");
  h203_MC->SetTitleSize(0.04,"Y");
  h203_MC->SetTitleOffset(1.2,"X");
  h203_MC->SetTitleOffset(1.8,"Y");
  h203_MC->SetLabelSize(0.05,"X");
  h203_MC->SetLabelSize(0.05,"Y");
  h203_MC->SetStats(false);
  h203_MC->SetLineColor(kRed);
  h203_MC->Draw("hist");
  //h203_MC->Print("all");
  catLeg0->SetBorderSize(0);
  catLeg0->SetFillStyle(0);
  catLeg0->SetTextSize(0.04);
  catLeg0->AddEntry(h203_MC,"206","l");
  catLeg0->Draw();
  
  TH1F* h205_pulses = NULL;
  h205_pulses = makePulseHisto(s205,69+shiftpulses+phase,"205");
  h205_pulses->SetLineColor(kGreen+2);
  h205_pulses->Draw("same hist");
  catLeg0->AddEntry(h205_pulses,"205","l");

  TH1F* new205 = NULL;
  new205 = makeLandau(4.57,"New 205");
  new205->Draw("hist same");
  //new205->SetLineColor(kBlue);
  //catLeg0->AddEntry(new205,"New 205");

  TH1F* h203_pulses = NULL;
  h203_pulses = makePulseHisto(s203,65+shiftpulses+phase,"203");
  h203_pulses->SetLineColor(kBlack);
  h203_pulses->Draw("same hist");
  catLeg0->AddEntry(h203_pulses,"203","l");

  TH1F* h206_pulses = NULL;
  h206_pulses = makePulseHisto(s206,76+shiftpulses+phase,"206 pulse");
  h206_pulses->SetLineColor(kGreen);
  //h206_pulses->Draw("same hist");
  //catLeg0->AddEntry(h206_pulses,"206 pulse","l");

  TH1F* h207_pulses = NULL;
  h207_pulses = makePulseHisto(s207,70+shiftpulses+phase,"207");
  h207_pulses->SetLineColor(kBlue);
  h207_pulses->Draw("same hist");
  catLeg0->AddEntry(h207_pulses,"207","l");

  TCanvas *c1 = new TCanvas("PulseShape_Y11","",1000,800);  
  TLegend* catLeg1 = new TLegend(0.55,0.65,0.99,0.9);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.16);
  //gPad->SetLogy();

  TH1F* Y11_MC = NULL;
  Y11_MC = makeY11Histo("Y11_MC");
  Y11_MC->SetMinimum(0.00000001);
  Y11_MC->SetMaximum(1.1);
  Y11_MC->GetXaxis()->SetRange(0,100);
  Y11_MC->GetXaxis()->SetTitle("time (ns)");
  Y11_MC->GetYaxis()->SetTitle("A.U.");
  Y11_MC->SetTitle("Y11 Shape Comparison");
  Y11_MC->SetName("");
  Y11_MC->SetTitleSize(0.002);
  Y11_MC->SetTitleSize(0.05,"X");
  Y11_MC->SetTitleSize(0.04,"Y");
  Y11_MC->SetTitleOffset(1.2,"X");
  Y11_MC->SetTitleOffset(1.8,"Y");
  Y11_MC->SetLabelSize(0.05,"X");
  Y11_MC->SetLabelSize(0.05,"Y");
  Y11_MC->SetStats(false);
  Y11_MC->SetLineColor(kBlack);
  Y11_MC->Draw("hist");
  catLeg1->SetBorderSize(0);
  catLeg1->SetFillStyle(0);
  catLeg1->SetTextSize(0.04);
  catLeg1->AddEntry(Y11_MC,"Y11+Landau","l");
  catLeg1->Draw();

  TGraph* DeConvData = NULL;
  DeConvData = makeGraph(17,deConvData);
  DeConvData->SetMarkerColor(kBlue);
  DeConvData->SetLineColor(kBlue);
  DeConvData->Draw("P same");
  
  TF1 *fa2 = new TF1("fa2","[0]*(1-exp(-(x-9)/0.44064))*exp(-(x-9)/10.0186) + [1]*TMath::Landau(x,[2]-9,[3])",0,250);
  fa2->SetParLimits(0,0,15);
  fa2->SetParLimits(1,0,15);
  fa2->SetParLimits(2,0,15);
  fa2->SetParLimits(3,5,15);
  fa2->SetParameter(0,5.84053);
  fa2->FixParameter(1,0.0);
  fa2->SetParameter(2,2);
  fa2->SetParameter(3,5);
  DeConvData->Fit(fa2,"NQ","",0,60);
  fa2->SetLineColor(kRed);
  catLeg1->AddEntry(fa2,"Y11: Fit Only","l");
  fa2->Draw("same");
  catLeg1->AddEntry(DeConvData,"Deconvolved Data","P");
  DeConvData->Draw("P same");
  Y11_MC->Draw("hist same");
  
  //std::cout<<"203Root:  "<<h203_MC->GetMaximumBin()<<endl;
  //std::cout<<"203pulses.h:  "<<h203_pulses->GetMaximumBin()<<endl;
  //std::cout<<"new205:  "<<new205->GetMaximumBin()<<endl;
  //std::cout<<"Y11 New:    "<<Y11_MC->GetMaximumBin()<<endl;
  
  TCanvas *c2 = new TCanvas("PulseShape_Y11_components","",1000,800);  
  TLegend* catLeg2 = new TLegend(0.6,0.65,0.99,0.9);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.16);
  gPad->SetLogy();

  TH1F* htemp = new TH1F("","",251,-10.,250);
  TF1 *Y11 = new TF1("Y11","Y11TimePDF(x)",-10,250);
  TF1 *fitY11 = new TF1("fitY11","0.89*fitY11results(x)",-10,250);
  TF1 *CorTerm = new TF1("CorTerm","0.11*corTerm(x)",-10,250);
  std::cout<<"norm of overallY11  "<<1/Y11->Integral(0,250)<<std::endl;
  std::cout<<"norm of fitY11 "<<1/fitY11->Integral(0,250)<<std::endl;
  std::cout<<"norm of Landau "<<1/CorTerm->Integral(0,250)<<std::endl;
  htemp->SetMinimum(0.000001);
  htemp->SetMaximum(0.08);
  //htemp->GetXaxis()->SetRange(0,100);  
  htemp->GetXaxis()->SetTitle("time (ns)");
  htemp->GetYaxis()->SetTitle("A.U.");
  htemp->SetTitle("Scintillator+Y11 Components");
  htemp->SetName("Shape Comparison");
  htemp->SetTitleSize(0.002);
  htemp->SetTitleSize(0.05,"X");
  htemp->SetTitleSize(0.04,"Y");
  htemp->SetTitleOffset(1.2,"X");
  htemp->SetTitleOffset(1.8,"Y");
  htemp->SetLabelSize(0.05,"X");
  htemp->SetLabelSize(0.05,"Y");
  htemp->SetStats(false);
  htemp->SetLineColor(kRed);
  htemp->Draw();
  fitY11->Draw("Same");
  fitY11->SetLineColor(kRed);
  catLeg2->AddEntry(fitY11,"Y11: Fit Only","l");  
  CorTerm->Draw("Same");
  CorTerm->SetLineColor(kBlue);
  catLeg2->AddEntry(CorTerm,"Landau Correction","l");
  Y11->Draw("Same");
  Y11->SetLineColor(kBlack);
  catLeg2->AddEntry(Y11,"Corrected Y11","l");  
  catLeg2->SetBorderSize(0);
  catLeg2->SetFillStyle(0);
  catLeg2->SetTextSize(0.04);
  catLeg2->Draw();

  //std::cout<<Y11TimePDF(-0.07)<<std::endl;
  //std::cout<<Y11TimePDF(-0.07)<<std::endl;
  //std::cout<<TMath::Landau(2,0,20)<<std::endl;
    
  //hand203new->Print("all");
  //hand203old->Print("all");
  //hand205->Print("all");
  
  c0->SaveAs("203_205_new205_206.pdf");
  //c1->SaveAs("Y11Comparison.pdf");
  //c2->SaveAs("Y11components.pdf");
  //h203_MC->Print("all");
  //h203_pulses->Print("all");
  //h205_pulses->Print("all");
  //Y11_MC->Print("all");
  //analyticSiPMHE_MC->Print("all");
}
