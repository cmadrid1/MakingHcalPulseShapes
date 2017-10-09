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

int shiftMC  = 81; //shift MC Pulse
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
  TH1F* gtemp = new TH1F(name.c_str(),"",250,0,249);
  for(int i = 1; i < 250; i++){
    gtemp->SetBinContent(i,TMath::Landau(i,89+shiftpulses+phase,sigma));
  }
  gtemp->Scale(1.0/gtemp->Integral(0,250));
  return gtemp;
}

double Y11TimePDF(double t) {
  //double A, n, a, t0;
  //t/=0.4;
  //A=1/9868.9; n=2.528; a=-0.1518; t0=-0.0635;//original
  //return A*exp(a*t + t0)*pow(t, n);

  //return (1/8581.98)*exp(-t/3.1568003892614205)*pow(t,4.049441808115009);

  //Fit From Deconvolved Data
  double A,n,t0;
  A=0.104204; n=0.44064; t0=10.0186;
  return A*(1-exp(-t/n))*exp(-t/t0);
    
  //Fit from New data
  //return (1/20.3334)*5.51092*TMath::Landau(t,11.8357,3.75244);

  //Best One yet
  //t/=0.4;
  //return TMath::Landau(t,5,4.786);

  //
  //double wd1,ts1,wd2,ts2,wd3,ts3;
  //wd1=2.0, ts1=8.0, wd2=0.7, ts2=10.0, wd3=1.0, ts3=22.3; //125
  //return wd1*exp(-t/ts1) + wd2*exp(-t/ts2) + wd3*exp(-t/ts3);
}

double Y11Old(double t) {
  double A, n, a, t0;
  A=1/2485.9; n=2.528; a=-0.1518; t0=-0.0635;//original
  return A*exp(a*t + t0)*pow(t, n);
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

TH1F* hand_int(TH1F* hraw, string name){
  TH1F* htemp = new TH1F(name.c_str(),"",10,-0.5,9.5);
  for(unsigned int j = 1; j<=11; j++){
    htemp->SetBinContent(j,hraw->Integral(j*25,(j+1)*25));
  }
  //htemp->Scale(1.0/htemp->Integral(-0.5,9.5));
  return htemp;
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

  TH1F* htemp = new TH1F(name.c_str(),"",nBinsSiPM_,0.,nBinsSiPM_);
  for(unsigned int i = 1; i<nBinsSiPM_; i++){
    htemp->SetBinContent(i+shiftMC+phase,siPMShape2017_[i]);
  }
  htemp->Scale(nor/htemp->Integral(0,htemp->GetNbinsX()+1));
  return htemp;
}

TH1F* makeY11Histo(string name){
  TH1F* htemp = new TH1F(name.c_str(),"",nBinsSiPM_,0.,nBinsSiPM_);
  for(unsigned int i = 1; i<nBinsSiPM_; i++){
    htemp->SetBinContent(i+shift,Y11TimePDF(i));
  }
  htemp->Scale(1.0/htemp->Integral(0,htemp->GetNbinsX()+1));
  return htemp;
}

TH1F* makeOldY11Histo(string name){
  TH1F* htemp = new TH1F(name.c_str(),"",nBinsSiPM_,0.,nBinsSiPM_);
  for(unsigned int i = 1; i<nBinsSiPM_; i++){
    htemp->SetBinContent(i+shift,Y11Old(i));
  }
  htemp->Scale(1.0/htemp->Integral(0,htemp->GetNbinsX()+1));
  return htemp;
}

TH1F* makeAnalyticSiPMHE(string name){
  TH1F* htemp = new TH1F(name.c_str(),"",2*nBinsSiPM_-1,0,nBinsSiPM_);
  for(unsigned int i = 1; i<2*nBinsSiPM_-1; i++){
    htemp->SetBinContent(i+2*shift,analyticPulseShapeSiPMHE(i));
  }
  htemp->Scale(1.0/htemp->Integral(0,htemp->GetNbinsX()+1));
  return htemp;
}

TH1F* makePulseHisto(vector<double>& pulse, int bstart, string name){
	//make pulse template histo (from convolution)
	TH1F* htempT = new TH1F(name.c_str(),"",250,0,250);
	for(unsigned b = 0; b < pulse.size(); ++b){
		htempT->SetBinContent(b+bstart,pulse[b]);
	}
	htempT->Scale(nor/htempT->Integral(0,htempT->GetNbinsX()+1));
	return htempT;
}

void FixedSiPM(){
  TCanvas *c0 = new TCanvas("PulseShape_203","",1000,800);  
  TLegend* catLeg0 = new TLegend(0.73,0.65,0.99,0.9);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.16);
  //gPad->SetLogy();

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
  catLeg0->SetBorderSize(0);
  catLeg0->SetFillStyle(0);
  catLeg0->SetTextSize(0.04);
  catLeg0->AddEntry(h203_MC,"New MC","l");
  catLeg0->Draw();

  //TGraph* DeConvData = NULL;
  //DeConvData = makeGraph(17,deConvData);
  //DeConvData->SetMarkerColor(kGreen+2);
  //DeConvData->SetLineColor(kGreen+2);
  //DeConvData->Draw("P same");
  //catLeg0->AddEntry(DeConvData,"Deconvolved Data","P");
  //
  //TF1 *fa2 = new TF1("fa2","[0]*(1-exp(-(x-9)/0.44064))*exp(-(x-9)/10.0186) + [1]*TMath::Landau(x,[2]-9,[3])",0,250);
  //fa2->SetParLimits(0,0,15);
  //fa2->SetParLimits(1,1,15);
  //fa2->SetParLimits(2,0,15);
  //fa2->SetParLimits(3,5,15);
  //fa2->SetParameter(0,5.84053);
  //fa2->SetParameter(1,10.3613);
  //fa2->SetParameter(2,2);
  //fa2->SetParameter(3,2);
  ////NewData->Fit(fa2,"N","",0,50);
  //DeConvData->Fit(fa2,"N","",0,60);
  ////DeConvDataSpline->Fit(fa2,"N","",0,60);
  //fa2->SetLineColor(kBlue);
  //catLeg0->AddEntry(fa2,"Deconvolved Fit: Exp + Landau","l");
  //fa2->Draw("same");
  
  TH1F* Y11_MC = NULL;
  Y11_MC = makeY11Histo("Y11_MC");
  Y11_MC->SetLineColor(4);
  //Y11_MC->Draw("same hist");
  //catLeg0->AddEntry(Y11_MC,"Y11 (New FSU)","l");

  TH1F* Y11_MCOld = NULL;
  Y11_MCOld = makeOldY11Histo("Y11_MCOld");
  Y11_MCOld->SetLineColor(kGreen);
  //Y11_MCOld->Draw("same hist");
  //catLeg0->AddEntry(Y11_MCOld,"Y11","l");

  TH1F* analyticSiPMHE_MC = NULL;
  analyticSiPMHE_MC = makeAnalyticSiPMHE("analyticSiPMHE_MC");
  analyticSiPMHE_MC->SetLineColor(6);
  //analyticSiPMHE_MC->Draw("same hist");
  //catLeg0->AddEntry(analyticSiPMHE_MC,"analytic SiPM","l");

  TH1F* h205_pulses = NULL;
  h205_pulses = makePulseHisto(s205,69+shiftpulses+phase,"205");
  h205_pulses->SetLineColor(kBlue);
  //h205_pulses->Draw("same hist");
  //catLeg0->AddEntry(h205_pulses,"205","l");

  TH1F* new205 = NULL;
  new205 = makeLandau(4.57,"new 205");
  new205->Draw("hist same");
  new205->SetLineColor(kBlue);
  catLeg0->AddEntry(new205,"Data");

  TH1F* lowCharge = NULL;
  lowCharge = makeLandau(5.122,"lowCharge");
  lowCharge->SetLineColor(kOrange-2);
  //lowCharge->Draw("hist same");
  //catLeg0->AddEntry(lowCharge,"Testbeam 29,000fC");

  TH1F* highCharge = NULL;
  highCharge = makeLandau(4.324,"highCharge");
  highCharge->SetLineColor(kGreen+2);
  //highCharge->Draw("hist same");
  //catLeg0->AddEntry(highCharge,"Testbeam 336,000fC");

  TH1F* h203_pulses = NULL;
  h203_pulses = makePulseHisto(s203,65+shiftpulses+phase,"203");
  h203_pulses->SetLineColor(kBlack);
  h203_pulses->Draw("same hist");
  catLeg0->AddEntry(h203_pulses,"Current MC","l");
  
  //std::cout<<"203Root:  "<<h203_MC->GetMaximumBin()<<endl;
  //std::cout<<"203pulses.h:  "<<h203_pulses->GetMaximumBin()<<endl;
  //std::cout<<"205pulses.h:  "<<h205_pulses->GetMaximumBin()<<endl;
  //std::cout<<"Testbeam:    "<<highCharge->GetMaximumBin()<<endl;
  //std::cout<<"Y11 New:    "<<Y11_MC->GetMaximumBin()<<endl;

  /*
  TCanvas *c1 = new TCanvas("Hand_PulseShape_203","",1000,800);  
  TLegend* catLeg1 = new TLegend(0.64,0.65,0.99,0.9);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.16);

  TH1F* hand203new = NULL;
  hand203new = hand_int(h203_MC,"hand203new");
  hand203new->SetLineColor(1);  
  TH1F* hand203old = NULL;
  hand203old = hand_int(h203_pulses,"hand203old");
  hand203old->SetLineColor(2);
  TH1F* hand205 = NULL;
  hand205 = hand_int(h205_pulses,"hand205");
  hand205->SetLineColor(11);

  hand203new->Draw("hist"); 
  hand203new->SetMinimum(0.);
  hand203new->SetMaximum(0.75);
  hand203new->GetXaxis()->SetTitle("Time Step");
  hand203new->GetYaxis()->SetTitle("a.u.");
  hand203new->SetTitle("Hand Integrated");
  hand203new->SetName("Hand Integrated");
  hand203new->SetTitleSize(0.002);
  hand203new->SetTitleSize(0.07,"X");
  hand203new->SetTitleSize(0.08,"Y");
  hand203new->SetTitleOffset(1.0,"X");
  hand203new->SetTitleOffset(1.0,"Y");
  hand203new->SetLabelSize(0.05,"X");
  hand203new->SetLabelSize(0.05,"Y");
  hand203new->SetName("");
  hand203new->SetStats(false);
  
  catLeg1->SetTextSize(0.05);
  catLeg1->AddEntry(hand203new, "203 Corrected" ,"l");  
  catLeg1->Draw();

  hand203old->Draw("same hist");
  catLeg1->AddEntry(hand203old,"203 Current", "l");

  hand205->Draw("same hist");
  catLeg1->AddEntry(hand205,"205", "l");
  */

  //hand203new->Print("all");
  //hand203old->Print("all");
  //hand205->Print("all");

  //c0->SaveAs("log203comparison.pdf");
  //c1->SaveAs("hand_integrated.pdf");
  //h203_MC->Print("all");
  //h203_pulses->Print("all");
  //h205_pulses->Print("all");
  //Y11_MC->Print("all");
  //analyticSiPMHE_MC->Print("all");
}

