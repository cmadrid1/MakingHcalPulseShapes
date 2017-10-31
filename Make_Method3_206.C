#include "TMath.h"
#include "TGraph.h"

static const int nBinsSiPM_ = 251;
std::vector<float> siPMShape(nBinsSiPM_,0.0);
static const double frac_ = 0.13;

template <class F1, class F2>
static std::vector<double> convolve(unsigned nbin, F1 f1, F2 f2){
  std::vector<double> result(2*nbin-1,0.);
  for(unsigned i = 0; i < 2*nbin-1; ++i){
    for(unsigned j = 0; j < std::min(i+1,nbin); ++j){
      double tmp = f1(j)*f2(i-j);
      if(std::isnan(tmp) or std::isinf(tmp)) continue;
      result[i] += tmp;
    }
  }
  return result;
}

double fitY11results(double t){
  //Fit From Deconvolved Data
  double A,n,t0;
  A=0.104204; n=0.44064; t0=10.0186;
  return A*(1-exp(-t/n))*exp(-t/t0);
}

double corTerm(double t){
  double norm,mpv,sigma;
  norm=0.0806123; mpv=0; sigma=20;
  return norm*TMath::Landau(t,mpv,sigma);  
}

double Y11TimePDF(double t) {
  return (1-frac_)*fitY11results(t) + frac_*corTerm(t);
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

TH1F* compute206Shape(string name){
  //numerical convolution of SiPM pulse + WLS fiber shape
  std::vector<double> nt = convolve(nBinsSiPM_,analyticPulseShapeSiPMHE,Y11TimePDF);

  double norm = 0.;
  for (unsigned int j = 1; j <= nBinsSiPM_; ++j) {
    norm += (nt[j]>0) ? nt[j] : 0.;
  }
  for (unsigned int j = 1; j <= nBinsSiPM_; ++j) {
    nt[j] /= norm;
    siPMShape[j]=nt[j];
  }

  TH1F* htemp = new TH1F(name.c_str(),"",nBinsSiPM_+1,0.,nBinsSiPM_);
  for(unsigned int i = 1; i<=nBinsSiPM_+1; i++){
    htemp->SetBinContent(i,siPMShape[i]);
  }
  htemp->Scale(1.0/htemp->Integral(0,htemp->GetNbinsX()+1));
  return htemp;
}

void Make_Method3_206(){
  TH1F *method2_Shape = NULL;
  method2_Shape = compute206Shape("206");
  method2_Shape->SetLineColor(kMagenta);

  float thenewShapeSum=0;
  vector<float> thenewShape;
  float shapeMax=0;
  for (int i=0; i<250; i++) {
    thenewShape.push_back(method2_Shape->Integral(i,i+25));
    thenewShapeSum+=method2_Shape->Integral(i,i+25);
    if (method2_Shape->Integral(i,i+25)>shapeMax) shapeMax=method2_Shape->Integral(i,i+25);
  }

  TGraph *method3_Shape = new TGraph();

  float sum2=0;
  for (int i=0; i<4; i++) {
    method3_Shape->SetPoint(i,i,0);
  }

  for (int i=4; i<thenewShape.size(); i++) {
    method3_Shape->SetPoint(i,i+100,thenewShape.at(i-4)*method2_Shape->GetMaximum()/shapeMax);
    sum2+=thenewShape.at(i)*method2_Shape->GetMaximum()/shapeMax;
  }
  
  method2_Shape->GetYaxis()->SetRangeUser(0,0.04);
  method2_Shape->Draw("hist");
  method3_Shape->SetLineColor(kBlack);
  method3_Shape->Draw("same l");

  cout << endl << endl << "Method 3 Shape: "<< endl;
  for (int i=0; i<125; i++) {
    double tmp, val;
    method3_Shape->GetPoint(i,tmp,val);
    cout << val << ", ";
    if (i%10==1 && i>1) cout << endl;
  }
  cout << endl;
  
}
