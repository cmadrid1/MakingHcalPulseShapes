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

void Jay_Make_Method3_206(){
  std::vector<float> landauFrac = {
    0, 7.6377e-05, 0.000418655, 0.00153692, 0.00436844, 0.0102076, 
    0.0204177, 0.0360559, 0.057596, 0.0848493, 0.117069, 0.153152, 0.191858, 0.23198, 0.272461, 0.312438, 
    0.351262, 0.388476, 0.423788, 0.457036, 0.488159, 0.517167, 0.54412, 0.569112, 0.592254, 0.613668, 
    0.633402, 0.651391, 0.667242, 0.680131, 0.688868, 0.692188, 0.689122, 0.67928, 0.662924, 0.64087, 
    0.614282, 0.584457, 0.552651, 0.51997, 0.487317, 0.455378, 0.424647, 0.395445, 0.367963, 0.342288, 
    0.318433, 0.29636, 0.275994, 0.257243, 0.24, 0.224155, 0.2096, 0.196227, 0.183937, 0.172635, 
    0.162232, 0.15265, 0.143813, 0.135656, 0.128117, 0.12114, 0.114677, 0.108681, 0.103113, 0.0979354, 
    0.0931145, 0.0886206, 0.0844264, 0.0805074, 0.0768411, 0.0734075, 0.0701881, 0.0671664, 0.0643271, 
    0.0616564, 0.0591418, 0.0567718, 0.054536, 0.0524247, 0.0504292, 0.0485414, 0.046754, 0.0450602, 
    0.0434538, 0.041929, 0.0404806, 0.0391037, 0.0377937, 0.0365465, 0.0353583, 0.0342255, 0.0331447, 
    0.032113, 0.0311274, 0.0301854, 0.0292843, 0.0284221, 0.0275964, 0.0268053, 0.0253052, 0.0238536, 
    0.0224483, 0.0210872, 0.0197684, 0.0184899, 0.01725, 0.0160471, 0.0148795, 0.0137457, 0.0126445, 
    0.0115743, 0.0105341, 0.00952249, 0.00853844, 0.00758086, 0.00664871,0.00574103, 0.00485689, 0.00399541, 
    0.00315576, 0.00233713, 0.00153878, 0.000759962, 0 };

  TH1F* LandauFrac = new TH1F("landauFrac","",126,0,125);
  for(unsigned int i = 1; i<=125; i++){
    LandauFrac->SetBinContent(i,landauFrac[i]);
  }  
  std::cout<<"Landau:  "<<LandauFrac->Integral(0,250+1)<<std::endl;

  std::vector<float> siPM205Frac = {
    0, 0, 0, 0, 0.00133129, 0.00444633, 0.0115, 0.0243992, 0.0443875, 0.0716386, 0.105298, 0.143832,
    0.185449, 0.228439, 0.271367, 0.31315, 0.353041, 0.390587, 0.425555, 0.45788, 0.487604, 0.514843,
    0.539752, 0.562504, 0.583282, 0.602263, 0.619612, 0.635457, 0.649765, 0.66208, 0.671249, 0.675509,
    0.673048, 0.662709, 0.644394, 0.619024, 0.588194, 0.55375, 0.517448, 0.480768, 0.444831, 0.410418,
    0.378015, 0.347879, 0.320103, 0.294667, 0.271474, 0.250391, 0.231257, 0.213907, 0.198178, 0.183914,
    0.170967, 0.159205, 0.148505, 0.138758, 0.129864, 0.121737, 0.114299, 0.107478, 0.101214, 0.0954507,
    0.0901402, 0.0852385, 0.0807069, 0.0765108, 0.0726194, 0.0690052, 0.0656435, 0.0625123, 0.0595916, 0.0568637,
    0.0543125, 0.0519236, 0.0496838, 0.0475815, 0.0456058, 0.0437472, 0.0419966, 0.0403463, 0.0387887, 0.0373173,
    0.0359259, 0.034609, 0.0333615, 0.0321786, 0.0310561, 0.02999, 0.0289767, 0.0280127, 0.0270951, 0.0262209,
    0.0253875, 0.0245923, 0.0238333, 0.0231082, 0.022415, 0.021752, 0.0211174, 0.0205097, 0.0199274, 0.0193692,
    0.0188336, 0.0183196, 0.017826, 0.0173518, 0.0168959, 0.0164575, 0.0160356, 0.0156296, 0.0152385, 0.0148617,
    0.0144984, 0.0141482, 0.0138103, 0.0134842, 0.0131693, 0.0128652, 0.0125714, 0.0122873, 0.0120127, 0.011747,
    0.01149, 0.0112412, 0.0110002 };

  TH1F* SiPM205Frac = new TH1F("siPM205Frac","",126,0,125);
  for(unsigned int i = 1; i<=125; i++){
    SiPM205Frac->SetBinContent(i,siPM205Frac[i]);
  }
  std::cout<<"205:  "<<SiPM205Frac->Integral(0,250+1)<<std::endl;

  std::vector<float> siPM206Frac = {
    0,0,0, 7.99045e-05, 0.0033365, 0.0157494, 0.0387623, 0.0701868, 0.10687, 0.146053, 0.185725, 0.224641, 0.262309, 0.298594, 
    0.33347, 0.366936, 0.398996, 0.429654, 0.45892, 0.486807, 0.513336, 0.538531, 0.562426, 0.585055, 
    0.606461, 0.626688, 0.64578, 0.663788, 0.680759, 0.696664, 0.708454, 0.710198, 0.700502, 0.681598, 
    0.656686, 0.628567, 0.599294, 0.570153, 0.541672, 0.514023, 0.487267, 0.461437, 0.43656, 0.41266, 
    0.389755, 0.367858, 0.346973, 0.327095, 0.308216, 0.290318, 0.273379, 0.257373, 0.242269, 0.228033, 
    0.214633, 0.20203, 0.190189, 0.179072, 0.168642, 0.158863, 0.1497, 0.141116, 0.133079, 0.125557, 
    0.118518, 0.111932, 0.105772, 0.100009, 0.09462, 0.0895793, 0.0848645, 0.0804542, 0.0763281, 0.0724674, 
    0.0688541, 0.0654717, 0.0623044, 0.0593378, 0.0565582, 0.0539527, 0.0515095, 0.0492176, 0.0470665, 0.0450467, 
    0.0431491, 0.0413655, 0.039688, 0.0381094, 0.036623, 0.0352226, 0.0339022, 0.0326567, 0.0314808, 0.03037, 
    0.0293199, 0.0283265, 0.027386, 0.0264951, 0.0256504, 0.0248489, 0.0240879, 0.0233648, 0.0226771, 0.0220227, 
    0.0213994, 0.0208053, 0.0202386, 0.0196976, 0.0191809, 0.0186868, 0.0182141, 0.0177616, 0.0173279, 0.0169122, 
    0.0165132, 0.0161301, 0.0157621, 0.0154082, 0.0150676, 0.0147398, 0.0144239, 0.0141195, 0.0138258, 0.0135423, 
    0.0132685};//, 0.0130039, 0.0127482};

  TH1F* SiPM206Frac = new TH1F("siPM206Frac","",126,0,125);
  for(unsigned int i = 1; i<=125; i++){
    SiPM206Frac->SetBinContent(i,siPM206Frac[i]);
  }  
  std::cout<<"206:  "<<SiPM206Frac->Integral(0,250+1)<<std::endl;
  
  TCanvas *canvas = new TCanvas("canvas","canvas",1000,800);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.16);

  TLegend* catLeg1 = new TLegend(0.55,0.7,0.9,0.9);
  
  TH1F *method2_Shape = NULL;
  method2_Shape = compute206Shape("206");
  method2_Shape->SetLineColor(kRed);

  float thenewShapeSum=0;
  vector<float> thenewShape;
  float shapeMax=0;
  //i range used to align the start of the pulse with 0
  for (int i=-25; i<225; i++) {
    thenewShape.push_back(method2_Shape->Integral(i,i+25));
    thenewShapeSum+=method2_Shape->Integral(i,i+25);
    if (method2_Shape->Integral(i,i+25)>shapeMax) shapeMax=method2_Shape->Integral(i,i+25);
  }

  std::cout << thenewShapeSum << std::endl;

  TGraph *method3_Shape = new TGraph();

  float sum2=0;

  for (int i=0; i<thenewShape.size(); i++) {
    //method3_Shape->SetPoint(i,i+100,thenewShape.at(i)*method2_Shape->GetMaximum()/shapeMax);
    method3_Shape->SetPoint(i,i,thenewShape.at(i));
    sum2+=thenewShape.at(i)*method2_Shape->GetMaximum()/shapeMax;
  }
  
  method2_Shape->SetMinimum(0.000001);
  method2_Shape->SetMaximum(0.04);
  method2_Shape->GetXaxis()->SetRange(0,120);
  method2_Shape->GetXaxis()->SetTitle("time (ns)");
  method2_Shape->GetYaxis()->SetTitle("a.u.");
  method2_Shape->SetTitle("206 Shapes");
  method2_Shape->SetName("206 Shapes");
  method2_Shape->SetTitleSize(0.002);
  method2_Shape->SetTitleSize(0.07,"X");
  method2_Shape->SetTitleSize(0.08,"Y");
  method2_Shape->SetTitleOffset(1.0,"X");
  method2_Shape->SetTitleOffset(1.0,"Y");
  method2_Shape->SetLabelSize(0.05,"X");
  method2_Shape->SetLabelSize(0.05,"Y");
  method2_Shape->SetStats(false);
  method2_Shape->Draw("hist");
  method2_Shape->Scale(shapeMax/method2_Shape->GetMaximum());
  catLeg1->SetTextSize(0.03);
  catLeg1->AddEntry(method2_Shape,"M2: 206 Shape","l");
  catLeg1->Draw();

  method3_Shape->SetLineColor(kBlack);
  //method3_Shape->Draw("same el");
  LandauFrac->SetLineColor(kBlue);
  LandauFrac->Draw("same hist");
  SiPM205Frac->SetLineColor(kGreen+2);
  SiPM205Frac->Draw("same hist");
  SiPM206Frac->SetLineColor(kGray);
  SiPM206Frac->Draw("same hist");
  
  //catLeg1->AddEntry(method3_Shape,"M3: 206 Shape","l");
  catLeg1->AddEntry(LandauFrac,"M3: Landau Shape","l");
  catLeg1->AddEntry(SiPM205Frac,"M3: 205 Shape","l");
  catLeg1->AddEntry(SiPM206Frac,"M3: 206 Shape","l");
  //method3_Shape->Draw("apel");

  cout << endl << endl << "Method 3 Shape: "<< endl;
  for (int i=0; i<125; i++) {
    double tmp, val;
    method3_Shape->GetPoint(i,tmp,val);
    cout << val << ", ";
    if (i%10==1 && i>1) cout << endl;
  }
  cout << endl;
  canvas->SaveAs("Method3Shape_206.pdf");  
}
