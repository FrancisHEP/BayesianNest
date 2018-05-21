#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TFile.h>

#include "NestBModel.h"
#include "Mhisto.h"
using namespace std;

// plot the hist2d of monte carlo at a given pair of parameters
// ./manualScan 1.8 0.014 20000
void drawHistMC(double alpha, double gamma, int SampleNum)
{
  NestBModel m("NestBModel");
  m.SetData();
  m.SetSampleNum(SampleNum);
  m.SetHist(42,3,45,300,0,3);
  
  vector<double> pars = {alpha,gamma};
  double likelihood = 0;

  cout<<"initialized successfully..."<<endl;

  m.drawMCandData(pars,likelihood);

  cout<<"likelihood is "<<likelihood<<endl;
}

void drawLikelihoodHist(int stepX, double alpha_min, double alpha_max, int stepY, double gamma_min, double gamma_max, int SampleNum)
{
  // int stepX = 32;
  // int stepY = 32;
  // int alpha_min = 0.5, alpha_max = 3;
  // int gamma_min = 0.005, gamma_max = 0.035;
  double widthX = (alpha_max - alpha_min)/stepX;
  double widthY = (gamma_max - gamma_min)/stepY;
  TH2D * hist_likelihood = new TH2D("hist_likelihood","hist_likelihood",
				    stepX, alpha_min, alpha_max, stepY, gamma_min, gamma_max);

  NestBModel m("NestBModel");
  m.SetData();
  m.SetSampleNum(SampleNum);
  m.SetHist(42,3,45,100,1,3);

  for (int i=0;i<stepX;++i) {
    for (int j=0;j<stepY;++j) {
      double alpha = alpha_min + (double(i)+0.5)*widthX;
      double gamma = gamma_min + (double(j)+0.5)*widthY;
      vector<double> pars = {alpha, gamma};
      cout<<"alpha "<<pars[0]<<" gamma "<<pars[1]<<endl;
      double L = m.LogLikelihood(pars);
      hist_likelihood->SetBinContent(i+1,j+1,L);
    }
  }

  TCanvas * c = new TCanvas("hist_likelihood","hist_likelihood",1);
  hist_likelihood->Draw("colz");
  c->Print("hist_likelihood.pdf");

  TFile * histFile = new TFile("hist_likelihood.root","recreate");
  hist_likelihood->Write();
  histFile->Close();

  delete histFile;
  delete c;
}

void fluctuationTest()
{
  NestBModel m("NestBModel");
  m.SetData();
  m.SetSampleNum(300000);
  m.SetHist(42,3,45,100,0,3);
  vector<double> pars = {2.61, 0.00961};

  auto likelihoodHist = new TH1D("likelihoodHist","",10000,-30e3,-10e3);
  auto LMean2 = new TH1D("LMean2","",10000,-30e3,-10e3);
  auto LMean3 = new TH1D("LMean3","",10000,-30e3,-10e3);
  auto LMean5 = new TH1D("LMean5","",10000,-30e3,-10e3);
  auto LMean10 = new TH1D("LMean10","",10000,-30e3,-10e3);

  double sum2 = 0, sum3 = 0, sum5 = 0, sum10 = 0;
  for (size_t i=1;i<=10;++i) {
    double L = m.LogLikelihood(pars);
    cout<<L<<endl;
    likelihoodHist->Fill(L);
    sum2 += L; sum3 += L; sum5 += L; sum10 += L;
    if (i%2 == 0) {
     LMean2->Fill(sum2/2);
     sum2 = 0;
    }
    if (i%3 == 0) {
     LMean3->Fill(sum3/3);
     sum3 = 0;
    }
    if (i%5 == 0) {
     LMean5->Fill(sum5/5);
     sum5 = 0;
    }
    if (i%10 == 0) {
     LMean10->Fill(sum10/10);
     sum10 = 0;
    }
  }

  /*

mean/sigma value

1,0000 times
-0.0144496  1
-0.00881559 2
-0.00843399 3
-0.00479853 5
-0.00334083 10 

10,0000 times 
-0.00370492 1
-0.00295318 2
-0.00245378 3
-0.00188652 5
-0.00155471 10

100,0000 times
-0.000716168 1
-0.000598745 2
-0.000501434 3
-0.000393784 5


  */
  cout<<endl;
  cout<<likelihoodHist->GetMean()<<endl;
  cout<<likelihoodHist->GetStdDev()<<endl;
  cout<<likelihoodHist->GetStdDev()/likelihoodHist->GetMean()<<endl;
  cout<<LMean2->GetStdDev()/LMean2->GetMean()<<endl;
  cout<<LMean3->GetStdDev()/LMean3->GetMean()<<endl;
  cout<<LMean5->GetStdDev()/LMean5->GetMean()<<endl;
  cout<<LMean10->GetStdDev()/LMean10->GetMean()<<endl;
  
  TCanvas * c = new TCanvas("hist","",1);
  likelihoodHist->Draw();
  c->Print("likelihoodFlucTest.pdf");

  TFile * histFile = new TFile("LFlucTest.root","recreate");
  likelihoodHist->Write();
  histFile->Close();
  
  return;
}

int main(int argc, char * argv[])
{
  /*
    test the fluctuation of the likelihood function
    $ ./manualScan 1
  */
  if (atoi(argv[1]) == 1) {
    fluctuationTest();
    return 0;
  }

  
  if (argc==4) {
    drawHistMC(atof(argv[1]),atof(argv[2]),atoi(argv[3]));
    return 0;
  }

  /*
    draw likelihood in parameter space
    $ ./manualScan draw_likelihood 10000 
  */
  if (argc==3) {
    drawLikelihoodHist(32,0.5,3,32,0.005,0.035,atoi(argv[2]));
    return 0;
  }
  
  return 0;
}
