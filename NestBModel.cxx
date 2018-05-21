// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************
#include "NestBModel.h"
#include "QuantaNest.h"
// #include "runGenMC.h"

#include <vector>
#include <cmath>
#include <memory>
#include <thread>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <fstream>

#include <TF1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TUUID.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TProfile.h>

#include <BAT/BCParameter.h>
#include <BAT/BCMath.h>
#include "TMath.h"
#include "Mhisto.h"

using namespace std;

NestBModel::NestBModel(const std::string& name)
    : BCModel(name)
{
  // alpha [0.5, 3.2] -> [1.8, 3]
  // gamma [0.005, 0.035] -> [0.006, 0.014]
  BCParameter * alpha = new BCParameter("alpha", 1.8,   3,     "alpha");
  BCParameter * gamma = new BCParameter("gamma", 0.006, 0.014, "gamma");
  AddParameter(*alpha);
  AddParameter(*gamma);
  SetPriorConstantAll();
  //  GetParameter("gamma").Fix(0.01385);
}

NestBModel::~NestBModel() {}

void NestBModel::SetData()
{
  // read position and energy info
  TFile filePositionX("energy_position/positionX.root","READ");
  TFile filePositionY("energy_position/positionY.root","READ");
  TFile filePositionZ("energy_position/positionZ.root","READ");
  TFile fileEnergySpec("energy_position/energySpecNr.root","READ");
  TFile fileData("energy_position/data_hist_rough.root","READ");

  positionX = (TH1F *)filePositionX.Get("h1");
  positionY = (TH1F *)filePositionY.Get("h1");
  positionZ = (TH1F *)filePositionZ.Get("h1");
  energySpecNr = (TH1F *)fileEnergySpec.Get("h1");
  data = (TH2D *)fileData.Get("hdata2D");
  
  positionX->SetDirectory(0);
  positionY->SetDirectory(0);
  positionZ->SetDirectory(0);
  energySpecNr->SetDirectory(0);
  data->SetDirectory(0);

  // read correction map 
  TFile fmap("energy_position/pandaxii_correction_map_run9_sm2_LRF_run9.root","READ");
  hcorrection_cy = (TH2F *) fmap.Get("hcorrection_cy");
  hcorrection_cy->SetDirectory(0);
  correctionMapS2 = hcorrection_cy;
  fvCorrectionMeanS2 = getMapMeanValue2D(correctionMapS2);

  // initialize fermi-dirac functions
  initEfficiencyFunctions();

  ifstream file("energy_position/data_points.txt");
  double xtemp, ytemp;
  while (1) {
    if (file.eof()!=0) break;
    file >> xtemp >> ytemp;
    dataX.push_back(xtemp);
    dataY.push_back(ytemp);
  }
  file.close();
  
  SampleNum = 10000;
  
  hist_nbinsx = 42;
  hist_xmin = 3;
  hist_xmax = 45;
  hist_nbinsy = 50;
  hist_ymin = 0;
  hist_ymax = 5;
}

void NestBModel::initEfficiencyFunctions ()
{
  qS1eff = new TF1("qS1eff", "1-1/(TMath::Exp((x-[0])/[1])+1)", 0, 50);
  qS2eff = new TF1("qS2eff", "1-1/(TMath::Exp((x-[0])/[1])+1)", 0, 10000);
  qS1Bdteff = new TF1("qS1Bdteff", "0.9716*(1-1/(TMath::Exp((x-[0])/[1])+1))", 0, 50);
  qS1eff->SetParameters(3.2, 0.44);
  qS2eff->SetParameters(72, 19);
  qS1Bdteff->SetParameters(2.567,0.9268);
}

//const int samplingNumber = 100000;

void NestBModel::SetSampleNum(int n)
{
  SampleNum = n;
}

void NestBModel::SetHist(int nbinsx,double xmin,double xmax,int nbinsy,double ymin,double ymax)
{
  hist_nbinsx = nbinsx;
  hist_xmin = xmin;
  hist_xmax = xmax;
  hist_nbinsy = nbinsy;
  hist_ymin = ymin;
  hist_ymax = ymax;
}

double NestBModel::LogLikelihood(const std::vector<double>& pars)
{
  // 42,3,45,500,0,5
  Mhisto histMC(hist_nbinsx, hist_xmin, hist_xmax, hist_nbinsy, hist_ymin, hist_ymax);
  GenMC(pars, histMC);
  double L = 0;

  // unbinned likelihood function
  // scale factor, including acceptance and normalize, is wrote in the GenMC()
  for (size_t i=0;i!=dataX.size();++i) {
    double prob = histMC.GetPointContent(dataX[i],dataY[i]);
    if (prob<=0) prob = 1e-8; // it was 1e-31 before. it seems too strict. 
    L += log(prob);
  }

  cout<<"likelihood L = "<<L<<endl;
  return L;
}

// derive TGraph(data)
TGraph * NestBModel::getDataTGraph()
{
  auto dataGraph = new TGraph(dataX.size(), &dataX[0], &dataY[0]);
  return dataGraph;
}

// derive TH2D(data)
TH2D * NestBModel::getDataTH2D()
{
  auto dataTH2D = new TH2D("dataTH2D","the histogram of data points",
			  hist_nbinsx, hist_xmin, hist_xmax, hist_nbinsy, hist_ymin, hist_ymax);
  for (size_t i=0;i<dataX.size();++i) {
    dataTH2D->Fill(dataX[i],dataY[i]);
  }
  return dataTH2D;
}

// get the Histogram(alpha,gamma), and return corresponding likelihood via reference..
TH2D * NestBModel::getHistMC(const std::vector<double> & pars, double & L)
{
  Mhisto histMC(hist_nbinsx, hist_xmin, hist_xmax, hist_nbinsy, hist_ymin, hist_ymax);
  GenMC(pars, histMC);
  L = 0;
  
  for (size_t i=0;i!=dataX.size();++i) {
    double prob = histMC.GetPointContent(dataX[i],dataY[i]);
    if (prob<=0) prob =  1e-8; // it was 1e-31 before. it seems too strict. 
    L += log(prob);
  }
  
  TH2D * histCopy = histMC.CopytoTH2("histMC","histMC");
  
  return histCopy;
}

// 
void NestBModel::drawMCandData(const std::vector<double> & pars, double & L)
{
  Mhisto histMC(hist_nbinsx, hist_xmin, hist_xmax, hist_nbinsy, hist_ymin, hist_ymax);
  GenMC(pars, histMC);
  L = 0;

  for (size_t i=0;i!=dataX.size();++i) {
    double prob = histMC.GetPointContent(dataX[i],dataY[i]);
    if (prob<=0) prob = 1e-8; // it was 1e-31 before. it seems too strict. 
    L += log(prob);
  }
  
  TH2D * histCopy = histMC.CopytoTH2("histMC","histMC");
  TH2D * histData = getDataTH2D();
  TGraph * graphData = getDataTGraph();
  
  TH1D * profileMC = histCopy->ProfileX();
  TH1D * profileData = histData->ProfileX();
  delete histData;
  
  TCanvas * histCanvas = new TCanvas("histCopy_canvas","",1);
  histCopy->GetYaxis()->SetRangeUser(0.8,3);
  histCopy->Draw("colz");
  graphData->SetMarkerColor(2);
  graphData->SetMarkerStyle(20);
  graphData->SetMarkerSize(0.3);
  graphData->SetMarkerColorAlpha(kRed, 0.6);
  graphData->Draw("p same");
  profileMC->SetLineColor(kAzure+10);
  profileMC->SetLineWidth(3);
  profileMC->Draw("hist same");
  profileData->SetLineColor(kRed-4);
  //  profileData->SetLineColorAlpha(kRed-4,0.8);
  profileData->SetLineWidth(3);
  profileData->Draw("hist same");

  histCanvas->Print("MCandData.pdf");

  TFile * file = new TFile("MCandData.root","recreate");
  histCopy->Write();
  graphData->Write();
  profileMC->Write();
  profileData->Write();
  file->Close();

  delete histCanvas;
  delete profileData;
  delete profileMC;
  delete graphData;
  delete histCopy;
  
}

// compare the reconstructed energy 
void NestBModel::drawReconstructedEnergy(TH2D * histCopy)
{
  // energyReconstruction
  TH1D * MCReconstructedEnergy = new TH1D("mc_energy","MC reconstructed energy histogram",50,0,10);
  for (int i=0;i<hist_nbinsx;++i) {
    for (int j=0;j<hist_nbinsy;++j) {
      double x = (hist_xmax-hist_xmin)*(double(i)+0.5)/hist_nbinsx + hist_xmin;
      double y = x*pow(10,(hist_ymax-hist_ymin)*(double(j)+0.5)/hist_nbinsy + hist_ymin);
      // 0.0137 keV pde 0.1114 eee 0.5450 seg 24.4
      double E = 0.0137*(x/0.1114+y/0.5450/24.4);
      MCReconstructedEnergy->Fill(E,histCopy->GetBinContent(i,j));
    }
  }

  TH1D * dataReconstructedEnergy = new TH1D("data_energy","data reconstructed energy histogram",50,0,10);
  for (int i=0;i<4331;++i) {
    double x = dataX[i];
    double y = x*pow(10,dataY[i]);
    double E = 0.0137*(x/0.1114+y/0.5450/24.4);
    dataReconstructedEnergy->Fill(E);
  }
  MCReconstructedEnergy->Scale(1/MCReconstructedEnergy->Integral());
  dataReconstructedEnergy->Scale(1/dataReconstructedEnergy->Integral());

  TCanvas * energy_canvas = new TCanvas("energy_canvas","",1);
  MCReconstructedEnergy->SetFillStyle(3144);
  MCReconstructedEnergy->GetYaxis()->SetRangeUser(0,0.12);
  MCReconstructedEnergy->SetLineColor(kAzure+1);
  MCReconstructedEnergy->SetFillColor(kAzure+1);
  MCReconstructedEnergy->Draw("hist");
  dataReconstructedEnergy->SetLineColor(kRed);
  dataReconstructedEnergy->SetLineWidth(3);
  dataReconstructedEnergy->Draw("b same");
  energy_canvas->Print("reconstructedEnergy.pdf");

  delete MCReconstructedEnergy;
  delete dataReconstructedEnergy;
  delete energy_canvas;
}

void NestBModel::draw1DSignalHist(TH2D * histCopy)
{
    // energyReconstruction
  TH1D * S1SignalHist = new TH1D("S1 signal","S1 signal spectrum from MC",47,3,50);
  TH1D * S2SignalHist = new TH1D("S2 signal","S2 signal spectrum from MC",50,0,3500);
  
  for (int i=0;i<hist_nbinsx;++i) {
    for (int j=0;j<hist_nbinsy;++j) {
      double x = (hist_xmax-hist_xmin)*(double(i)+0.5)/hist_nbinsx + hist_xmin;
      double y = x*pow(10,(hist_ymax-hist_ymin)*(double(j)+0.5)/hist_nbinsy + hist_ymin);
      S1SignalHist->Fill(x,histCopy->GetBinContent(i,j));
      S2SignalHist->Fill(y,histCopy->GetBinContent(i,j));
    }
  }

  S1SignalHist->Scale(1/S1SignalHist->Integral());
  S2SignalHist->Scale(1/S2SignalHist->Integral());

  TCanvas * canvas = new TCanvas("1D signal canvas","",1);
  S1SignalHist->SetFillColor(kAzure+6);
  S1SignalHist->GetYaxis()->SetRangeUser(0,0.12);
  S1SignalHist->Draw("hist");
  canvas->Print("S1Signal.pdf");
  S2SignalHist->SetFillColor(kAzure+1);
  S2SignalHist->GetYaxis()->SetRangeUser(0,0.12);
  S2SignalHist->Draw("hist");
  canvas->Print("S2Signal.pdf");

  delete S1SignalHist;
  delete S2SignalHist;
  delete canvas;

}

void NestBModel::GenMC(const std::vector<double> & pars, Mhisto & histMC)
{

  double xc, yc, zc;
  double energyNr;
  double qS1;
  double qS2, qS2raw;
  
  TRandom3 ttt;
  ttt.SetSeed(0);

  // prepare for variables
  double field = 317.5; // run9 400V/cm run10 317.5V/cm 
  double pde = 0.1134;  // run9 0.1114 run10 0.1134 
  double eee = 0.5759;  // run9 0.5450 run10 0.5759
  double seg = 23.9;    // run9 24.4   run10 23.9
  bool recombFluc = true;
  QuantaNest qn;
  qn.setPars(pars);
  qn.setField(field);
  qn.setPDE(pde);
  qn.setEEE(eee);
  qn.setGasGain (seg);
  qn.setRecombFluctuation(recombFluc);

  int i = 0;
  int accept_counts = 0; // events that finally fall into our observe window
  // it will undergo various fermi-dirac functions, cuts, and so on.
  while (i!=SampleNum) {
    ++i;
    
    energyNr = energySpecNr->GetRandom();
    xc = positionX->GetRandom();
    yc = positionY->GetRandom();
    zc = positionZ->GetRandom();
        
    if (zc<7245.7||zc>7549.7||xc*xc+yc*yc>7.2e4) continue;

    // previous mistakes
    // electron life time 758.35us -> 739.6us
    // electron drift velocity 1.714 -> 1.69 cm/s
    // double driftTime = (7583.5 - zc)/1.714;
    double driftTime = (7396 - zc)/1.69;
    qn.setDriftTime(driftTime);
    double electronLifeTime = 739.6;
    qn.setElectronLifeTime(electronLifeTime);
    
    qS1 = 0;
    qS2 = 0;

    if (energyNr>0) {
      qn.setEnergy(energyNr);
      qn.setType(0);
      qn.calculate();
      qS1 += qn.getLightPE();
      qS2 += qn.getChargePE();
    }
    if (qS1<=0 || qS2<=0) continue;
    if (qS1<3||qS1>45||qS2>10000) continue;
    
    // S2 position correction
    // find bin may change the value of correction map?
    double s2CorrectionZFactor = 1.0;
    double s2CorrectionTimeFactor = TMath::Exp(-(7583.5-zc)/1.69/electronLifeTime);
    double s2CorrectionRadiusFactor = correctionMapS2->GetBinContent(correctionMapS2->FindBin(xc,yc))/fvCorrectionMeanS2;
    qS2raw = qS2*s2CorrectionRadiusFactor;

    qS2 = qS2/s2CorrectionTimeFactor*s2CorrectionZFactor;
    
    if (qS2raw<100) continue;
    
    // S1 Fermi-Dirac Function cut
    double s1e = qS1eff->Eval(qS1);
    if (ttt.Uniform()>s1e) continue;
    // S1 boosted decision tree cut
    double bdte = qS1Bdteff->Eval(qS1);
    if (ttt.Uniform()>bdte) continue;
    // S2 Fermi-Dirac Function cut
    double s2e = qS2eff->Eval(qS2raw);
    if (ttt.Uniform()>s2e) continue;

    if (qS1<3||qS1>45||qS2raw<100||qS2>10000) continue;
    if (log10(qS2/qS1)<(1.074-0.595*exp(-qS1/6.38))) continue; 

    int succeed = histMC.Fill(qS1,log10(qS2/qS1));
    if (succeed==1) ++accept_counts;
    
    // ++i;
    cout<<100*i/SampleNum<<"\r";
    // if(i%(N/20)==0) cout <<100*i/N<<"%" << endl;
  }

  // cout<<"accept_counts: "<<accept_counts<<endl;
  cout<<"finished... HID: "<<int(qS1*1000)<<endl;

  double binArea = (hist_xmax - hist_xmin)*(hist_ymax - hist_ymin)/hist_nbinsx/hist_nbinsy;
  double acceptance = double(accept_counts)/double(SampleNum);

  cout<<"binarea = "<<binArea<<" acceptance = "<<acceptance<<endl;
  
  histMC.Scale(1/double(SampleNum)/binArea);
}

// the version that save ttree file.
void NestBModel::GenMC(const std::vector<double> & pars, Mhisto & histMC, const char * ofname)
{

  double xc, yc, zc;
  double energyNr;
  double qS1;
  double qS2, qS2raw;
  
  TRandom3 ttt;
  ttt.SetSeed(0);

  TFile fout(ofname, "RECREATE");
  TTree * tree = new TTree("histMC_tree","histMC");
  tree->Branch("energyNr",&energyNr,"energyNr/D");
  tree->Branch("qS1",&qS1,"qS1/D");
  tree->Branch("qS2",&qS2,"qS2/D");

  // prepare for variables
  double field = 317.5; // run9 400V/cm run10 317.5V/cm 
  double pde = 0.1134;  // run9 0.1114 run10 0.1134 
  double eee = 0.5759;  // run9 0.5450 run10 0.5759
  double seg = 23.9;    // run9 24.4   run10 23.9
  bool recombFluc = true;
  QuantaNest qn;
  qn.setPars(pars);
  qn.setField(field);
  qn.setPDE(pde);
  qn.setEEE(eee);
  qn.setGasGain (seg);
  qn.setRecombFluctuation(recombFluc);

  int i = 0;
  int accept_counts = 0; // events that finally fall into our observe window
  // it will undergo various fermi-dirac functions, cuts, and so on.
  while (i!=SampleNum) {
    ++i;
    
    energyNr = energySpecNr->GetRandom();
    xc = positionX->GetRandom();
    yc = positionY->GetRandom();
    zc = positionZ->GetRandom();
        
    if (zc<7245.7||zc>7549.7||xc*xc+yc*yc>7.2e4) continue;

    // previous mistakes
    // electron life time 758.35us -> 739.6us
    // electron drift velocity 1.714 -> 1.69 cm/s
    // double driftTime = (7583.5 - zc)/1.714;
    double driftTime = (7396 - zc)/1.69;
    qn.setDriftTime(driftTime);
    double electronLifeTime = 739.6;
    qn.setElectronLifeTime(electronLifeTime);
    
    qS1 = 0;
    qS2 = 0;

    if (energyNr>0) {
      qn.setEnergy(energyNr);
      qn.setType(0);
      qn.calculate();
      qS1 += qn.getLightPE();
      qS2 += qn.getChargePE();
    }
    if (qS1<=0 || qS2<=0) continue;
    if (qS1<3||qS1>45||qS2>10000) continue;
    
    // S2 position correction
    // find bin may change the value of correction map?
    double s2CorrectionZFactor = 1.0;
    double s2CorrectionTimeFactor = TMath::Exp(-(7583.5-zc)/1.69/electronLifeTime);
    double s2CorrectionRadiusFactor = correctionMapS2->GetBinContent(correctionMapS2->FindBin(xc,yc))/fvCorrectionMeanS2;
    qS2raw = qS2*s2CorrectionRadiusFactor;

    qS2 = qS2/s2CorrectionTimeFactor*s2CorrectionZFactor;
    
    if (qS2raw<100) continue;
    
    // S1 Fermi-Dirac Function cut
    double s1e = qS1eff->Eval(qS1);
    if (ttt.Uniform()>s1e) continue;
    // S1 boosted decision tree cut
    double bdte = qS1Bdteff->Eval(qS1);
    if (ttt.Uniform()>bdte) continue;
    // S2 Fermi-Dirac Function cut
    double s2e = qS2eff->Eval(qS2raw);
    if (ttt.Uniform()>s2e) continue;

    if (qS1<3||qS1>45||qS2raw<100||qS2>10000) continue;
    if (log10(qS2/qS1)<(1.074-0.595*exp(-qS1/6.38))) continue; 

    int succeed = histMC.Fill(qS1,log10(qS2/qS1));
    if (succeed==1) ++accept_counts;

    tree->Fill();

    // ++i;
    cout<<100*i/SampleNum<<"\r";
    // if(i%(N/20)==0) cout <<100*i/N<<"%" << endl;
  }

  // cout<<"accept_counts: "<<accept_counts<<endl;
  cout<<"finished... HID: "<<int(qS1*1000)<<endl;

  double binArea = (hist_xmax - hist_xmin)*(hist_ymax - hist_ymin)/hist_nbinsx/hist_nbinsy;
  double acceptance = double(accept_counts)/double(SampleNum);

  cout<<"binarea = "<<binArea<<" acceptance = "<<acceptance<<endl;
  
  histMC.Scale(1/double(SampleNum)/binArea);

  tree->Write();
  fout.Close();
  delete tree;
}

double NestBModel::getMapMeanValue2D(TH2F *map)
{
  Double_t sum = 0;
  Int_t nSamples = 0;
  Double_t x, y;
  for (Int_t iX=0; iX<map->GetNbinsX(); ++iX) {
    for (Int_t iY=0; iY<map->GetNbinsY(); ++iY) {
      x = map->GetXaxis()->GetBinCenter(iX);
      y = map->GetYaxis()->GetBinCenter(iY);
      if ((x*x+y*y)<72000) {
	sum += map->GetBinContent(iX, iY);
	++nSamples;
      }
    }
  }
  return sum/nSamples;
}
