// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************
// It's a combination of QuantaNest and BAT.
// Wenbo Ma
// ***************************************************************

#ifndef __BAT__NESTBMODEL__H
#define __BAT__NESTBMODEL__H

#include <BAT/BCModel.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TGraph.h>
#include <string>
#include <vector>

#include "Mhisto.h"

class NestBModel : public BCModel
{

public:

  NestBModel(const std::string& name);
  ~NestBModel();

  double LogLikelihood(const std::vector<double>& pars);
  void SetData();

private:
	 
  TH1F * positionX;
  TH1F * positionY;
  TH1F * positionZ;
  TH1F * energySpecNr;
	   
  TH2D * data;
  
  TH2F * hcorrection_cy;
  TH2F * correctionMapS2;
  double fvCorrectionMeanS2;

  TF1 *qS1eff;
  TF1 *qS2eff;
  TF1 *qS1Bdteff;

  // read real data points
  std::vector<double> dataX;
  std::vector<double> dataY;
  // ===
  
  int SampleNum;

  int hist_nbinsx, hist_nbinsy;
  double hist_xmin, hist_xmax, hist_ymin, hist_ymax;

public:
  
  void initEfficiencyFunctions();
  
  void SetSampleNum(int n);
  void SetHist(int nbinsx,double xmin,double xmax,int nbinsy,double ymin,double ymax);
  
  TGraph * getDataTGraph();
  TH2D * getDataTH2D();

  TH2D * getHistMC(const std::vector<double> & pars, double & L);

  void drawMCandData(const std::vector<double> & pars, double & L);
  void drawReconstructedEnergy(TH2D * histCopy);
  void draw1DSignalHist(TH2D * histCopy);

  void GenMC(const std::vector<double> & pars, Mhisto & histMC);
  void GenMC(const std::vector<double> & pars, Mhisto & histMC, const char * ofname);
  double getMapMeanValue2D(TH2F *map);


};

#endif
