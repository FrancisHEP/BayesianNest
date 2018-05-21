// ***************************************************************
// This file was created using the bat-project script
// for project NestB.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
// #include <BAT/BCModelOutput.h> can not found this header file
#include <vector>
#include <fstream>
#include <stdio.h>
#include "NestBModel.h"

int main()
//int runNestB()
{
  BCAux::SetStyle();
  BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);
	 
  NestBModel m("NestBModel");
  m.SetData();
  m.SetSampleNum(300000);
  m.SetHist(42,3,45,100,1,3);
  
  m.SetPrecision(BCEngineMCMC::kMedium);
  m.SetNChains(4); // switch 5 chains to 1 chain

  // // expect to write to root file
  // BCModelOutput(m,"output.root");
  // BCModelOutput.WriteMarginalizedDistributions();
  // BCModelOutput.WriteMarkovChain(true);

  BCLog::OutSummary("Test model created");

  //  m.Normalize();
  m.MarginalizeAll(BCIntegrate::kMargMetropolis);
  m.FindMode(m.GetBestFitParameters());

  m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");
  m.PrintParameterPlot(m.GetSafeName() + "_parameters.pdf");
  m.PrintCorrelationPlot(m.GetSafeName() + "_correlation.pdf");
  m.PrintCorrelationMatrix(m.GetSafeName() + "_correlationMatrix.pdf");
  m.PrintKnowledgeUpdatePlots(m.GetSafeName() + "_update.pdf");
  m.PrintSummary();

  BCLog::OutSummary("Exiting");
  BCLog::CloseLog();

  return 0;
}
