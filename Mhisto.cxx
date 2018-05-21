#include <vector>
#include <iostream>

#include <TH2D.h>

#include "Mhisto.h"

using std::vector;

Mhisto::Mhisto(int nbinsx_, double xmin_, double xmax_, int nbinsy_, double ymin_, double ymax_)
{
  nbinsx = nbinsx_;
  xmin = xmin_;
  xmax = xmax_;
  nbinsy = nbinsy_;
  ymin = ymin_;
  ymax = ymax_;
  
  xwidth = (xmax-xmin)/nbinsx;
  ywidth = (ymax-ymin)/nbinsy;

  vector<double> tempY(nbinsy,0);
  vector<vector<double>> tempXY(nbinsx,tempY);
  value = tempXY;
}

Mhisto::~Mhisto(){}

int Mhisto::Fill(double x,double y)
{
  int i = int((x-xmin)/xwidth);
  int j = int((y-ymin)/ywidth);
  if (i>=nbinsx||j>=nbinsy||i<0||j<0) {
    std::cout<<"Dimension exceeding in Mhisto::Fill()."<<std::endl;
    return -1;
  }
  value[i][j] += 1;
  return 1;
}  

double Mhisto::GetBinContent(int i, int j)
{
  if (i>=nbinsx||j>=nbinsy||i<0||j<0) {
    std::cout<<"Dimension exceeding in Mhisto::GetBinContent."<<std::endl;
    return -1;
  }
  return value[i][j];
}


double Mhisto::GetPointContent(double x,double y)
{
  int i = int((x-xmin)/xwidth);
  int j = int((y-ymin)/ywidth);
  if (i>=nbinsx||j>=nbinsy||i<0||j<0) {
    std::cout<<"Dimension exceeding in Mhisto::GetPointContent()."<<std::endl;
    return -1;
  }
  return GetBinContent(i,j);
}  

Mhisto::Mhisto(TH2D * hist)
{
  std::cout<<"This function is not available yet :)"<<std::endl;
}

TH2D * Mhisto::CopytoTH2(const char * name, const char * title)
{
  TH2D * histCopy = new TH2D(name,title,nbinsx,xmin,xmax,nbinsy,ymin,ymax);
  for (int i=0;i<nbinsx;++i) {
    for (int j=0;j<nbinsy;++j) {
      histCopy->SetBinContent(i+1,j+1,value[i][j]);
    }
  }
  std::cout<<"copy to TH2 finished"<<std::endl;
  return histCopy;
}

void Mhisto::Scale(double factor)
{
  for (int i=0;i<nbinsx;++i) 
    for (int j=0;j<nbinsy;++j)
      {
	//	std::cout<<value[i][j]<<std::endl;
	value[i][j] *= factor;
	//std::cout<<value[i][j]<<std::endl;
      }
}
