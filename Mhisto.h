#ifndef MHISTO_H
#define MHISTO_H

#include <vector>
#include <TH2D.h>

class Mhisto
{

 public:

  Mhisto(int nbinsx_, double xmin_, double xmax_, int nbinsy_, double ymin_, double ymax_);
  Mhisto(TH2D * hist);
  ~Mhisto();
  
  int Fill(double x, double y);
  double GetBinContent(int i, int j);
  double GetPointContent(double x, double y);
  TH2D * CopytoTH2(const char * name, const char * title);
  void Scale(double factor);
  
 public:

  std::vector<std::vector<double>> value;

  double nbinsx;
  double nbinsy;

 private:
  
  double xmin;
  double xmax;
  double xwidth;

  double ymin;
  double ymax;
  double ywidth;
  
};

#endif
