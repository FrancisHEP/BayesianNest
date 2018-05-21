#ifndef QUANTANEST_HH
#define QUANTANEST_HH

#include <TRandom3.h>

class QuantaNest {

public:
  QuantaNest ();
  ~QuantaNest () {}

  void calculate();

  int getNElectron();
  int getNPhoton();

  double getLightPE();
  double getChargePE();

  void setEnergy(double e);

  void setField(double f);

  void setType (int t);

  void setDensity (double d);

  double getLindhardFactor ();

  int getNQuanta ();

  int getNIonization ();

  int getNExcitation ();

  double getRatio ();

  double getRecombinationRate();

  double getRecombinationRateT();

  void setPmtResolution (double r);

  void setDpeFraction (double f);

  void setPDE (double p);

  void setEEE (double e);

  void setGasGain (double g);

  void setRecombFluctuation (bool t);

  void setElectronLifeTime (double t);

  void setDriftTime (double t);
private:

public:
  void setPars(const std::vector<double> & pars);
    
  void calculateNQuanta();
  void calculateNexNi ();
  void calculateNexNiRatio ();
  void calculateRecombinationRate();
  void calculateErRecombinationRate();
  void calculatePhotonElectron();
  void calculateCharge();
  void calculateF0F1F2();
  // parameters for input/output
  double energy;
  double field;

  int type;

  int nElectron;
  int nPhoton;

  double lightPE;
  double chargePE;

  double density;

  // parameters used in the calculation
  double lindhardFactor;
  double epsilon;
  double kappa;

  double resolution;            // for quanta generation
  int nQuanta;                  // number of quanta generated
  double ratio;                 // ratio of intial Nex/Ni
  int nIonization;              // number of initial ionization
  int nExcitation;              // number of initial excitation

  double recombinationRate;     // rate of recombination
  double recombinationRate_t;	// rate of combination without fluctuation
  double lightQuenching;        // quenching factor for light

  double pmtResolution;		// resolution of PMT
  double dpeFraction;		// fraction of double photon emission
  double pde;			// photon detection efficiency
  double eee;			// electron extraction efficiency
  double f0;			// detection inefficiency of photon
  double gasGain;		// gain of gas

  bool recombFluctuation;       // control whether to use the fluctuation of recombination
  bool enableElectronLifeTime;  // enable the lifetime
  double electronLifeTime;      // value of electron life time.
  double driftTime;             // drift time

  TRandom3 tr;

  // NEST free parameters; 8+1+1
  double alpha,gamma;
  // gVar: gas-gain-fluctuation;
  // temporarily fixed in calculateCharge() 
  // gVar = 0.165; OR gVar = Gaus(1,0.165)
  double gVar;
  // zeta is considered not to be a real parameter.
  // zeta = function(field|gamma,delta) in my opinion.
  // zeta is declared and computed in constructor.
  double zeta;
  double W = 13.7; // work function, 13.7eV

};

#endif // QUANTANEST_HH
