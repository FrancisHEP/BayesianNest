#include "QuantaNest.h"

#include <TMath.h>
#include <iostream>
#include <vector>

QuantaNest::QuantaNest ()
{
  // set default values
  energy = 5;                   // keV
  field = 1.03;                 // V/cm F_0
  type = 0;                     // 0: NR; 1: ER.
  nElectron = -1;
  nPhoton = -1;
  lightPE = -1;
  chargePE = -1;
  lindhardFactor = -1;
  epsilon = 11.5 * energy * TMath::Power(54, -7./3.);
  gamma = -1;
  kappa = -1;
  density = 2.9;
  resolution = (0.12724-0.032152*density-0.0013492*TMath::Power(density,2.))*1.5;
  zeta = gamma * TMath::Power(field, -0.062);
  nQuanta = -1;
  ratio = -1;
  lightQuenching = 1./(1 + 3.3 * TMath::Power(epsilon, 1.14));
  dpeFraction = 0.215;
  pmtResolution = 0.6;
  pde = 0.1;
  eee = 0.5;
  gasGain = 20;
  recombFluctuation = false;
  enableElectronLifeTime = false;
  calculateF0F1F2();
  tr.SetSeed(0);
}

void QuantaNest::setPars(const std::vector<double> & pars)
{
  alpha = pars[0];
  gamma = pars[1];
  zeta = gamma * pow(field, -0.062);
  lightQuenching = 1./(1 + 3.3 * pow(epsilon,1.14));
}

void QuantaNest::calculate()
{
  if (type==0) {
    getLindhardFactor();
  }
  calculateNQuanta();
  calculateNexNi();
  calculateRecombinationRate();
  calculatePhotonElectron();
  calculateCharge();
}

int QuantaNest::getNElectron()
{
  if (nElectron<0) {
    calculate();
  }
  return nElectron;
}

int QuantaNest::getNPhoton()
{
  if (nPhoton<0) {
    calculate();
  }
  return nPhoton;
}

double QuantaNest::getLightPE()
{
  // if (lightPE<0) {
  //   calculate();
  // }
  return lightPE;
}

double QuantaNest::getChargePE()
{
  // if (chargePE<0) {
  //   calculate();
  // }
  return chargePE;
}

void QuantaNest::setEnergy (double e)
{
  energy = e;
  epsilon = 11.5 * energy * TMath::Power(54, -7./3.);
  lightQuenching = 1./(1 + 3.3 * TMath::Power(epsilon, 1.14));
}

void QuantaNest::setField (double f)
{
  field = f;
  zeta = gamma * TMath::Power(field, -0.062);
}

void QuantaNest::setType (int t)
{
  type = t;
}

void QuantaNest::setDensity (double d)
{
  density = d;
  resolution = (0.12724-0.032152*density-0.0013492*TMath::Power(density,2.))*1.5;
}

double QuantaNest::getLindhardFactor ()
{
  gamma = 3 * TMath::Power(epsilon, 0.15) + 0.7 * TMath::Power(epsilon, 0.6) + epsilon;
  kappa = 0.1394;
  lindhardFactor = kappa * gamma / (1 + kappa * gamma);
  return lindhardFactor;
}

int QuantaNest::getNQuanta ()
{
  return nQuanta;
}

int QuantaNest::getNIonization ()
{
  return nIonization;
}

int QuantaNest::getNExcitation ()
{
  return nExcitation;
}

double QuantaNest::getRatio ()
{
  return ratio;
}

double QuantaNest::getRecombinationRate ()
{
  return recombinationRate;
}

double QuantaNest::getRecombinationRateT ()
{
  return recombinationRate_t;
}

void QuantaNest::calculateNQuanta ()
{
  double meanQuanta = energy*1000/W;
  double sigma = TMath::Sqrt(resolution*meanQuanta);
  nQuanta = (int) (tr.Gaus(meanQuanta, sigma));
  if (type == 0) {
    double smearedLF = tr.Gaus(lindhardFactor, 0.25*lindhardFactor);
    nQuanta = tr.Binomial(nQuanta, smearedLF);
  }
  if (nQuanta < 0) {
    nQuanta = 0;
  }
}

void QuantaNest::calculateNexNi ()
{
  if (nQuanta<=0) {
    nIonization = 0;
    nExcitation = 0;
    return;
  }
  calculateNexNiRatio();
  nExcitation = tr.Binomial(nQuanta, ratio/(1+ratio));
  nIonization = nQuanta - nExcitation;
}

void QuantaNest::calculateNexNiRatio ()
{
  // only for nr.
  if (type==0) {
    //ratio = 1.24 * TMath::Power(field, -0.0472)*(1 - TMath::Exp(-239*epsilon)) * 1.26;//for pandax-II
    //ratio = 1.24 * TMath::Power(field, -0.0472)*(1 - TMath::Exp(-239*epsilon)) * 1.19032258064516139;
    ratio = alpha * TMath::Power(field, -0.0472)*(1 - TMath::Exp(-239*epsilon)) * 1.286;
   // ratio = 1.24 * TMath::Power(field, -0.0472)*(1 - TMath::Exp(-239*epsilon)) * 1.0; // orignal parameter
    //ratio = 1.24 * TMath::Power(field, -0.0472)*(1 - TMath::Exp(-239*epsilon)) * 1.5; // for pandax-II engineering run
  } else if (type==1) {
    ratio = 0.059813 + 0.031228 * density;
  }
}

void QuantaNest::calculateRecombinationRate ()
{
  if (nIonization == 0) {
    recombinationRate = 0;
  } else {
    if (type == 0) {
      recombinationRate = 1 - TMath::Log(1+nIonization*zeta)/(nIonization*zeta);
    } else {
      calculateErRecombinationRate();
    }
    recombinationRate_t = recombinationRate;
    if (recombFluctuation) {
      while (1) {
        double fac = 0.1;
        if (type == 1) {
          fac = 0.127;
        }
        recombinationRate = tr.Gaus(recombinationRate_t, fac * TMath::Sqrt(recombinationRate_t));
        if (recombinationRate<1) {
//          std::cout << recombinationRate << std::endl;
          break;
        }
      }
    }
    if (recombinationRate<0) {
      recombinationRate = 0;
    }
    //    std::cout << recombinationRate << std::endl;
  }
}

void QuantaNest::calculateErRecombinationRate ()
{
  double tibF = 0.6347 * TMath::Exp(-0.00014 * field);
  double tibE = -0.373 * TMath::Exp(-(field*0.001)/tibF) + 1.5;
  double tibCurlA = 10 * TMath::Power(field, -0.04) * TMath::Exp(18/field);
  double tibCurlZ = 1 - TMath::Power(field, 0.2147) + 3;
  double tomasImel = tibF * TMath::Power(energy, -tibE) * (1 - TMath::Exp(-TMath::Power((energy - tibCurlZ)/tibCurlA, 0.188 * TMath::Power(field, 1./3)))) * TMath::Power(density/2.888, 0.3);
  recombinationRate = 1 - TMath::Log(1+tomasImel/4.*nIonization)/(tomasImel/4.*nIonization);
}

void QuantaNest::calculatePhotonElectron ()
{
  nPhoton = nExcitation + tr.Binomial(nIonization, recombinationRate);
  nElectron = nQuanta - nPhoton;
  if (type==0) {
    nPhoton = tr.Binomial(nPhoton, lightQuenching);
  }
}

void QuantaNest::calculateCharge ()
{
  // calculate s1
  // 1. detected photon
  int iPhoton = tr.Binomial(nPhoton, 1-f0);
  if (iPhoton>0) {
    // 2. generated pe
    int iPhoton2 = tr.Binomial(iPhoton, dpeFraction);
    iPhoton += iPhoton2;
    // 3. smearing to non-negative value
    lightPE = -1;
    while (lightPE <= 0) {
      lightPE = tr.Gaus(iPhoton, pmtResolution * TMath::Sqrt(iPhoton));
    }
  } else {
    lightPE = 0;
  }

  // calculate s2
  // 1. extracted electrons
  int iElectron;
  if (enableElectronLifeTime) {
    double surviveRate = TMath::Exp(-driftTime/electronLifeTime);
    int nsElectron = tr.Binomial(nElectron, surviveRate);
    iElectron = tr.Binomial(nsElectron, eee);
  } else {
    iElectron = tr.Binomial(nElectron, eee);
  }
  if (iElectron > 0) {
    double gVar = -1;
    while (gVar <= 0) {
      gVar = tr.Gaus(1, 0.165);
    }
    chargePE = -1;
    while (chargePE <=0 ) {
      chargePE = tr.Gaus(iElectron * gasGain * gVar, pmtResolution * TMath::Sqrt(iElectron * gasGain * gVar));
      // iElectron = iElectron*gVar*gasGain;
      // chargePE = tr.Gaus(iElectron,pmtResolution*TMath::Sqrt(iElectron));
    }
  } else {
    chargePE = 0;
  }
}

void QuantaNest::setPmtResolution (double r)
{
  pmtResolution = r;
}

void QuantaNest::setDpeFraction (double f)
{
  dpeFraction = f;
  calculateF0F1F2();
}

void QuantaNest::setPDE (double p)
{
  pde = p;
  calculateF0F1F2();
}

void QuantaNest::setEEE (double e)
{
  eee = e;
}

void QuantaNest::setGasGain (double g)
{
  gasGain = g;
}

void QuantaNest::setRecombFluctuation (bool t)
{
  recombFluctuation = t;
}

void QuantaNest::setElectronLifeTime (double t)
{
  if (t>0) {
    electronLifeTime = t;
    enableElectronLifeTime = true;
  }
}

void QuantaNest::setDriftTime (double t)
{
  if (t>0) {
    driftTime = t;
  }
}

void QuantaNest::calculateF0F1F2 ()
{
  f0 = 1 - pde/(1 + dpeFraction);
}
