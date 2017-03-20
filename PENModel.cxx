// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "PENModel.h"
#include <cmath>
#include <iostream>
#include <BAT/BCMath.h>

PENModel_combined::PENModel_combined(std::vector<PENdata> *adata, const char * name) : BCModel(name), data(adata){

}

double PENModel_combined::LogLikelihood(const std::vector<double> & parameters) {
  double LLH = 0;

  double tau = parameters[0];
  double protdeteff = parameters[1];
  unsigned int nentries = data->size();
  for (unsigned int i = 0; i < nentries; i++){
	PENdata *d = &((*data)[i]);
    double N0 = parameters[2+i*3];
    double noisen = parameters[3+i*3];
    double noisep = parameters[4+i*3];
    //LLH += BCMath::LogPoisson(d->monitorcount, N0 + noisen*d->monitortime);
	double decay = std::exp(-d->storagetime/tau); // fraction of stored neutrons that have survived after storage time
	LLH += BCMath::LogPoisson(d->neutroncount, N0*decay + noisen*d->neutrontime);
	for (unsigned int j = 0; j < d->protonbins; j++){
	  double t = d->pbinwidth*j;
	  double rate = std::exp(-t/tau) - std::exp(-(t + d->pbinwidth)/tau); // integrate decay rate over bin width
	  //      std::cout << pd.protoncounts[j] << "n vs. " << Np*decay << '\n';
	  LLH += BCMath::LogPoisson(d->protoncounts[j], protdeteff*N0*rate + noisep*d->pbinwidth);
	}
	LLH += BCMath::LogPoisson(d->neutronbackground, noisen*d->bgtime);
	LLH += BCMath::LogPoisson(d->protonbackground, noisep*d->bgtime);
  }
  return LLH;
}


PENModel_neutrons::PENModel_neutrons(std::vector<PENdata> *adata, const char * name) : BCModel(name), data(adata){

}

double PENModel_neutrons::LogLikelihood(const std::vector<double> & parameters) {
  double LLH = 0;

  double tau = parameters[0];
  double deteff = parameters[1];
  unsigned int nentries = data->size();
  for (unsigned int i = 0; i < nentries; i++){
	PENdata *d = &((*data)[i]);
    double N0 = parameters[2+i*2];
    double noisen = parameters[3+i*2];
    LLH += BCMath::LogPoisson(d->monitorcount, N0 + noisen*d->monitortime);
	double decay = std::exp(-d->storagetime/tau); // fraction of stored neutrons that have survived after storage time
	LLH += BCMath::LogPoisson(d->neutroncount, deteff*N0*decay + noisen*d->neutrontime);
	LLH += BCMath::LogPoisson(d->neutronbackground, noisen*d->bgtime);
  }
//  std::cout << Nn << "n, " << Np << "p, " << tau << "s, " << LLH << "logp\n";

	return LLH;
}


PENModel_protons::PENModel_protons(std::vector<PENdata> *adata, const char * name) : BCModel(name), data(adata){

}

double PENModel_protons::LogLikelihood(const std::vector<double> & parameters) {
  double LLH = 0;

  double tau = parameters[0];
  unsigned int nentries = data->size();
  for (unsigned int i = 0; i < nentries; i++){
	PENdata *d = &((*data)[i]);
    double N0 = parameters[1+i*2];
    //double noisen = parameters[2+i*3];
    double noisep = parameters[2+i*2];
    //LLH += BCMath::LogPoisson(d->monitorcount, N0 + noisen*d->monitortime);
	for (unsigned int j = 0; j < d->protonbins; j++){
	  double t = d->pbinwidth*j;
	  double rate = std::exp(-t/tau) - std::exp(-(t + d->pbinwidth)/tau); // integrate decay rate over bin width
	  //      std::cout << pd.protoncounts[j] << "n vs. " << Np*decay << '\n';
	  LLH += BCMath::LogPoisson(d->protoncounts[j], N0*rate + noisep*d->pbinwidth);
	}
	//LLH += BCMath::LogPoisson(d->neutronbackground, noisen*d->bgtime);
	LLH += BCMath::LogPoisson(d->protonbackground, noisep*d->bgtime);
  }
//  std::cout << Nn << "n, " << Np << "p, " << tau << "s, " << LLH << "logp\n";

	return LLH;
}

