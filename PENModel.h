// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__PENMODEL__H
#define __BAT__PENMODEL__H

#include <BAT/BCModel.h>

// This is a PENModel header file.
// Model source code is located in file PENBAT/PENModel.cxx

const int MaxT = 10000; // max. possible storage time
struct PENdata{
  double storagetime;
  double pbinwidth;
  double monitortime;
  double neutrontime;
  double bgtime;
  unsigned int monitorcount;
  unsigned int neutroncount;
  unsigned int neutronbackground;
  unsigned int protonbackground;
  unsigned int protonbins;
  unsigned int protoncounts[MaxT];
};


class PENModel_combined : public BCModel {
 public:
	// Constructor and destructor
	PENModel_combined(std::vector<PENdata> *adata, const char * name = "PENModel_combined");

	// Methods to overload, see file PENModel.cxx
	double LogLikelihood(const std::vector<double> & parameters);
	// double LogAPrioriProbability(const std::vector<double> & parameters);

	std::vector<PENdata> *data;
};


class PENModel_neutrons : public BCModel {
 public:
	// Constructor and destructor
	PENModel_neutrons(std::vector<PENdata> *adata, const char * name = "PENModel_neutrons");

	// Methods to overload, see file PENModel.cxx
	double LogLikelihood(const std::vector<double> & parameters);
	// double LogAPrioriProbability(const std::vector<double> & parameters);

	std::vector<PENdata> *data;
};


class PENModel_protons : public BCModel {
 public:
	// Constructor and destructor
	PENModel_protons(std::vector<PENdata> *adata, const char * name = "PENModel_protons");

	// Methods to overload, see file PENModel.cxx
	double LogLikelihood(const std::vector<double> & parameters);
	// double LogAPrioriProbability(const std::vector<double> & parameters);

	std::vector<PENdata> *data;
};



#endif
