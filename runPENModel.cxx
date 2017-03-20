// ***************************************************************
// This file was created using the bat-project script
// for project PENBAT.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <iostream>
#include <sstream>
#include <thread>

#include "PENModel.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"


void RunModel(BCModel *m){
	  m -> MarginalizeAll(BCIntegrate::kMargMetropolis);

	  m -> FindMode();
	  std::cout << "tau = " << m->GetBestFitParameters()[0] << " +/- " << m->GetBestFitParameterErrors()[0] << '\n';
	  delete m;
}


int main(int argc, char **argv)
{

  for (int i= 0; i < argc; i ++)
	  std::cout << i << ": " << argv[i] << '\n';

  if (argc < 2)
	  std::cout << "You have to give me a data file!\n";

  // set nicer style for drawing than the ROOT default
  BCAux::SetStyle();
  BCLog::SetLogLevelScreen(BCLog::summary);

  std::cout << "Opening " << argv[1] << "\n";
  TFile f(argv[1], "READ");
  f.cd();
  TTree *t = (TTree*)f.Get("PENdata");
  std::cout << "Found tree " << t << '\n';

  PENdata pd, pd_ideal;
  double tau, neutdeteff, monitoreff, protdeteff,/* Nn0,*/ noise_n, noise_p/*, reactorfluctuation, cleaningfluctuation, neutdetbgfluctuation, protdetbgfluctuation, lossfraction, losstau*/;
  t->SetBranchAddress("PENdata", &pd);
  t->SetBranchAddress("PENdata_ideal", &pd_ideal);
  t->SetBranchAddress("tau", &tau);
  t->SetBranchAddress("neutdeteff", &neutdeteff);
  t->SetBranchAddress("monitoreff", &monitoreff);
  t->SetBranchAddress("protdeteff", &protdeteff);
//  t->SetBranchAddress("Nn0", &Nn0);
  t->SetBranchAddress("noise_n", &noise_n);
  t->SetBranchAddress("noise_p", &noise_p);
//  t->SetBranchAddress("reactorfluctuation", &reactorfluctuation);
//  t->SetBranchAddress("cleaningfluctuation", &cleaningfluctuation);
//  t->SetBranchAddress("neutdetbgfluctuation", &neutdetbgfluctuation);
//  t->SetBranchAddress("protdetbgfluctuation", &protdetbgfluctuation);
//  t->SetBranchAddress("lossfraction", &lossfraction);
//  t->SetBranchAddress("losstau", &losstau);


  int nentries = t->GetEntries();
  std::cout << "Found data from " << nentries << " runs\n";
  std::vector<PENdata> data, data_ideal;
  for (int i = 0; i < nentries; i++){
	t->GetEntry(i);
	data.push_back(pd);
	data_ideal.push_back(pd_ideal);
  }

  f.Close();

//  std::thread *thr[3];
  for (int model = 0; model < 6; model++){
	  if (argc > 2)
		  model = atoi(argv[2]) % 6;
	  BCModel *m;
	  if (model == 0)
		  m = new PENModel_combined(&data);
	  else if (model == 1)
		  m = new PENModel_neutrons(&data);
	  else if (model == 2)
		  m = new PENModel_protons(&data);
	  if (model == 3)
		  m = new PENModel_combined(&data_ideal, "PENModel_combined_ideal");
	  else if (model == 4)
		  m = new PENModel_neutrons(&data_ideal, "PENModel_neutrons_ideal");
	  else if (model == 5)
		  m = new PENModel_protons(&data_ideal, "PENModel_protons_ideal");

	  //m->SetNIterationsPreRunMax(20000);
	  m->SetNIterationsRun(25000);
	  m->SetProposeMultivariate(false);

	  m -> AddParameter("tau", tau-20, tau+20);
	  m -> SetPriorConstant("tau");

	  if (model == 1 || model == 4){
		  m -> AddParameter("deteff", 0, neutdeteff/monitoreff*2);
		  m -> SetPriorConstant("deteff");
	  }
	  if (model == 0 || model == 3){
		  m -> AddParameter("protdeteff", 0, protdeteff/neutdeteff*2);
		  m -> SetPriorConstant("protdeteff");
	  }
	  for (unsigned int i = 0; i < data.size(); i++){
		TString pname = TString::Format("N%03d", i);
		m->AddParameter(pname.Data(), 1e4, 1e9);
		m->SetPriorConstant(pname.Data());

		if (model != 2 && model != 5){
			pname = TString::Format("noisen%03d", i);
			m -> AddParameter(pname.Data(), 0, noise_n*2);
			m -> SetPriorConstant(pname.Data());
		}

		if (model != 1 && model != 4){
			pname = TString::Format("noisep%03d", i);
			m -> AddParameter(pname.Data(), 0, noise_p*2);
			m -> SetPriorConstant(pname.Data());

		}
	  }
	  m -> WriteMarkovChain(argv[1], "UPDATE", true, false);
	  RunModel(m);

	  if (argc > 2)
		  break;

//	  thr[model] = new std::thread(RunModel, m);
  }
/*  for (int i = 0; i < 3; i++){
	  thr[i]->join();
	  delete thr[i];
  }
*/
  return 0;
}
