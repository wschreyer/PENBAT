#include <vector>
#include <iostream>
#include <chrono>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include <boost/math/special_functions/expint.hpp>

#include "PENModel.h"

double protdetrate(double *x, double *p){ // params: e_p * N0, tau, lossfrac, losstau, deadtime, noise_p
	double decayrate = p[0]/p[1]*(exp(-x[0]/p[1]) + p[2]*exp(-x[0]/p[3]));
	return exp(-p[4]*(decayrate + p[5]))*(decayrate + p[5]);
}

int main(int argc, char **argv){
	if (argc < 2){
		std::cout << "You have to give me a file name!\n";
		return -1;
	}
	std::cout << "Opening file...\n";
	TFile f(argv[1], "RECREATE");
	f.mkdir("histograms");

	std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
	int64_t seed = start.time_since_epoch().count();
	std::cout << "Creating RNG, seed: " << seed << "\n";
	TRandom3 r(seed);

	const unsigned int N = 3; // number of different storage times
	double StorageTimes[N] = {900, 1800, 2700}; // list of different storage times
	unsigned int NumberOfRuns[N] = { 10,   10,   10}; // number of runs for each storage time
	const double monitortime = 20;
	const double neutrontime = 300;
	const double bgtime = 200;

	std::cout << "Creating tree...\n";
	f.cd();
	TTree t("PENdata", "PENdata");
	std::map<std::string, double> params;
	params["tau"] = 880;
	params["neutdeteff"] = 0.3;
	params["monitoreff"] = 0.3;
	params["protdeteff"] = 0.15;
	params["Nn0"] = 2e7;
	params["noise_n"] = 1;
	params["noise_p"] = 40;
	params["reactorfluctuation"] = 0.0;
	params["cleaningfluctuation"] = 0.0;
	params["neutdetbgfluctuation"] = 0.0;
	params["protdetbgfluctuation"] = 0.0;
	params["protdeadtime"] = 0;
	params["pbinwidth"] = 100;
	params["lossfraction"] = 0; //pow(10, r.Uniform(-5, -1));
	params["losstau"] = 880; //r.Uniform(0, 880);

	params["neutrontime"] = 300;
	params["monitortime"] = 20;
	params["bgtime"];
	for (int i = 0; i < N; i++){
		std::ostringstream s;
		s << "storagetime" << i;
		params[s.str()] = StorageTimes[i];
		std::ostringstream n;
		n << "numberofruns" << i;
		params[n.str()] = NumberOfRuns[i];
	}

	int job = 0;
	if (argc > 2)
		job = atoi(argv[2]);
	if (job >= 1 && job < 20)
		params["neutdeteff"] = 0.05001*(job);
	else if (job >= 20 && job < 40)
		params["monitoreff"] = 0.05001*(job - 19);
	else if (job >= 40 && job < 50)
		params["protdeteff"] = 0.05001*(job - 39);
	else if (job >= 50 && job < 60)
		params["Nn0"] = 2.00001e7*pow(2., job - 55);
	else if (job >= 60 && job < 70)
		params["noise_n"] = 1.00001*pow(2., job - 60);
	else if (job >= 70 && job < 80)
		params["noise_p"] = 40.0001*pow(2., job - 75);
	else if (job >= 80 && job < 90)
		params["reactorfluctuation"] = 1e-9 + 0.025*(job - 80);
	else if (job >= 90 && job < 100)
		params["neutdetbgfluctuation"] = 1e-9 + 0.025*(job - 90);
	else if (job >= 100 && job < 110)
		params["protdetbgfluctuation"] = 1e-9 + 0.025*(job - 100);
	else if (job >= 110 && job < 120)
		params["cleaningfluctuation"] = 1e-9 + 0.0001*(job - 110);
	else if (job >= 120 && job < 130)
		params["protdeadtime"] = 10e-6*pow(2., job - 123)/400; // dead time (400 channels)
	else if (job >= 130 && job < 140)
		params["pbinwidth"] = 45*(job - 129);
	else if (job >= 140){
		params["lossfraction"] = pow(10., r.Uniform(-5,-1));
		params["losstau"] = 1./(1./880. + pow(10., r.Uniform(-8, -1)));
	}

	std::vector<double> vparams(params.size());
	for (std::map<std::string, double>::iterator i = params.begin(); i != params.end(); i++){
		int idx = std::distance(params.begin(), i);
		vparams[idx] = i->second;
		t.Branch(i->first.c_str(), &vparams[idx]);
	}

	PENdata pd, pd_ideal;
	TBranch *b = t.Branch("PENdata", &pd, "storagetime/D:pbinwidth/D:monitortime/D:neutrontime/D:bgtime/D:monitorcount/i:neutroncount/i:neutronbackground/i:protonbackground/i:protonbins/i:protoncounts[protonbins]/i");
	TBranch *b_ideal = t.Branch("PENdata_ideal", &pd_ideal, "storagetime/D:pbinwidth/D:monitortime/D:neutrontime/D:bgtime/D:monitorcount/i:neutroncount/i:neutronbackground/i:protonbackground/i:protonbins/i:protoncounts[protonbins]/i");
	for (unsigned int i = 0; i < N; i++){
		pd.storagetime = pd_ideal.storagetime = StorageTimes[i];
		pd.monitortime = pd_ideal.monitortime = monitortime;
		pd.neutrontime = pd_ideal.neutrontime = neutrontime;
		pd.bgtime = pd_ideal.bgtime = bgtime;
		std::cout << NumberOfRuns[i] << " x Storage time " << pd.storagetime << "\n";
		for (unsigned int j = 0; j < NumberOfRuns[i]; j++){
			double noise_n = r.Uniform(params["noise_n"]*(1 - params["neutdetbgfluctuation"]), params["noise_n"]*(1 + params["neutdetbgfluctuation"]));
			double noise_p = r.Uniform(params["noise_p"]*(1 - params["protdetbgfluctuation"]), params["noise_p"]*(1 + params["protdetbgfluctuation"])); // neutron and proton noise fluctuates between storage cycles
			unsigned int Nn = r.Uniform(params["Nn0"]*(1 - params["reactorfluctuation"]), params["Nn0"]*(1 + params["reactorfluctuation"])); // mean number of neutrons fluctuates with UCN source intensity

			pd.monitorcount = r.Poisson(params["monitoreff"]*Nn) + r.Poisson(noise_n*monitortime); // seen by monitor detector
			pd_ideal.monitorcount = params["monitoreff"]*Nn + noise_n*monitortime;

			pd.neutronbackground = r.Poisson(noise_n*bgtime); // seen during background measurement
			pd_ideal.neutronbackground = noise_n*bgtime;

			pd.protonbackground = r.Poisson(noise_p*bgtime*exp(-noise_p*params["protdeadtime"])); // detected proton counts during background measurement with reduction due to dead time
			pd_ideal.protonbackground = noise_p*bgtime*exp(-noise_p*params["protdeadtime"]);

			//Nn = r.Uniform(Nn*(1 - params["cleaningfluctuation"]), Nn*(1 + params["cleaningfluctuation"])); // stored density fluctuates during cleaning (not seen by monitor detector)
			if (StorageTimes[i] == 900)
				Nn *= 1 - params["cleaningfluctuation"];
			else if (StorageTimes[i] == 2700)
				Nn *= 1 + params["cleaningfluctuation"];

			double decay = exp(-pd.storagetime/params["tau"]) + params["lossfraction"]*exp(-pd.storagetime/params["losstau"]);

			pd.neutroncount = r.Poisson(params["neutdeteff"]*Nn*decay) + r.Poisson(noise_n*neutrontime); // seen by lifetime detector
			pd_ideal.neutroncount = params["neutdeteff"]*Nn*decay + noise_n*neutrontime;

			if (pd.storagetime > 0){
				pd.protonbins = pd_ideal.protonbins = 10; //pd.storagetime/pd.pbinwidth;
				pd.pbinwidth = pd_ideal.pbinwidth = pd.storagetime/pd.protonbins;//params["pbinwidth"];

				std::ostringstream name, title;
				name << pd.storagetime << "s_" << j;
				title << name.str() << " (M: " << pd.monitorcount << ", N: " << pd.neutroncount << ")";
				TH1I *h = new TH1I(name.str().c_str(), title.str().c_str(), pd.protonbins, 0, pd.storagetime);

				//std::cout << pd.storagetime << "s: " << pd.neutroncount << "n, " << pd.protonbins << "bins, ";
				TF1 rate("rate", &protdetrate, 0, pd.storagetime, 6);
				rate.SetParameters(Nn*params["protdeteff"], params["tau"], params["lossfraction"], params["losstau"], params["protdeadtime"], noise_p);
				for (unsigned int k = 0; k < pd.protonbins; k++){
					double time = pd.pbinwidth*k;
					double counts = rate.Integral(time, time + pd.pbinwidth);
					pd.protoncounts[k] = r.Poisson(counts);
					pd_ideal.protoncounts[k] = counts;

					h->Fill(time + pd.pbinwidth/2, pd.protoncounts[k]);
				}
				//      cout << '\n';

				//      cout << "Drawing histogram...\n";
				f.cd("histograms");
				h->Write();
				delete h;
			}
			//      cout << "Filling tree...\n";
			f.cd();
			t.Fill();
		}
	}
	f.Write();
	f.Close();
}

