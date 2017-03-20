#include <iostream>

#include "TSystemDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TGraph2DErrors.h"

int main(int argc, char **argv){
	if (argc < 2){
		std::cout << "You have to give me a data folder!\n";
		return -1;
	}

	TSystemDirectory dir(argv[1], argv[1]); // open given folder
	TList *files = dir.GetListOfFiles(); // get all files in that folder

	std::string modelname[6] = {"neutrons", "protons", "combined", "neutrons_ideal", "protons_ideal", "combined_ideal"};
	double tau_n[6], dtau_n[6];
	TFile outfile("results.root", "RECREATE");
	TTree outtree("PENresults", "PENresults");

	outtree.Branch("tau_n", &tau_n[0], "tau_n/D");
	outtree.Branch("dtau_n", &dtau_n[0], "dtau_n/D");
	outtree.Branch("tau_p", &tau_n[1], "tau_p/D");
	outtree.Branch("dtau_p", &dtau_n[1], "dtau_p/D");
	outtree.Branch("tau_c", &tau_n[2], "tau_c/D");
	outtree.Branch("dtau_c", &dtau_n[2], "dtau_c/D");
	outtree.Branch("tau_ni", &tau_n[3], "tau_ni/D");
	outtree.Branch("dtau_ni", &dtau_n[3], "dtau_ni/D");
	outtree.Branch("tau_pi", &tau_n[4], "tau_pi/D");
	outtree.Branch("dtau_pi", &dtau_n[4], "dtau_pi/D");
	outtree.Branch("tau_ci", &tau_n[5], "tau_ci/D");
	outtree.Branch("dtau_ci", &dtau_n[5], "dtau_ci/D");

	std::vector<std::string> paramnames;

	for (Int_t i = 0; ; i++){ // for loop incrementing index i
		//std::cout << "File " << i << std::endl;
		TObject *f = files->At(i); // get file from folder with index i
		//std::cout << f << std::endl;
		if (f){ // if next file was found
			TString filename = f->GetName(); // get filename
			std::cout << "Found " << filename << std::endl;
			if (filename.BeginsWith("PENdata_")){
				TFile datafile(filename, "READ");
				TTree *data = (TTree*)datafile.Get("PENdata");
				if (!data){
					std::cout << "Tree PENdata not found in " << filename << '\n';
					datafile.Close();
					continue;
				}
				TObjArray *parambranches = data->GetListOfBranches();
				int nparams = parambranches->GetEntries();
				//std::cout << "Found " << nparams << " parameters\n";
				std::vector<double> params(nparams);
				for (int j = 0; j < nparams; j++){
					datafile.cd();
					TBranch *parambranch = (TBranch*)parambranches->At(j);
					TString pname = parambranch->GetName();
					if (pname.BeginsWith("PENdata"))
						continue;
					//std::cout << "Found parameter " << pname << std::endl;
					data->SetBranchAddress(pname.Data(), &params[j]);
					outfile.cd();
					TBranch *outparam = outtree.GetBranch(pname.Data());
					if (!outparam){
						outparam = outtree.Branch(pname.Data(), &params[j], (pname + "/D").Data());
						//std::cout << "Did not find out-branch - creating " << outparam->GetName() << std::endl;
						paramnames.push_back(pname.Data());
					}
					else{
						outtree.SetBranchAddress(pname.Data(), &params[j]);
						//std::cout << "Found out-branch\n";
					}
				}
				datafile.cd();
				data->GetEntry();
				//for (int j = 0; j < nparams; j++)
				//	std::cout << params[j] << " ";
				//std::cout << std::endl;
				for (int i = 0; i < 6; i++){
					tau_n[i] = 880;
					dtau_n[i] = 0;
					TString treename("PENModel_");
					treename += modelname[i] + "_mcmc";
					TTree *results = (TTree*)datafile.Get(treename.Data());
					if (!results){
						std::cout << "Tree " << treename.Data() << " not found in " << filename << std::endl;
						continue;
					}
					//std::cout << "Tree " << treename.Data() << " found in " << filename << std::endl;
					TH1I result("tau", "tau", 100, 870, 890);
					results->Draw("tau>>tau");
					//TFitResultPtr fit = result.Fit("gaus", "SQ");
					//std::cout << "Fit result " << fit << std::endl;
					//if ((int)fit == 0){
						tau_n[i] = result.GetMean();//fit->Parameter(1);
						dtau_n[i] = result.GetRMS();//fit->Parameter(2);
						std::cout << tau_n[i] << " +/- " << dtau_n[i] << '\n';
					//}
				}
				//std::cout << "Closing" << std::endl;
				datafile.Close();

				outfile.cd();
				//std::cout << "Filling" << std::endl;
				outtree.Fill();

				//std::cout << "Done with " << filename << std::endl;
			}
		}
		else // if no more files were found
			break;
	}

	outfile.Write();

	for (auto i = paramnames.begin(); i != paramnames.end(); i++){
		std::cout << i->c_str() << std::endl;

		std::string sel;
		bool logx = false, logy = false;
		if (*i == "Nn0"){
			sel = "2e7";
			logx = logy = true;
		}
		if (*i == "cleaningfluctuation") sel = "0";
		if (*i == "lossfraction"){
			sel = "0";
			logx = true;
		}
		if (*i == "losstau") sel = "880";
		if (*i == "monitoreff") sel = "0.3";
		if (*i == "neutdeadtime") continue; //sel = "0";
		if (*i == "neutdetbgfluctuation") sel = "0";
		if (*i == "neutdeteff") sel = "0.3";
		if (*i == "noise_n"){
			sel = "1";
			logx = true;
		}
		if (*i == "noise_p"){
			sel = "40";
			logx = true;
		}
		if (*i == "protdeadtime"){
			sel = "0";
			logx = true;
		}
		if (*i == "protdetbgfluctuation") sel = "0";
		if (*i == "protdeteff") sel = "0.15";
		if (*i == "pbinwidth") sel = "100";
		if (*i == "reactorfluctuation") sel = "0";
		if (*i == "tau") continue; //sel = "880";

		TMultiGraph mg(i->c_str(), i->c_str());
		TGraph *gr;
		TCanvas *c;

		c = new TCanvas("c","c");
		std::string str = "dtau_c:" + *i;
		outtree.Draw(str.c_str(), (*i + "!=" + sel).c_str());
		gr = (TGraph*)c->GetPrimitive("Graph")->Clone();
		gr->SetMarkerStyle(20);
		gr->SetMarkerColor(kBlack);
		mg.Add(gr);
		delete c;

		c = new TCanvas("c","c");
		str = "dtau_p:" + *i;
		outtree.Draw(str.c_str(), (*i + "!=" + sel).c_str());
		gr = (TGraph*)c->GetPrimitive("Graph")->Clone();
		gr->SetMarkerStyle(20);
		gr->SetMarkerColor(kRed);
		mg.Add(gr);
		delete c;

		c = new TCanvas("c","c");
		str = "dtau_n:" + *i;
		outtree.Draw(str.c_str(), (*i + "!=" + sel).c_str());
		gr = (TGraph*)c->GetPrimitive("Graph")->Clone();
		gr->SetMarkerStyle(20);
		gr->SetMarkerColor(kBlue);
		mg.Add(gr);
		delete c;

		mg.Write();

		c = new TCanvas("c", "c");
		if (logx) c->SetLogx();
		if (logy) c->SetLogy();
		mg.Draw("AP");
		c->Print((*i + ".pdf").c_str());
		c->Print((*i + ".C").c_str());
		delete c;


		TMultiGraph mgi(i->c_str(), i->c_str());

		c = new TCanvas("c","c");
		str = "abs(tau_ci):" + *i;
		outtree.Draw(str.c_str(), (*i + "!=" + sel).c_str());
		gr = (TGraph*)c->GetPrimitive("Graph")->Clone();
		gr->SetMarkerStyle(20);
		gr->SetMarkerColor(kBlack);
		mgi.Add(gr);
		delete c;

		c = new TCanvas("c","c");
		str = "abs(tau_pi):" + *i;
		outtree.Draw(str.c_str(), (*i + "!=" + sel).c_str());
		gr = (TGraph*)c->GetPrimitive("Graph")->Clone();
		gr->SetMarkerStyle(20);
		gr->SetMarkerColor(kRed);
		mgi.Add(gr);
		delete c;

		c = new TCanvas("c","c");
		str = "abs(tau_ni):" + *i;
		outtree.Draw(str.c_str(), (*i + "!=" + sel).c_str());
		gr = (TGraph*)c->GetPrimitive("Graph")->Clone();
		gr->SetMarkerStyle(20);
		gr->SetMarkerColor(kBlue);
		mgi.Add(gr);
		delete c;

		mg.Write();

		c = new TCanvas("c", "c");
		if (logx) c->SetLogx();
		if (logy) c->SetLogy();
		mgi.Draw("AP");
		c->Print((*i + "i.pdf").c_str());
		c->Print((*i + "i.C").c_str());
		delete c;

	}

	outfile.Close();

	return 0;
}
