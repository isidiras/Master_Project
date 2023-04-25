//With this script I am presenting my analysis at E_CM = 14 TeV,  E_CM 500 GeV and E_CM = 200 GeV
//and following similar analysis with  Shushu Shi et al(2015) https://arxiv.org/pdf/1507.00614.pdf 
//but using SoftQCD settings and MONASH tune. The script is designed that way so the analysis
//can be done for other hadrons.


//C++ libraries
#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
//Root Libraries
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"

#define PI 3.14159265
using namespace std;
using namespace TMath;

Double_t DeltaPhi(Double_t phi1, Double_t phi2){
	//returns delta phi in range -pi/2 3pi/2
	return fmod(phi1-phi2+2.5*PI,2*PI)-0.5*PI;
	}

void correlation_plot(Int_t id_trigger,Int_t id_associate, TString filename, const char* title){
	//This function takes the pdg ids of the triger and associate particle as inputs and produces all kinds of
	//correlation plots concerning the DeltaPhi and DeltaEta as a function of pT.
	//filename: filename of root file that the histograms are going to be stored
	//title: Part of the histograms title concerning the particle for example it can be something like "D^{+}D^{+}"
	
	//Define the TChain
	TChain *ch1 = new TChain("tree");
	TFile *output = new TFile(filename,"RECREATE");
	
	int ntrees = 30; //Number of trees we want to add to the TChain
	
	//We put the trees to the chain
	for( int i = 1; i < ntrees+1;  i++){
	//One can change the file path accordingly to his/hers machine set up
	ch1->Add(Form("/home/isidiras/university_staff/ccbar_MONASH/MONASH_BATCH1/Group%i/output.root",i));
	}
	
	//Now we define the vectors and set branch addresses remember each vector corresponds to an event
	vector<Int_t>* vID = 0;
	vector<Double_t>* vPt = 0;//Carefull of null pointers.
	vector<Double_t>* vEta = 0;
	vector<Double_t>* vPhi = 0;
	vector<Double_t>* vCharge = 0;
	//Setting up chain addresses
	ch1->SetBranchAddress("ID",&vID);
	ch1->SetBranchAddress("PT",&vPt);
	ch1->SetBranchAddress("ETA",&vEta);
	ch1->SetBranchAddress("PHI",&vPhi);
	ch1->SetBranchAddress("CHARGE",&vCharge);
	
	//Here we define the variables for each particle
	Int_t pid,aid;
	int nTrigger = 0;
	Double_t pPt,pEta,pPhi,pCharge;//For triger
	Double_t aPt,aEta,aPhi,aCharge;//For associate
	Int_t trigL,trigI,trigH,pairLL,pairIL, pairII, pairHL, pairHI, pairHH;//Trigger particle counter variable
	
	int nEvents = ch1->GetEntries();
	
	cout<<"The number of events for this analysis is: "<<nEvents<<endl; 
	
	//Setting up histograms
	//General simulation variables distribution
	TH1D *hPt = new TH1D("hPt","Transverse momentum distribution for D^{0};(Gev/c);Counts",100,0,50);
	TH1D *hEta = new TH1D("hEta","#eta distribution of all the charm particle  produced;#eta;Counts",40,-4,4);
	//Multidimensional correlation plot
	TH2D *hDPhiDEta = new TH2D("hDPhiDEta",Form("%s #Delta#phi and #Delta#eta;#Delta#phi (rad);#Delta#eta",title),100,-PI/2,3*PI/2,30,-3,3);
	//DeltaPhi correlation plots
	TH1D* hDeltaPhi = new TH1D("hDeltaPhi",Form("%s correlations;#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	
	//Begining of event loop
	for(int iEvent = 0; iEvent < nEvents; iEvent++){
		ch1->GetEntry(iEvent);
		int nparticles = vID->size();
		//Starting particle loop
		for(int iparticle = 0; iparticle < nparticles; iparticle++){
			pid = (*vID)[iparticle];
			pPhi = (*vPhi)[iparticle];
			pPt = (*vPt)[iparticle];
			pEta = (*vEta)[iparticle];
			pCharge = (*vCharge)[iparticle];
			//Trigger loop
			if( pid == id_trigger && pPt >= 2. && abs(pEta) < 1 ){
				nTrigger++;
				hEta->Fill(pEta);
				hPt->Fill(pPt);
				for(int jparticle = 0; jparticle < nparticles; jparticle++){
					if (jparticle == iparticle) continue;//So we do not correlate with itself.
						aid = (*vID)[jparticle];
						aPhi = (*vPhi)[jparticle];
						aEta = (*vEta)[jparticle];
						aPt = (*vPt)[jparticle];
						if(aid == id_associate && abs(aEta) < 1){
							//Filling general correlationplots
							hDeltaPhi->Fill(DeltaPhi(pPhi,aPhi));
							hDPhiDEta->Fill(DeltaPhi(pPhi,aPhi),pEta-aEta);
					}//End of DD correlation
				}//End of associate particle loop D+D+
			}//D+ triger
			
		}//End of  particle loop
		
	}//End of event loop
	if(nTrigger == 0){
		cout<<"Have not found any trigger particle with id: "<<id_trigger<<endl;
		output->Close();
		return;
	}
	output->Write();
	output->Close();
	cout<<"The total number of triggers is: "<<nTrigger<<endl;
	
	cout<<"File: "<<filename<<" has been created!"<<endl;
	
}//End of correlation plot function


TTree* read_tree(TString input_file){
	TFile *input = new TFile(input_file,"READ");
	TTree *tree = (TTree*)input->Get("tree");
	return tree;
}

void correlation_tree(Int_t id_trigger,Int_t id_associate,TString filename, const char* title){
//This function does exacly what the above one does but it reads the TTree created specific
//for this analysis.


	TTree* tree = read_tree("/home/isidiras/university_staff/hf291122/D0200.root");
	//Now we define the vectors and set branch addresses remember each vector corresponds to an event
	vector<Int_t>* vID = 0;
	vector<Double_t>* vPt = 0;//Carefull of null pointers.
	vector<Double_t>* vEta = 0;
	vector<Double_t>* vPhi = 0;
	vector<Double_t>* vCharge = 0;
	//Setting up chain addresses
	tree->SetBranchAddress("ID",&vID);
	tree->SetBranchAddress("PT",&vPt);
	tree->SetBranchAddress("ETA",&vEta);
	tree->SetBranchAddress("PHI",&vPhi);
	tree->SetBranchAddress("CHARGE",&vCharge);
	
	Int_t pid,aid;
	int nTrigger = 0;
	Double_t pPt,pEta,pPhi,pCharge;//For triger
	Double_t aPt,aEta,aPhi,aCharge;//For associate
	Int_t trigL,trigI,trigH,pairLL,pairIL, pairII, pairHL, pairHI, pairHH;//Trigger particle counter variable
	
	int nEvents = tree->GetEntries();
	
	cout<<"The number of events for this analysis is: "<<nEvents<<endl; 
	TFile *output = new TFile(filename,"RECREATE");
	//Setting up histograms
	//General simulation variables distribution
	TH1D *hPt = new TH1D("hPt","Transverse momentum distribution for D^{0};(Gev/c);Counts",100,0,50);
	TH1D *hEta = new TH1D("hEta","#eta distribution of all the charm particle  produced;#eta;Counts",40,-4,4);
	//Multidimensional correlation plot
	TH2D *hDPhiDEta = new TH2D("hDPhiDEta",Form("%s #Delta#phi and #Delta#eta;#Delta#phi (rad);#Delta#eta",title),100,-PI/2,3*PI/2,30,-3,3);
	//DeltaPhi correlation plots
	TH1D* hDeltaPhi = new TH1D("hDeltaPhi",Form("%s correlations;#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	
	//Begining of event loop
	for(int iEvent = 0; iEvent < nEvents; iEvent++){
		tree->GetEntry(iEvent);
		int nparticles = vID->size();
		//Starting particle loop
		for(int iparticle = 0; iparticle < nparticles; iparticle++){
			pid = (*vID)[iparticle];
			pPhi = (*vPhi)[iparticle];
			pPt = (*vPt)[iparticle];
			pEta = (*vEta)[iparticle];
			pCharge = (*vCharge)[iparticle];
			
			//Trigger loop
			if( pid == id_trigger && pPt > 2.){
				nTrigger++;
				hEta->Fill(pEta);
				hPt->Fill(pPt);
				for(int jparticle = 0; jparticle < nparticles; jparticle++){
					if (jparticle == iparticle) continue;//So we do not correlate with itself.
						aid = (*vID)[jparticle];
						aPhi = (*vPhi)[jparticle];
						aEta = (*vEta)[jparticle];
						aPt = (*vPt)[jparticle];
						if(aid == id_associate){
							//Filling general correlationplots
							hDeltaPhi->Fill(DeltaPhi(pPhi,aPhi));
							hDPhiDEta->Fill(DeltaPhi(pPhi,aPhi),pEta-aEta);
					}//End of DD correlation
				}//End of associate particle loop D+D+
			}//D+ triger
			
		}//End of  particle loop
		
	}//End of event loop
	if(nTrigger == 0){
		cout<<"Have not found any trigger particle with id: "<<id_trigger<<endl;
		output->Close();
		return;
	}
	output->Write();
	output->Close();
	cout<<"The total number of triggers is: "<<nTrigger<<endl;
	
	cout<<"File: "<<filename<<" has been created!"<<endl;
	
}//End of correlation plot function 



TH1D *Read_Hist(TString filename, TString hist_name){
	//This function read a histogram from .root files.
	TH1D *hist;
	TFile *input = TFile::Open(filename);
	if(!input->IsOpen() || !input){
		cout<<"File: "<<filename<<" is not found!"<<endl;
		cout<<"Terminating!"<<endl;
		return 0;
	}
	
	hist = (TH1D*)input->Get(hist_name);
	if(!hist){
		cout<<"Histogram: "<<hist_name<<" not found in "<<filename<<endl;
		cout<<"Terminating!"<<endl;
		return 0;
	}
	return hist;
}

void normalize(TH1D *hist, TH1D *trig){
	hist->Scale(1/trig->Integral());
}

void plot_results(TString filename1, TString filename2, TString filename3){
	//Setting up histogram arrays
	TH1D *hist[3];
	TH1D *trig[3];
	TLegend *leg = new TLegend(0.89,0.7,0.99,0.9);
	leg->SetBorderSize(0);
	//Dictionaries used
	const char* filenames[3] = {filename1,filename2,filename3};
	const char* legend[3] = {"E_{CM} = 14 TeV", "E_{CM} = 500 GeV", "E_{CM} = 200 GeV"};
	
	//Imorting histograms
	for(int i = 0; i < 3; i++){
		hist[i] = Read_Hist(filenames[i],"hDeltaPhi");
		trig[i] = Read_Hist(filenames[i],"hPt");
	}
	//Normalizing and plot cosmetics and Legend
	for(int i = 0; i < 3; i++){
		normalize(hist[i],trig[i]);
		hist[i]->SetMarkerColor(i+2);
		hist[i]->SetMarkerStyle(7);
		hist[i]->SetStats(0);
		hist[i]->GetYaxis()->SetTitle("#frac{1}{N_{tr}} #frac{dN_{assoc}}{d#Delta#phi}");
		leg->AddEntry(hist[i],legend[i],"p");
	}
	TCanvas* c1 = new TCanvas();
	hist[0]->Draw("p");
	hist[1]->Draw("p SAME");
	hist[2]->Draw("p SAME");
	leg->Draw();
	
}

void paper_verification(){
	//correlation_plot(421,-421,"D014000.root","D^{0}#barD^{0}");
	//correlation_tree(421,-421,"D02.root","D^{0}#barD^{0}");
	plot_results("D014000.root","D05.root","D02.root");
	return 0;
}

