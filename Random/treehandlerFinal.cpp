//This script contains functions that can be used to create store and plot correlation plots for varius particles.
//C++ libraries
#include <iostream>
#include <vector>
#include <cstring>
//Root Libraries
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TChain.h"

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
	
	int ntrees = 100; //Number of trees we want to add to the TChain
	
	//We put the trees to the chain
	for( int i = 1; i < ntrees+1;  i++){
	//One can change the file path accordingly to his/hers machine set up
	ch1->Add(Form("/home/isidiras/university_staff/ccbar_MONASH/output_MONASH/Group%i/output.root",i));
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
	
	int nEvents = ch1->GetEntries();
	
	cout<<"The number of events for this analysis is: "<<nEvents<<endl; 
	
	//Setting up histograms
	//General simulation variables distribution
	TH1D *hPt = new TH1D("hPt","Transverse momentum distribution of all the charm partile produced;p_{T} (Gev/c);Counts",100,0,50);
	TH1D *hEta = new TH1D("hEta","#eta distribution of all the charm particle  produced;#eta;Counts",40,-4,4);
	TH1D *hPhi = new TH1D("hPhi","Azimuthial angle distribution of all charm particle produced;#phi (rad);Counts",60,-PI,PI);
	TH1D *hID = new TH1D("hID","ID distribution of all the charm particle produced;ID;Counts",12001,-6000.5,6000.5);
	TH1D *hCharge = new TH1D("hCharge","Charge distribution of charm particles;q (e);Counts",100,-4.5,4.5);
	//Multidimensional correlation plots
	TH3D *hDPhiDEtaPt = new TH3D("hDPhiDEtaPt",Form("#Delta#phi #Delta#eta p_{T} for %s;#Delta#phi (rad);#Delta#eta;p_{T} (GeV/c)",title),100,-PI/2,3*PI/2,80,-8,8,100,0,50);
	TH2D *hDeltaPhiPt = new TH2D("hDeltaPhiPt",Form("%s #Delta#phi and trigger p_{T};p_{T} (GeV/c); #Delta#phi (rad)",title),100,0,50,100,-PI/2,3*PI/2);
	TH2D *hDEtaPt = new TH2D("hDEtaPt",Form("%s #Delta#eta and trigger p_{T};p_{T} (GeV/c);#Delta#eta",title),100,0,50,80,-8,8);
	TH2D *hDPhiDEta = new TH2D("hDPhiDEta",Form("%s #Delta#phi and #Delta#eta;#Delta#phi (rad);#Delta#eta",title),100,-PI/2,3*PI/2,80,-8,8);
	//DeltaPhi correlation plots
	TH1D* hDeltaPhi = new TH1D("hDeltaPhi",Form("%s correlations;#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhill = new TH1D("hDPhill",Form("%s #Delta#phi correlation for low p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiiil = new TH1D("hDPhiiil",Form("%s #Delta#phi correlation for intermidiate p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhih = new TH1D("hDPhih",Form("%s #Delta#phi correlation for high p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	//DeltaEta correlation plots
	TH1D* hDEta = new TH1D("hDEta",Form("%s #Delta#eta correlation;#Delta#eta;Counts",title),80,-8,8);
	TH1D* hDEtall = new TH1D("hDEtall",Form("%s #Delta#eta correlation for low p_{T};#Delta#eta;Counts",title),80,-8,8);
	TH1D* hDEtaiil = new TH1D("hDEtaiil",Form("%s #Delta#eta correlation for intermidiate p_{T};#Delta#eta;Counts",title),80,-8,8);
	TH1D* hDEtah= new TH1D("hDEtah",Form("%s #Delta#eta correlation for high p_{T};#Delta#eta;Counts",title),80,-8,8);
	
	
	//Begining of event loop
	for(int iEvent = 0; iEvent < nEvents; iEvent++){
		ch1->GetEntry(iEvent);
		int nparticles = vID->size();
		for(int iparticle = 0; iparticle < nparticles; iparticle++){
			pid = (*vID)[iparticle];
			pPhi = (*vPhi)[iparticle];
			pPt = (*vPt)[iparticle];
			pEta = (*vEta)[iparticle];
			pCharge = (*vCharge)[iparticle];
			//Filling general histograms
			hID->Fill((Double_t) pid);
			hPhi->Fill(pPhi);
			hPt->Fill(pPt);
			hEta->Fill(pEta);
			hCharge->Fill(pCharge);
			if(pid == id_trigger){
				nTrigger++;
				for(int jparticle = 0; jparticle < nparticles; jparticle++){
					if (iparticle == jparticle) continue;//So we do not correlate with itself.
						aid = (*vID)[jparticle];
						aPhi = (*vPhi)[jparticle];
						aEta = (*vEta)[jparticle];
						aPt = (*vPt)[jparticle];
						if(aid == id_associate){
							//Filling general correlationplots
							hDeltaPhi->Fill(DeltaPhi(pPhi,aPhi));
							hDeltaPhiPt->Fill(pPt,DeltaPhi(pPhi,aPhi));
							hDPhiDEtaPt->Fill(pPt,DeltaPhi(pPhi,aPhi),pEta-aEta);
							hDPhiDEta->Fill(DeltaPhi(pPhi,aPhi),pEta-aEta);
							hDEtaPt->Fill(pPt,pEta-aEta);
							hDEta->Fill(pEta-aEta);
							//Filling triger momentum range correlations
							if(pPt >= 1 && pPt< 3 && aPt < 3){
								hDPhill->Fill(DeltaPhi(pPhi,aPhi));
								hDEtall->Fill(pEta-aEta);
							}//End of low transverse momentum
							if(pPt >= 3 && pPt < 8 && aPt < 8){
								hDPhiiil->Fill(DeltaPhi(pPhi,aPhi));
								hDEtaiil->Fill(pEta-aEta);								
							}//End of intermidiate transverse momentum
							if(pPt >= 8){
								hDPhih->Fill(DeltaPhi(pPhi,aPhi));
								hDEtah->Fill(pEta-aEta);
							}//End of high transverse momentum
							
					}//End of DD correlation
				}//End of associate particle loop D+D+
			}//D+ triger
		}//End of triger particle loop
	}//End of event loop
	if(nTrigger == 0){
		cout<<"Have not found any trigger particle with id: "<<id_trigger<<endl;
		output->Close();
		return;
	}
	output->Write();
	output->Close();
	
	cout<<"File: "<<filename<<" has been created!"<<endl;
	
}//End of correlation plot function

void draw_correlations(TString filename){
	//This function will take the root file with histograms we created using correlation_plot and draw the histograms
	TFile *input = TFile::Open(filename);
	if(!input->IsOpen() || !input){
		cout<<"File: "<<filename<<" is not found!"<<endl;
		cout<<"Terminating!"<<endl;
		return;
	}
	
	//Creating histogram arrays that I am going to use
	TH1D *OneDcorr[8];
	TH2D *TwoDcorr[3];
	
	//Arrays of histogram names
	const char* OneDnames[8] = {"hDeltaPhi","hDPhill","hDPhiiil","hDPhih","hDEta","hDEtall","hDEtaiil","hDEtah"};
	const char* TwoDnames[3] = {"hDeltaPhiPt","hDEtaPt","hDPhiDEta"};
	
	//Here we are importing the histograms
	for(int i = 0; i<8; i++){
		OneDcorr[i] = (TH1D*)input->Get(OneDnames[i]);
		if(!OneDcorr[i]){
			cout<<OneDnames[i]<<" not found in "<<filename<<endl;
			cout<<"Terminating!"<<endl;
			return;
		}
	}
	
	for(int i = 0; i<3; i++){
		TwoDcorr[i] = (TH2D*)input->Get(TwoDnames[i]);
		if(!TwoDcorr[i]){
			cout<<TwoDnames[i]<<" not found in "<<filename<<endl;
			cout<<"Terminating!"<<endl;
			return;
		}
	}
	
	TH3D* DPhiDEtaPt = (TH3D*)input->Get("hDPhiDEtaPt");
	if(!DPhiDEtaPt){
		cout<<"3D histogram not found in "<<filename<<endl;
		cout<<"Terminating"<<endl;
		return;
	}
	
	//Plot looks
	for(int i = 0; i<8;i++){
		OneDcorr[i]->SetMarkerStyle(7);
	}
	//Drawing in canvases
	//3D
	TCanvas *c3D = new TCanvas();
	DPhiDEtaPt->SetStats(0);
	DPhiDEtaPt->Draw();
	
	//2D
	TCanvas *c2D[3];
	gStyle->SetPalette(kRainBow);
	for(int i = 0; i < 3; i++){
		c2D[i] = new TCanvas(TwoDcorr[i]->GetName(),TwoDcorr[i]->GetName(),1920,1080);
		gStyle->SetPadLeftMargin(0.08); gStyle->SetPadRightMargin(0.2);
		TwoDcorr[i]->GetZaxis()->SetTitle("Counts");
		TwoDcorr[i]->SetContour(1000);
		TwoDcorr[i]->SetStats(0);
		TwoDcorr[i]->Draw("colz");
	}
	
	//1D
	TCanvas *c1D[8];
	for(int i = 0; i < 8; i++){
		c1D[i] = new TCanvas(OneDcorr[i]->GetName(),OneDcorr[i]->GetName(),1920,1080);
		gStyle->SetPadLeftMargin(0.1); gStyle->SetPadRightMargin(0.01);
		OneDcorr[i]->SetStats(0);
		OneDcorr[i]->Draw("E1");
	}
	
	
	 
}//End of draw_correlations draw functions

void treehandlerFinal(){
	//correlation_plot(411,411,"DD_correlations.root","D^{+}D^{+}");
	draw_correlations("DD_correlations_MONASH.root");
	
}
