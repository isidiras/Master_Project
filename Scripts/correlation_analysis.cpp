//This script contains functions that can be used to create various plots for correlation analysis
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
	
	int ntrees = 3; //Number of trees we want to add to the TChain
	
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
	Double_t trigL,trigI,trigH;//Trigger particle counter variable
	
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
	TH3D *hPtrPaDPhi = new TH3D("hDPtrPaDPhi",Form("#Delta#phi p_{T} Triger p_{T} Associate for %s;p_{T} triger (Gev/c); p_{T} associate (GeV/c); #Delta#phi (rad)",title),100,0,50,100,0,50,100,-PI/2,3*PI/2);
	TH2D *hDeltaPhiPt = new TH2D("hDeltaPhiPt",Form("%s #Delta#phi and trigger p_{T};p_{T} (GeV/c); #Delta#phi (rad)",title),100,0,50,100,-PI/2,3*PI/2);
	TH2D *hDEtaPt = new TH2D("hDEtaPt",Form("%s #Delta#eta and trigger p_{T};p_{T} (GeV/c);#Delta#eta",title),100,0,50,80,-8,8);
	TH2D *hDPhiDEta = new TH2D("hDPhiDEta",Form("%s #Delta#phi and #Delta#eta;#Delta#phi (rad);#Delta#eta",title),100,-PI/2,3*PI/2,80,-8,8);
	TH2D *htrIDasID = new TH2D("htrIDasID",Form("%s Triger and Associate ID; trigger ID;associate ID",title),12001,-6000.5,6000.5,12001,-6000.5,6000.5);
	//DeltaPhi correlation plots
	TH1D* hDeltaPhi = new TH1D("hDeltaPhi",Form("%s correlations;#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiLL = new TH1D("hDPhiLL",Form("%s #Delta#phi correlation for low-low p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiIL = new TH1D("hDPhiIL",Form("%s #Delta#phi correlation for intermediate-low p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiII = new TH1D("hDPhiII",Form("%s #Delta#phi correlation for intermediate-intermediate p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiHL = new TH1D("hDPhiHL",Form("%s #Delta#phi correlation for high-low p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiHI = new TH1D("hDPhiHI",Form("%s #Delta#phi correlation for high-intermediate p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiHH = new TH1D("hDPhiHH",Form("%s #Delta#phi correlation for high-high p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	//DeltaEta correlation plots
	TH1D* hDEta = new TH1D("hDEta",Form("%s #Delta#eta correlation;#Delta#eta;Counts",title),80,-8,8);
	TH1D* hDEtaLL = new TH1D("hDEtaLL",Form("%s #Delta#eta correlation for low-low p_{T};#Delta#eta;Counts",title),80,-8,8);
	TH1D* hDEtaIL = new TH1D("hDEtaIL",Form("%s #Delta#eta correlation for intermediate-low p_{T};#Delta#eta;Counts",title),80,-8,8);
	TH1D* hDEtaII = new TH1D("hDEtaII",Form("%s #Delta#eta correlation for intermediate-intermediate p_{T};#Delta#eta;Counts",title),80,-8,8);
	TH1D* hDEtaHL = new TH1D("hDEtaHL",Form("%s #Delta#eta correlation for high-low p_{T};#Delta#eta;Counts",title),80,-8,8);
	TH1D* hDEtaHI = new TH1D("hDEtaHI",Form("%s #Delta#eta correlation for high-intermediate p_{T};#Delta#eta;Counts",title),80,-8,8);
	TH1D* hDEtaHH = new TH1D("hDEtaHH",Form("%s #Delta#eta correlation for high-high p_{T};#Delta#eta;Counts",title),80,-8,8);
	//Trigger counter plots
	TH1D* hTrigL = new TH1D("hTrigL",Form("%s Triggers for low p_{T};#frac{Triggers}{Event};Counts",title),10,0.5,10.5);
	TH1D* hTrigI = new TH1D("hTrigI",Form("%s Triggers for intermediate p_{T};#frac{Triggers}{Event};Counts",title),10,0.5,10.5);
	TH1D* hTrigH = new TH1D("hTrigH",Form("%s Triggers for high p_{T};#frac{Triggers}{Event};Counts",title),10,0.5,10.5);
	
	//Begining of event loop
	for(int iEvent = 0; iEvent < nEvents; iEvent++){
		ch1->GetEntry(iEvent);
		int nparticles = vID->size();
		//Initializing trigger counters
		trigL = 0;
		trigI = 0;
		trigH = 0;

		//Starting particle loop
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
			//Trigger loop
			if(pid == id_trigger){
				nTrigger++;
				if(pPt >= 1 && pPt < 3) trigL++;
				if(pPt >= 3 && pPt < 8) trigI++;
				if(pPt >= 8) trigH++;
				for(int jparticle = 0; jparticle < nparticles; jparticle++){
					if (iparticle == jparticle) continue;//So we do not correlate with itself.
						aid = (*vID)[jparticle];
						aPhi = (*vPhi)[jparticle];
						aEta = (*vEta)[jparticle];
						aPt = (*vPt)[jparticle];
						if(aid == id_associate){
							//Filling general correlationplots
							hPtrPaDPhi->Fill(pPt,aPt,DeltaPhi(pPhi,aPhi));
							hDeltaPhi->Fill(DeltaPhi(pPhi,aPhi));
							hDeltaPhiPt->Fill(pPt,DeltaPhi(pPhi,aPhi));
							hDPhiDEta->Fill(DeltaPhi(pPhi,aPhi),pEta-aEta);
							hDEtaPt->Fill(pPt,pEta-aEta);
							hDEta->Fill(pEta-aEta);
							//Filling triger momentum range correlations
							if(pPt >= 1 && pPt< 3 && aPt < 3 && aPt<pPt){
								hDPhiLL->Fill(DeltaPhi(pPhi,aPhi));
								hDEtaLL->Fill(pEta-aEta);
								htrIDasID->Fill(pid,aid);
								
							}//End of low transverse momentum
							if(pPt >= 3 && pPt < 8 && aPt < 3){
								hDPhiIL->Fill(DeltaPhi(pPhi,aPhi));
								hDEtaIL->Fill(pEta-aEta);
								htrIDasID->Fill(pid,aid);
																
							}//End of intermidiate-low transverse momentum
							
							if(pPt >= 3 && pPt < 8 && aPt >=3 && aPt < 8 && aPt<pPt){
								hDPhiII->Fill(DeltaPhi(pPhi,aPhi));
								hDEtaII->Fill(pEta-aEta);
								htrIDasID->Fill(pid,aid);
								
							}
							
							if(pPt >= 8 && aPt < 3){
								hDPhiHL->Fill(DeltaPhi(pPhi,aPhi));
								hDEtaHL->Fill(pEta-aEta);
								htrIDasID->Fill(pid,aid);
								
							}//End of high-low transverse momentum
							
							if(pPt >= 8 && aPt >=3 && aPt < 8){
								hDPhiHI->Fill(DeltaPhi(pPhi,aPhi));
								hDEtaHI->Fill(pEta-aEta);
								htrIDasID->Fill(pid,aid);
								
							}//End of high-intermediate transverse momentum
							
							if(pPt >= 8 && aPt >= 8 && aPt<pPt){
								hDPhiHH->Fill(DeltaPhi(pPhi,aPhi));
								hDEtaHH->Fill(pEta-aEta);
								htrIDasID->Fill(pid,aid);
							
							}//End of high-high transverse momentum
					}//End of DD correlation
				}//End of associate particle loop D+D+
			}//D+ triger
			
		}//End of  particle loop
		hTrigL->Fill(trigL);
		hTrigI->Fill(trigI);
		hTrigH->Fill(trigH);	
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

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

void draw_correlations(TString filename){
	//This function will take the root file with histograms we created using correlation_plot and draw the histograms
	TFile *input = TFile::Open(filename);
	if(!input->IsOpen() || !input){
		cout<<"File: "<<filename<<" is not found!"<<endl;
		cout<<"Terminating!"<<endl;
		return;
	}
	
	//Creating histogram arrays that I am going to use
	TH1D *general_hist[5];
	TH1D *DPhi_Hist[7];
	TH1D *DEta_Hist[7];
	TH1D *hTriggers[3];
	TH2D *TwoDcorr[3];
	
	//Arrays of histogram names
	const char* TwoDnames[3] = {"hDeltaPhiPt","hDEtaPt","hDPhiDEta"};
	const char* DPhi_Names[7] = {"hDeltaPhi","hDPhiLL","hDPhiIL","hDPhiII","hDPhiHL","hDPhiHI","hDPhiHH"};
	const char* DEta_Names[7] = {"hDEta","hDEtaLL","hDEtaIL","hDEtaII","hDEtaHL","hDEtaHI","hDEtaHH"};
	const char* General_Names[5] = {"hPt","hEta","hPhi","hID","hCharge"};
	const char* Triggers_Names[3] = {"hTrigL","hTrigI","hTrigH"};
	
	
	//Here we are importing the histograms
	for(int i = 0; i<5; i++){
		general_hist[i] = (TH1D*)input->Get(General_Names[i]);
		if(!general_hist[i]){
			cout<<General_Names[i]<<" not found in "<<filename<<endl;
			cout<<"Terminating!"<<endl;
			return;
		}
	}
	//Importing Delta Phi
	for(int i = 0; i<7; i++){
		DPhi_Hist[i] = (TH1D*)input->Get(DPhi_Names[i]);
		if(!DPhi_Hist[i]){
			cout<<DPhi_Names[i]<<" not found in "<<filename<<endl;
			cout<<"Terminating"<<endl;
		}
	}
	
	//Importing Delta Eta
	for(int i = 0; i<7; i++){
		DEta_Hist[i] = (TH1D*)input->Get(DEta_Names[i]);
		if(!DEta_Hist[i]){
			cout<<DEta_Names[i]<<" not found in "<<filename<<endl;
			cout<<"Terminating"<<endl;
		}
	}
	
	//Importing 2D histograms
	for(int i = 0; i<3; i++){
		TwoDcorr[i] = (TH2D*)input->Get(TwoDnames[i]);
		if(!TwoDcorr[i]){
			cout<<TwoDnames[i]<<" not found in "<<filename<<endl;
			cout<<"Terminating!"<<endl;
			return;
		}
	}
	//Importing trigger histogram
	for(int i = 0; i<3; i++){
		hTriggers[i] = (TH1D*)input->Get(Triggers_Names[i]);
		if(!hTriggers[i]){
			cout<<Triggers_Names[i]<<" not found in "<<filename<<endl;
			cout<<"Terminating!"<<endl;
		}
	}
	
	//Importing 3D histogram
	TH3D* PtrPaDPhi = (TH3D*)input->Get("hDPtrPaDPhi");
	if(!PtrPaDPhi){
		cout<<"3D histogram not found in "<<filename<<endl;
		cout<<"Terminating"<<endl;
		return;
	}
	
	//Plot looks
	for(int i = 0; i<7;i++){
		DEta_Hist[i]->SetMarkerStyle(7);
		DPhi_Hist[i]->SetMarkerStyle(8);
		DPhi_Hist[i]->SetLineWidth(0);
		if(i<5){
			general_hist[i]->SetMarkerStyle(7);
		}
		
		if(i>0){
		DPhi_Hist[i]->SetMarkerColor(i);
		}
	}
	
	//Normalozing histograms to the number of events
	
	/*for(int i = 0; i<7; i++){
		if(i == 0){
			DEta_Hist[i]->Scale(1. /(hTriggers[0]->Integral()+hTriggers[1]->Integral()+hTriggers[2]->Integral()));
			DPhi_Hist[i]->Scale(1. /(hTriggers[0]->Integral()+hTriggers[1]->Integral()+hTriggers[2]->Integral()));
		}
		if(i == 1){
			DEta_Hist[i]->Scale(1. /hTriggers[0]->Integral());
			DPhi_Hist[i]->Scale(1. /hTriggers[0]->Integral());
		}
		if(i == 2 || i == 3){
			DEta_Hist[i]->Scale(1. /hTriggers[1]->Integral());
			DPhi_Hist[i]->Scale(1. /hTriggers[1]->Integral());
		}
		if(i == 4 || i == 5 || i == 6){
			DEta_Hist[i]->Scale(1. /hTriggers[2]->Integral());
			DPhi_Hist[i]->Scale(1. /hTriggers[2]->Integral());
			
		}
		DEta_Hist[i]->GetYaxis()->SetTitle("#frac{1}{N_{tr}} #frac{dN_{assoc}}{d#Delta#eta}");
		DPhi_Hist[i]->GetYaxis()->SetTitle("#frac{1}{N_{tr}} #frac{dN_{assoc}}{d#Delta#phi}");
	}
	*/
	//Normalizing 2D histograms
	const char*  TwoD_AxisNames[3] = {"#frac{1}{N_{tr}} #frac{d^{2}N_{assoc}}{d#Delta#phi dp_{T}}","#frac{1}{N_{tr}} #frac{d^{2}N_{assoc}}{d#Delta#eta dp_{T}}","#frac{1}{N_{tr}} #frac{d^{2}N_{assoc}}{d#Delta#phi d#Delta#eta}"};
	for(int i = 0; i < 3; i++){
		TwoDcorr[i]->Scale(1./(hTriggers[0]->Integral()+hTriggers[1]->Integral()+hTriggers[2]->Integral()));
		TwoDcorr[i]->GetZaxis()->SetTitle(TwoD_AxisNames[i]);
	}
	
	TCanvas *c1 = new TCanvas();
	TLegend *leg = new TLegend(0.8,0.7,0.9,0.9);
	for(int i = 1; i<7; i++){
	leg->AddEntry(DPhi_Hist[i],DPhi_Names[i],"p");
	}
	DPhi_Hist[1]->SetStats(0);
	DPhi_Hist[1]->Draw("P E1");
	for (int i = 2; i < 7;i++){
		DPhi_Hist[i]->Draw("P E1 same");
	}
	leg->Draw();
	
	
	//Drawing in canvases
	//3D
	/*TCanvas *c3D = new TCanvas();
	PtrPaDPhi->SetStats(0);
	PtrPaDPhi->Draw("BOX2");
	*/
	/*
	//2D
	TCanvas *c2D[3];
	gStyle->SetPalette(kRainBow);
	for(int i = 0; i < 3; i++){
		c2D[i] = new TCanvas(TwoDcorr[i]->GetName(),TwoDcorr[i]->GetName(),1920,1080);
		//gStyle->SetPadLeftMargin(0.08); gStyle->SetPadRightMargin(0.2);
		//TwoDcorr[i]->GetZaxis()->SetTitle("Counts");
		TwoDcorr[i]->SetContour(80);
		TwoDcorr[i]->SetStats(0);
		TwoDcorr[i]->Draw("SURF2");
	}
	
	//1D
	
	TCanvas *c_general[5];
	for(int i = 0; i < 5; i++){
		c_general[i] = new TCanvas(general_hist[i]->GetName(),general_hist[i]->GetName(),1920,1080);
		gStyle->SetPadLeftMargin(0.1); gStyle->SetPadRightMargin(0.01);
		//general_hist[i]->SetStats(0);
		general_hist[i]->Draw("E1");
	}
	*/
	/*
	
	TCanvas *c_DPhi[7];
	for(int i = 0; i < 7; i++){
		gStyle->SetPadLeftMargin(0.15); gStyle->SetPadRightMargin(0.05);
		c_DPhi[i] = new TCanvas(DPhi_Hist[i]->GetName(),DPhi_Hist[i]->GetName(),1920,1080);
		DPhi_Hist[i]->SetStats(0);
		DPhi_Hist[i]->Draw("E1");
		c_DPhi[i]->Print(Form("%s.png",DPhi_Hist[i]->GetName()));
		
	}
	
	TCanvas *c_Eta[7];
	for(int i = 0; i < 7; i++){
		gStyle->SetPadLeftMargin(0.15); gStyle->SetPadRightMargin(0.05);
		c_Eta[i] = new TCanvas(DEta_Hist[i]->GetName(),DEta_Hist[i]->GetName(),1920,1080);
		DEta_Hist[i]->SetStats(0);
		DEta_Hist[i]->Draw("E1");
		c_Eta[i]->Print(Form("%s.png",DEta_Hist[i]->GetName()));
		}
	
	TCanvas *c_Triggers[3];
	for(int i = 0; i<3; i++){
		c_Triggers[i] = new TCanvas(hTriggers[i]->GetName(),hTriggers[i]->GetName(),1920,1080);
		gStyle->SetPadLeftMargin(0.1); gStyle->SetPadRightMargin(0.01);
		//hTriggers[i]->SetStats(0);
		hTriggers[i]->Draw("E1");
	}
	*/
	
	
	
}//End of draw_correlations draw functions
//+-++-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+---+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+

TGraphErrors *triger_PT(TString filename){
//This function takes the filename as input reads the file takes the trigger and momentum histogram
//and returns a graph of the number of triggers as a function of p_T

	//Opening file
	TFile *input = TFile::Open(filename);
	if(!input->IsOpen() || !input){
		cout<<"File: "<<filename<<" is not found!"<<endl;
		cout<<"Terminating!"<<endl;
		return 0;
	}
	//TFile *input2 = TFile::Open(general_hist_MONASH.root);

	//Allocation of graph we will return
	TGraphErrors *gr = new TGraphErrors();
	
	//Trigger histograms
	TH1D *hTriggers[3];
	const char* Triggers_Names[3] = {"hTrigL","hTrigI","hTrigH"};
	//Momentum histogram
	TH1D *hPT;
	//Importing trigger histogram
	for(int i = 0; i<3; i++){
		hTriggers[i] = (TH1D*)input->Get(Triggers_Names[i]);
		if(!hTriggers[i]){
			cout<<Triggers_Names[i]<<" not found in "<<filename<<endl;
			cout<<"Terminating!"<<endl;
			return 0;
		}
	}
	//Importing Pt distribution
	hPT = (TH1D*)input->Get("hPt");
	if(!hPT){
		cout<<"hPt not found in "<<filename<<endl;
		cout<<"Terminating!"<<endl;
		return 0;
	}
	
	//Input for x axis
	Double_t x[3];
	Double_t xerr[3];
	
	hPT->GetXaxis()->SetRangeUser(1.,3.);
	x[0] = hPT->GetMean();
	xerr[0] = hPT->GetStdDev();
	
	hPT->GetXaxis()->SetRangeUser(3.,8.);
	x[1] = hPT->GetMean();
	xerr[1] = hPT->GetStdDev();
	
	hPT->GetXaxis()->SetRangeUser(8.,50.);
	x[2] = hPT->GetMean();
	xerr[2] = hPT->GetStdDev();
	
	//Input for y axis
	Double_t y[3];
	Double_t yerr[3];
	for(int i = 0; i<3; i++){
		y[i] = hTriggers[i]->IntegralAndError(1,hTriggers[i]->GetNbinsX(),yerr[i],"");
		
	}
	//Checking if filled properly
	for(int i = 0; i<3; i++){
		cout<<x[i]<<" "<<xerr[i]<<endl;
		cout<<y[i]<<" "<<yerr[i]<<endl;
		cout<<"Next Point"<<endl;
	}
	
	//Filling graph
	
	Int_t n = 0;
	for(int i = 0; i<3;i++){
		n = gr->GetN();
		gr->SetPoint(n,x[i],y[i]);
		gr->SetPointError(n,xerr[i],yerr[i]);
	}
	
	gr->SetMarkerStyle(7);
	gr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gr->GetYaxis()->SetTitle("Count");
	gr->SetTitle("Number of Triggers as function of p_{T}");
	return gr;

}//End of triger PTfunction

//+-+-+-+-+-+-+-+-+-+-+-++-+-+-+-+-+-+-+-+

void plot_all_combinations(TString filename1, TString filename2, TString filename3){
//This functions takes the three filenames as input reads them takes the correlation and trigger
//histograms normalizes to the number of triggers and plot the correlations of each kind for a
//specific momentum range
//always put same sign correlations in filename2 and filename3
	//Allocations of arrays
	TH1D *DPhi_Hist1[6];
	TH1D *DPhi_Hist2[6];
	TH1D *DPhi_Hist3[6];
	TH1D *hTriggers1[3];
	TH1D *hTriggers2[3];
	TH1D *hTriggers3[3];
	TLegend *leg[6];
	TCanvas *c1[6];
	
	//Dictionaries used
	const char* DPhi_Names[6] = {"hDPhiLL","hDPhiIL","hDPhiII","hDPhiHL","hDPhiHI","hDPhiHH"};
	const char* Triggers_Names[3] = {"hTrigL","hTrigI","hTrigH"};
	const char* titles[6] = {"Low-Low p_{T} correlations","Intermediate-Low p_{T} correlations","Intermediate-Intermediate  p_{T} correlations","High-Low p_{T} correlations","High-Intermediate p_{T} correlations","High-High p_{T} correlations"};
	
	//Reading the files
	TFile *input1 = TFile::Open(filename1);
	if(!input1->IsOpen() || !input1){
		cout<<"File: "<<filename1<<" is not found!"<<endl;
		cout<<"Terminating!"<<endl;
		return 0;
	}
	
	TFile *input2 = TFile::Open(filename2);
	if(!input2->IsOpen() || !input2){
		cout<<"File: "<<filename2<<" is not found!"<<endl;
		cout<<"Terminating!"<<endl;
		return 0;
	}
	
	TFile *input3 = TFile::Open(filename3);
	if(!input3->IsOpen() || !input3){
		cout<<"File: "<<filename3<<" is not found!"<<endl;
		cout<<"Terminating!"<<endl;
		return 0;
	}
	
	//Importing histograms
	for(int i = 0; i<6; i++){
		DPhi_Hist1[i] = (TH1D*)input1->Get(DPhi_Names[i]);
		if(!DPhi_Hist1[i]){
			cout<<DPhi_Names[i]<<" not found in "<<filename1<<endl;
			cout<<"Terminating"<<endl;
		}
		DPhi_Hist2[i] = (TH1D*)input2->Get(DPhi_Names[i]);
		if(!DPhi_Hist2[i]){
			cout<<DPhi_Names[i]<<" not found in "<<filename2<<endl;
			cout<<"Terminating"<<endl;
		}
		DPhi_Hist3[i] = (TH1D*)input3->Get(DPhi_Names[i]);
		if(!DPhi_Hist3[i]){
			cout<<DPhi_Names[i]<<" not found in "<<filename3<<endl;
			cout<<"Terminating"<<endl;
		}
		if(i<3){
			hTriggers1[i] = (TH1D*)input1->Get(Triggers_Names[i]);
			if(!hTriggers1[i]){
				cout<<Triggers_Names[i]<<" not found in "<<filename1<<endl;
				cout<<"Terminating!"<<endl;
				return 0;
			}
			hTriggers2[i] = (TH1D*)input2->Get(Triggers_Names[i]);
			if(!hTriggers2[i]){
				cout<<Triggers_Names[i]<<" not found in "<<filename2<<endl;
				cout<<"Terminating!"<<endl;
				return 0;
			}
			hTriggers3[i] = (TH1D*)input3->Get(Triggers_Names[i]);
			if(!hTriggers3[i]){
				cout<<Triggers_Names[i]<<" not found in "<<filename3<<endl;
				cout<<"Terminating!"<<endl;
				return 0;
			}			
		
		}//End of Importing triggers
	}
	
	//Normalizing
	for(int i = 0; i<6; i++){
                
		if(i == 0){
			DPhi_Hist1[i]->Scale(1. /hTriggers1[0]->Integral(1,10));
			DPhi_Hist2[i]->Scale(1. /hTriggers2[0]->Integral(2,10));
			DPhi_Hist3[i]->Scale(1. /hTriggers3[0]->Integral(2,10));
			
		}
		if(i == 1 || i == 2){
			DPhi_Hist1[i]->Scale(1. /DPhi_Hist1[i]->GetEntries());
			DPhi_Hist2[i]->Scale(1. /DPhi_Hist2[i]->GetEntries());
			DPhi_Hist3[i]->Scale(1. /hTriggers3[1]->Integral(2,10));
		}
		if(i == 3 || i == 4 || i == 5){
			DPhi_Hist1[i]->Scale(1. /hTriggers1[2]->Integral());
			DPhi_Hist2[i]->Scale(1. /hTriggers2[2]->Integral(2,10));
			DPhi_Hist3[i]->Scale(1. /hTriggers3[2]->Integral(2,10));
		}
		
		DPhi_Hist1[i]->GetYaxis()->SetTitle("#frac{1}{N_{tr}} #frac{dN_{assoc}}{d#Delta#phi}");
		DPhi_Hist2[i]->GetYaxis()->SetTitle("#frac{1}{N_{tr}} #frac{dN_{assoc}}{d#Delta#phi}");
		DPhi_Hist3[i]->GetYaxis()->SetTitle("#frac{1}{N_{tr}} #frac{dN_{assoc}}{d#Delta#phi}");
		
		DPhi_Hist1[i]->SetMarkerStyle(7);
		DPhi_Hist2[i]->SetMarkerStyle(7);
		DPhi_Hist3[i]->SetMarkerStyle(7);
		DPhi_Hist1[i]->SetMarkerColor(2);
		DPhi_Hist2[i]->SetMarkerColor(3);
		DPhi_Hist3[i]->SetMarkerColor(4);
		DPhi_Hist1[i]->SetTitle(titles[i]);
		DPhi_Hist2[i]->SetTitle(titles[i]);
		DPhi_Hist1[i]->SetStats(0);
		DPhi_Hist2[i]->SetStats(0);
		DPhi_Hist3[i]->SetStats(0);
	}
	
	//Legends
	for(int i = 0; i<6; i++){
		leg[i]=new TLegend(0.89,0.7,0.99,0.9);
		leg[i]->AddEntry(DPhi_Hist1[i],"D^{+}D^{-}","P");
		leg[i]->AddEntry(DPhi_Hist2[i],"D^{+}D^{+}","P");
		leg[i]->AddEntry(DPhi_Hist3[i],"D^{-}D^{-}","P");
		leg[i]->SetBorderSize(0);
	}
	TCanvas *c2 = new TCanvas();
	c2->Divide(3,2);
	//Plotting
	for(int i = 0; i<6; i++){
		gStyle->SetPadLeftMargin(0.15); gStyle->SetPadRightMargin(0.01);
		//c1[i] = new TCanvas(DPhi_Hist1[i]->GetName(),DPhi_Hist1[i]->GetName(),1920,1080);
		//c1[i]->Divide(2,1);
		//c1[i]->cd(1);
		//DPhi_Hist2[i]->GetYaxis()->SetRangeUser(0.00001,0.0025);
		//DPhi_Hist2[i]->Draw("PE1");
		//DPhi_Hist3[i]->Draw("PE1 same");
		//DPhi_Hist1[i]->Draw("PE1 same");
		//c1[i]->cd(2);
		//DPhi_Hist2[i]->Draw("PE1");
		//DPhi_Hist3[i]->Draw("PE1 same");
		//leg[i]->Draw();
		c2->cd(i+1);
		DPhi_Hist2[i]->Draw("PE1");
		DPhi_Hist3[i]->Draw("PE1 same");
		DPhi_Hist1[i]->Draw("PE1 same");
		leg[i]->Draw();
	}
	

}//End of plot_all_combinations functions

void sign_ratio(TString filename1, TString filename2){
//This function gives the ratio of the correlation histograms
//filename1->same sign
//filename2->opposite sign
	//Histogram arrays
	TH1D *DPhi_Hist1[6];
	TH1D *DPhi_Hist2[6];
	TH1D *DPhi_Hist3[6];
	//Dictionaries
	const char* DPhi_Names[6] = {"hDPhiLL","hDPhiIL","hDPhiII","hDPhiHL","hDPhiHI","hDPhiHH"};
	const char* titles[6] = {"Low-Low p_{T} correlations","Intermediate-Low p_{T} correlations","Intermediate-Intermediate  p_{T} correlations","High-Low p_{T} correlations","High-Intermediate p_{T} correlations","High-High p_{T} correlations"};
	
	//Reading the files
	TFile *input1 = TFile::Open(filename1);
	if(!input1->IsOpen() || !input1){
		cout<<"File: "<<filename1<<" is not found!"<<endl;
		cout<<"Terminating!"<<endl;
		return 0;
	}
	
	TFile *input2 = TFile::Open(filename2);
	if(!input2->IsOpen() || !input2){
		cout<<"File: "<<filename2<<" is not found!"<<endl;
		cout<<"Terminating!"<<endl;
		return 0;
	}
	
	for(int i = 0; i<6; i++){
		DPhi_Hist1[i] = (TH1D*)input1->Get(DPhi_Names[i]);
		if(!DPhi_Hist1[i]){
			cout<<DPhi_Names[i]<<" not found in "<<filename1<<endl;
			cout<<"Terminating"<<endl;
		}
		DPhi_Hist2[i] = (TH1D*)input2->Get(DPhi_Names[i]);
		if(!DPhi_Hist2[i]){
			cout<<DPhi_Names[i]<<" not found in "<<filename2<<endl;
			cout<<"Terminating"<<endl;
		}
	}
		
		for(int i = 0; i<6; i++){
			DPhi_Hist3[i] = (TH1D*)DPhi_Hist2[i]->Clone();
			DPhi_Hist3[i]->Divide(DPhi_Hist1[i]);
			DPhi_Hist3[i]->SetTitle(titles[i]);
			DPhi_Hist3[i]->GetYaxis()->SetTitle("#frac{opposite sign}{same sign}");
			DPhi_Hist3[i]->SetMarkerStyle(7);
			DPhi_Hist3[i]->SetMarkerColor(2);
		}
		
		TCanvas *c1 = new TCanvas();
		c1->Divide(3,2);
		for(int i = 0; i<6; i++){
			c1->cd(i+1);
			gStyle->SetPadLeftMargin(0.15); gStyle->SetPadRightMargin(0.01);
			DPhi_Hist3[i]->Draw("PE1");
		}
	
}


void correlation_analysis(){
	//correlation_plot(-411,-411,"DbarDbar_correlations.root","D^{-}D^{-}");
	//correlation_plot(-411,+411,"DbarD_correlations.root","D^{-}D^{+}");
	//draw_correlations("DDbar_corr_MONAH.root");
	
	//TGraphErrors *gr = triger_PT("DD_corr_MONASH.root");
	//TGraphErrors *gr1 = triger_PT("DbarDbar_corr_MONAH.root");
	/*
	gr->SetMarkerColor(2);
	gr->SetLineColor(2);
	gr1->SetMarkerColor(4);
	gr1->SetLineColor(4);
	gr1->SetMarkerStyle(6);
	
	TLegend *leg = new TLegend(0.5,0.6,0.8,0.8);
	leg->SetBorderSize(0);
	leg->AddEntry(gr,"D^{+} trigger","LP");
	leg->AddEntry(gr1,"D^{-} trigger","LP");
	
	TCanvas *c1 = new TCanvas();
	gr->Draw("ALP");
	gr1->Draw("SAME");
	leg->Draw();
	*/
	plot_all_combinations("DDbar_corr_MONAH.root","DD_corr_MONASH.root","DbarDbar_corr_MONAH.root");
	//sign_ratio("DD_corr_MONASH.root","DDbar_corr_MONAH.root");
}
