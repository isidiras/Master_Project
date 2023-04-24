//With this script I am doing a specific analysis concerning the azimuthial correlation plots of D^{+-} mesons from 
//production mechanism. This is layed out so it can be used for other particles as well by changing the IDs. The analysis
//is done for the higher momentum regimes were the correlation curve produced from MONASH is not smooth. But one can change the
//momentum constrains.
//C++ libraries
#include <iostream>
#include <vector>
//ROOT libraries
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TString.h"
#include "THStack.h"
#include "TAttMarker.h"
#include "TLegend.h"

#define PI 3.14159265
using namespace std;
using namespace TMath;

Double_t DeltaPhi(Double_t phi1, Double_t phi2){
	//returns delta phi in range -pi/2 3pi/2
	return fmod(phi1-phi2+2.5*PI,2*PI)-0.5*PI;
	}
	
void status_file(Int_t id_trigger,Int_t id_associate, TString filename, const char* title){
	//This functions takes the trigger and associate id and creates a .root file
	// with the name filename. The title is a string so the particles will be seen
	//in histogram titles.
	
	//Define the TChain
	TChain *ch1 = new TChain("tree");
	TFile *output = new TFile(filename,"RECREATE");
	
	int ntrees = 4; //Number of trees we want to add to the TChain This can be changed by the user.
	
	//We put the trees to the chain
	for( int i = 1; i < ntrees+1;  i++){
	//One can change the file path accordingly to his/hers machine set up
	ch1->Add(Form("/home/isidiras/university_staff/ccbar_MONASH/output_MONASH_STATUS/Group%i/output.root",i));//File is for my local set up!
	}
	
	//Now we define vectors that carry the information at event level.
	vector<Int_t>* vID = 0;
	vector<Double_t>* vPt = 0;
	vector<Double_t>* vPhi = 0;
	vector<Double_t>* vStatus = 0;
	//Setting up chain branch addresses to the vectors defined above
	ch1->SetBranchAddress("ID",&vID);
	ch1->SetBranchAddress("PT",&vPt);
	ch1->SetBranchAddress("PHI",&vPhi);
	ch1->SetBranchAddress("STATUS",&vStatus);
	
	//Definition of variables I am going to use
	Int_t aID,pID;
	Double_t pPt,pPhi,pStatus;//For triger
	Double_t aPt,aPhi,aStatus;//For associate
	int nTrigger = 0;
	
	//These two can be changed from the user according to the status histogram
	Double_t primary_status = 83;
	Double_t secondary_status = 91;
	
	//Each vector is an event number of events analyzed is the total number of vectors
	int nEvents = ch1->GetEntries();
	
	cout<<"The number of events for this analysis is: "<<nEvents<<endl;
	
	//Definition of produced histograms
	TH1D *hTrPt = new TH1D("hTrPt",Form("Trigger Transverse Momentum for %s;p_{T} GeV/c;Counts",title),100,0,50);
	TH1D *hAsPt = new TH1D("hAsPt",Form("Associate Transverse Momentum for %s;p_{T} GeV/c;Counts",title),100,0,50);
	TH1D *hStatus = new TH1D("hStatus",Form("Production mechanism of %s pair; Process ID; Counts",title),365,-182.5,182.5);
	TH1D *hDPhi = new TH1D("hDPhi",Form("#Delta#phi for all processes for %s pair;#Delta#Phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D *hDPhiPr = new TH1D("hDPhiPr",Form("#Delta#phi for production mechanism %f for %s pair;#Delta#Phi (rad);Counts",primary_status,title),100,-PI/2,3*PI/2);
	TH1D *hDPhiSe = new TH1D("hDphiSe",Form("#Delta#phi for production mechanism %f for %s pair;#Delta#Phi (rad);Counts",secondary_status,title),100,-PI/2,3*PI/2);
	
	
	
	//Event Loop
	for(int iEvent = 0; iEvent < nEvents; iEvent++){
		ch1->GetEntry(iEvent);
		int nparticles = vID->size();
		for(int ipart = 0; ipart < nparticles; ipart++){
			pID = (*vID)[ipart];
			pPhi = (*vPhi)[ipart];
			pPt = (*vPt)[ipart];
			pStatus = (*vStatus)[ipart];
			if(pID == id_trigger && pPt >= 8.){//Different trigger momentum regimes change here
				nTrigger++;
				hTrPt->Fill(pPt);
				hStatus->Fill(pStatus);
				if(pStatus == secondary_status) hAsPt->Fill(pPt);
				for(int jpart = 0; jpart < nparticles; jpart++){
					if(jpart == ipart) continue;//Do not correlate with it self
						aID = (*vID)[jpart];
						aPhi = (*vPhi)[jpart];
						aPt = (*vPt)[jpart];
						aStatus = (*vStatus)[jpart];
						if(aID == id_associate && aPt >= 8){//Different associate momentum change here
							hDPhi->Fill(DeltaPhi(pPhi,aPhi));
							if(pStatus == primary_status && aStatus == primary_status){
								hDPhiPr->Fill(DeltaPhi(pPhi,aPhi));
							}//Primary Status Codition
							if(pStatus == secondary_status && aStatus == secondary_status){
								hDPhiSe->Fill(DeltaPhi(pPhi,aPhi));
							}//Secondary Status codition
						}//Associate Codition
				}//Associate Loop
			}//Trigger Codition
		}//Trigger Loop
	}//End of event loop
	
	if(nTrigger == 0){
		cout<<"Have not found any trigger particle with id: "<<id_trigger<<endl;
		output->Close();
		return 0;
	}
	output->Write();
	output->Close();
	cout<<"The total number of triggers is: "<<nTrigger<<endl;
	
	cout<<"File: "<<filename<<" has been created!"<<endl;
}

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

void plot(TString filename){
	//This function reads the file and then plots the histogram inside of it.
	//Dictionaries used
	const char* Hist_names[6] = {"hTrPt","hAsPt","hStatus","hDPhi","hDPhiPr","hDphiSe"};
	TH1D *hist[6];
	TCanvas *c1[6];
	//Reading files
	for(int i = 0; i < 6; i++){
		hist[i] = Read_Hist(filename,Hist_names[i]);
	}
	//Plotting
	for(int i = 0; i < 6; i++){
		c1[i] = new TCanvas();
		hist[i]->SetMarkerStyle(7);
		hist[i]->Draw("P E1");
	}
	
}

void stacking_analysis(TString filename){
	//This function takes as input the filename and returns a plot
	//where the DeltaPhi distribution from all thr production mechanisms
	//is shown along with the distributions from primary and secondary mechanisms
	//and theit stack.
	//Dictionaries used
	const char* Hist_names[3] = {"hDPhi","hDPhiPr","hDphiSe"};
	TH1D *hist[3];
	//Importing
	for(int i = 0; i < 3; i++){
		hist[i] = Read_Hist(filename,Hist_names[i]);
	}
	//Adding histograms
	TH1D *sum = (TH1D*)hist[1]->Clone();
	sum->Add(hist[2]);
	
	//Cosmetics
	for(int i = 0; i<3; i++){
	hist[i]->SetMarkerStyle(7);
	hist[i]->SetMarkerColor(i+1);
	hist[i]->SetStats(0);
	}
	sum->SetMarkerStyle(7);
	sum->SetMarkerColor(4);
	sum->GetYaxis()->SetRangeUser(-2.,120.);
	
	//Legend
	TLegend *leg = new TLegend();
	leg->SetBorderSize(0);
	leg->AddEntry(hist[0],"From all production mechanisms","p");
	leg->AddEntry(hist[1],"From primary production mechanism","p");
	leg->AddEntry(hist[2],"From secondary production mechanism","p");
	leg->AddEntry(sum,"From the combination of primary and secondary","p");
	
	//Drawing
	TCanvas *c1 = new TCanvas();
	sum->Draw("P E1");
	sum->SetStats(0);
	for(int i = 0; i < 3; i++){
		hist[i]->Draw("P E1 SAME");
	}
	leg->Draw();
}
void status_analysis(){
	status_file(411,411,"DD_correlation_status.root","D^{+}D^{+}");
	plot("DD_correlation_status.root");
	stacking_analysis("DD_correlation_status_86.root");
	
	return 0;
}

