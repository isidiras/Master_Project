//This script creates ratio plots to
//compare between different outputs


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
#include "TLatex.h"

#define PI 3.14159265
using namespace std;
using namespace TMath;



TH1D *Divide(TH1D *hist1, TH1D *hist2){
	//Returns hist1/hist2
	TH1D *div = (TH1D*)hist1->Clone();
	div->Divide(hist2);
	return div;
}

TH1D *diff(TH1D *hist_1, TH1D *hist_2){
//This function subtracts histograms and returns the subtracted (1 -2)
TH1D *subhist = (TH1D*)hist_1->Clone();
subhist->Add(hist_2,-1);
return subhist;
}


void script_vector_ratio(TString filename1,TString filename2){
//filename1->script
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
	TH1D *DD_script = (TH1D*)input1->Get("hDeltaPhiDD");
	TH1D *DD_vector = (TH1D*)input2->Get("hDeltaPhi");
	TH1D *div = Divide(DD_script,DD_vector);
	div->SetTitle("Script Vector Ratio");
	div->GetYaxis()->SetTitle("#frac{Script}{Vector}");
	div->Draw();
	
	

}

void data_simulations(TString filename1,TString filename2, TString filename3){
//filename1 old simulation

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
	
	TH1D *DD_old = (TH1D*)input1->Get("hDeltaPhiDD");
	TH1D *DDbar_old = (TH1D*)input1->Get("hDeltaPhiDDbar");
	TH1D *DD_new = (TH1D*)input2->Get("hDeltaPhiDD");
	TH1D *DDbar_new = (TH1D*)input3->Get("hDeltaPhi");
	//TH1D *hTrigL = (TH1D*)input2->Get("hTrigL");
	//TH1D *hTrigI = (TH1D*)input2->Get("hTrigI");
	//TH1D *hTrigH = (TH1D*)input2->Get("hTrigH");
	
	
	TH1D *ratio_old = Divide(DDbar_old,DD_old);
	TH1D *ratio_new = Divide(DDbar_new, DD_new);
	
	ratio_old->SetTitle("Ratio from unnormalized distributions from old simulation");
	ratio_old->GetYaxis()->SetTitle("#frac{O.S.}{S.S.}");
	ratio_new->SetTitle("Ratio from unnormalized distributions from new simulation");
	ratio_new->GetYaxis()->SetTitle("#frac{O.S.}{S.S.}");
	
	/*//Normalizing to triggers
	DD_old->Scale(1./ (32.838 * pow(10,6)));
	DDbar_old->Scale(1. /(32.838 * pow(10,6)));
	DD_new->Scale(1. /(2.731*pow(10,8)));
	DDbar_new->Scale(1. /(2.731*pow(10,8)));
	
	DD_old->Scale(1./ DD_old->GetEntries());
	DDbar_old->Scale(1. /DDbar_old->GetEntries());
	DD_new->Scale(1. /DD_new->GetEntries());
	DDbar_new->Scale(1. /DDbar_new->GetEntries());
	*/
	
	TH1D* ratio_SS = Divide(DD_new,DD_old);
	TH1D* ratio_OS = Divide(DDbar_new, DDbar_old);
	
	ratio_SS->SetTitle("Same Sign Ratio of old (script) produced new (vector) produced");
	ratio_SS->GetYaxis()->SetTitle("#frac{new}{old}");
	
	ratio_OS->SetTitle("Opposite Sign Ratio of old (script) produced new (vector) produced");
	ratio_OS->GetYaxis()->SetTitle("#frac{new}{old}");
	
	
	TH1D* sub_ss = diff(DD_new,DD_old);
	TH1D* sub_os = diff(DDbar_new,DDbar_old);
	
	sub_ss->SetTitle("Same Sigh Subtraction");
	sub_ss->GetYaxis()->SetTitle("new-old");
	
	sub_os->SetTitle("Opposite Sigh Subtraction");
	sub_os->GetYaxis()->SetTitle("new-old");
	
	
	TCanvas* c1 = new TCanvas();
	c1->Divide(2,1);
	c1->cd(1);
	ratio_old->Draw("E1");
	c1->cd(2);
	ratio_new->Draw();
	
	TCanvas* c2 = new TCanvas();
	c2->Divide(2,1);
	c2->cd(1);
	ratio_SS->Draw("E1");
	c2->cd(2);
	ratio_OS->Draw();
	
	TCanvas* c3 = new TCanvas();
	c3->Divide(2,1);
	c3->cd(1);
	sub_ss->Draw("E1");
	c3->cd(2);
	sub_os->Draw();
	
	
	
		
}
void sign_ratio(TString filename1, TString filename2){
//This function gives the ratio of the correlation histograms
//filename1->same sign
//filename2->opposite sign
	//Histogram arrays
	TH1D *DPhi_Hist1[6];
	TH1D *DPhi_Hist2[6];
	TH1D *DPhi_Hist3[6];
	TH1D *hTriggers1[3];
	TH1D *hTriggers2[3];
	//Dictionaries
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
		}
	}
	
	
		
		for(int i = 0; i<6; i++){
			//Normalizing
		if(i == 0){
			DPhi_Hist1[i]->Scale(1. /hTriggers1[0]->Integral(1,10));
			DPhi_Hist2[i]->Scale(1. /hTriggers2[0]->Integral(1,10));
			//DPhi_Hist3[i]->Scale(1. /hTriggers3[0]->Integral(2,10));
			
			}
		if(i == 1 || i == 2){
			DPhi_Hist1[i]->Scale(1. /hTriggers1[1]->Integral(1,10));
			DPhi_Hist2[i]->Scale(1. /hTriggers2[1]->Integral(1,10));
			//DPhi_Hist3[i]->Scale(1. /hTriggers3[1]->Integral(2,10));
		}
		if(i == 3 || i == 4 || i == 5){
			DPhi_Hist1[i]->Scale(1. /hTriggers1[2]->Integral(1,10));
			DPhi_Hist2[i]->Scale(1. /hTriggers2[2]->Integral(1,10));
			//DPhi_Hist3[i]->Scale(1. /hTriggers3[2]->Integral(2,10));
		}
		
		
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

	


void ratio_simulations(){

	data_simulations("/home/isidiras/university_staff/hf291122/DD_old_MONASH.root","/home/isidiras/university_staff/hf291122/DD_MONASH.root","DDbar_correlation_MONASH.root");
	//sign_ratio("DD_correlation_MONASH.root","DDbar_correlation_MONASH.root");
	//script_vector_ratio("general_hist_MONASH.root","DD_corr_MONASH.root");
}
