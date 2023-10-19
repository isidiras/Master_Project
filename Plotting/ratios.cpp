//With this script we are calculating the ratio or
//the subtractions of the correlation histograms
//between the oppsite and same sign and see how this ratio 
//differs pair tune.

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

using namespace std;
using namespace TMath;

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

void ratios(){
    const char* DPhi_Names[6] = {"hDPhiLL","hDPhiIL","hDPhiII","hDPhiHL","hDPhiHI","hDPhiHH"};
	const char* DPhi_NamesPr[6] = {"hDPhiLLPr","hDPhiILPr","hDPhiIIPr","hDPhiHLPr","hDPhiHIPr","hDPhiHHPr"};
    const char* DPhi_NamesB[6] = {"hDPhiLLSc","hDPhiILSc","hDPhiIISc","hDPhiHLSc","hDPhiHISc","hDPhiHHSc"};
    const char* Ratio_Names[6] = {"Sign Ratio Low-Low","Sign Ratio Intermediate-Low","Sign Ratio Intermediate-Intermediate","Sign Ratio High-Low","Sign Ratio High-Intermediate","Sign Ratio High-High"};
    const char* leg_low[3] = {"pTHatmin = 1 GeV/c Total","pTHatmin = 1 GeV/c String Fragmented","pTHatmin = 1 GeV/c Decay"};
    const char* leg_high[3] = {"pTHatmin = 10 GeV/c Total","pTHatmin = 10 GeV/c String Fragmented","pTHatmin = 10 GeV/c Decay"};
    const char* Names[3][6];

    for(int i = 0; i < 6; i++){
        Names[0][i] = DPhi_Names[i];
        Names[1][i] = DPhi_NamesPr[i];
        Names[2][i] = DPhi_NamesB[i];
    }


    TH1D* histSSLow[3][6];
    TH1D* histOSLow[3][6];
    TH1D* histSSHigh[3][6];
    TH1D* histOSHigh[3][6];
    TH1D* ratioLow[3][6];
    TH1D* ratioHigh[3][6];

    TString filenameSSLow = "/home/isidiras/university_staff/ccbar_MONASH_Hard_low/complete_root/LplusDplus.root";
    TString filenameOSLow = "/home/isidiras/university_staff/ccbar_MONASH_Hard_low/complete_root/LplusDminus.root";
    TString filenameSSHigh = "/home/isidiras/university_staff/ccbar_MONASH_Hard_High/complete_root/LplusDplus.root";
    TString filenameOSHigh = "/home/isidiras/university_staff/ccbar_MONASH_Hard_High/complete_root/LplusDminus.root";


    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 6; j++){
            histSSLow[i][j] = Read_Hist(filenameSSLow,Names[i][j]);
            histOSLow[i][j] = Read_Hist(filenameOSLow,Names[i][j]);
            histSSHigh[i][j] = Read_Hist(filenameSSHigh,Names[i][j]);
            histOSHigh[i][j] = Read_Hist(filenameOSHigh,Names[i][j]);
            ratioLow[i][j] = diff(histOSLow[i][j],histSSLow[i][j]);
            ratioHigh[i][j] = diff(histOSHigh[i][j],histSSHigh[i][j]);
            ratioLow[i][j]->GetYaxis()->SetTitle("#frac{O.S.}{S.S.}");
            ratioHigh[i][j]->GetYaxis()->SetTitle("#frac{O.S.}{S.S.}");
            ratioLow[i][j]->SetTitle(Ratio_Names[j]);
            ratioHigh[i][j]->SetTitle(Ratio_Names[j]);
            ratioHigh[i][j]->SetLineColor(2);
            ratioLow[i][j]->GetYaxis()->SetRangeUser( ratioHigh[i][j]->GetMinimum(),ratioLow[i][j]->GetMaximum());
            ratioLow[i][j]->SetStats(0);
            ratioHigh[i][j]->SetStats(0);
        }
    }

    //Plotting
    TCanvas *c1[3];
    TLegend *leg[3];


    for(int i = 0; i < 3; i++){
        c1[i] = new TCanvas();
        c1[i]->Divide(3,2);
        leg[i] = new TLegend();
        leg[i]->AddEntry(ratioLow[i][3],leg_low[i],"l");
        leg[i]->AddEntry(ratioHigh[i][3],leg_high[i],"l");
        for(int j = 0; j < 6; j++){
            c1[i]->cd(j+1);
             ratioLow[i][j]->Draw();
             ratioHigh[i][j]->Draw("SAME");
             if(j == 3) leg[i]->Draw();
        }
    }

}