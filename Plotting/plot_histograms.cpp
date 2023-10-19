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

void normalize(TH1D *hist, TH1D *trig, Double_t xmin, Double_t xmax){
	Int_t bin_start = trig->FindBin(xmin);
	Int_t bin_finish = trig->FindBin(xmax);
	Double_t integral = trig->Integral(bin_start, bin_finish);
	hist->Scale(1./integral);
	hist->GetYaxis()->SetTitle("#frac{1}{N_{tr}} #frac{dN}{d#Delta#phi}");
}


void plot_histograms(){
    const char* DPhi_Names[6] = {"hDPhiLL","hDPhiIL","hDPhiII","hDPhiHL","hDPhiHI","hDPhiHH"};
    const char* Triggers_Names[3] = {"hTrPt","hTrPtPr","hTrPtSc"};
    const char* FileNames[3] = {"/home/isidiras/university_staff/ccbar_MONASH_Hard_low/complete_root/DplusDminus.root","/home/isidiras/university_staff/ccbar_MONASH_Hard_High/complete_root/DplusDminus.root","/home/isidiras/university_staff/ccbar_Junctions_Hard_Low/complete_root/DplusDminus.root"};

    TH1D *htrig[3];
    TH1D *hDPhi[3][6];

    //Importing_plot
    for(int i = 0; i < 3; i++){
        htrig[i] = Read_Hist(FileNames[i],"hTrPt");
        for(int j = 0; j < 6; j++){
            hDPhi[i][j] = Read_Hist(FileNames[i],DPhi_Names[j]);
            //Normalizing
            if(j == 0 ) normalize(hDPhi[i][j],htrig[i],1,3);
            if(j == 1 || j == 2) normalize(hDPhi[i][j],htrig[i],3,8);
            if(j == 3 || j == 4 || j == 5) normalize(hDPhi[i][j],htrig[i],8,50);
            hDPhi[i][j]->SetMarkerStyle(20);
            hDPhi[i][j]->SetMarkerColor(i+1);
            hDPhi[i][j]->SetStats(0);

        }
    }

    TLegend *leg = new TLegend();
  	//leg->SetHeader("Correlation Type","c");
  	leg->AddEntry(hDPhi[0][0],"MONASH pTHatmin = 1 Gev/c","p");
  	leg->AddEntry(hDPhi[1][0],"MONASH pTHatmin = 10 GeV/c","p");
  	leg->AddEntry(hDPhi[2][0],"Juctions pTHatmin = 1 Gev/c","p");
    gStyle->SetLegendTextSize(0.045);
  	leg->SetBorderSize(0);

    TCanvas *c1[6];

    for(int i = 0; i < 6; i++){
        c1[i] = new TCanvas();
        c1[i]->Divide(3,1);
        c1[i]->cd(1);
        hDPhi[0][i]->Draw("PE1");
        c1[i]->cd(2);
        hDPhi[1][i]->Draw("PE1");
        c1[i]->cd(3);
        hDPhi[2][i]->Draw("PE1");
        leg->Draw();
    }

}