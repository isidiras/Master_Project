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

void plot_all(TString filename){
	//Dictionaries used
	const char* DPhi_Names[6] = {"hDPhiLL","hDPhiIL","hDPhiII","hDPhiHL","hDPhiHI","hDPhiHH"};
    //Histogram arrays
    TH1D* hDPhi[6];

    TH1D* hTriggers = Read_Hist(filename,"hTrPt");
    for(int i = 0; i<6;i++){
        hDPhi[i] = Read_Hist(filename,DPhi_Names[i]);
    }
    //Normalizing
    for (int i = 0; i < 6; i++){
		if(i == 0){
			normalize(hDPhi[i],hTriggers,1.,3.);
        }
        if(i == 1 || i == 2){
			normalize(hDPhi[i],hTriggers,3.,8.);
        }
        if(i == 3 || i == 4 || i == 5){
			normalize(hDPhi[i],hTriggers,8.,50.);
        }

    
    }

    //Plot cosmetics
     for(int i=0;i<6;i++){
        hDPhi[i]->SetMarkerStyle(20);
        hDPhi[i]->SetMarkerSize(1);
        hDPhi[i]->SetMarkerColor(4);
        hDPhi[i]->GetYaxis()->SetLabelSize(0.055);
        hDPhi[i]->GetYaxis()->SetTitleSize(0.055);
        hDPhi[i]->GetYaxis()->SetTitleOffset(2.2);
        hDPhi[i]->GetXaxis()->SetLabelSize(0.055);
        hDPhi[i]->GetXaxis()->SetTitleSize(0.055);
        hDPhi[i]->GetXaxis()->SetTitleOffset(1.0);
        hDPhi[i]->SetStats(0);

     }
     //Drawing
     gStyle->SetTitleW(1.08); gStyle->SetTitleH(0.08);
    TCanvas *c1 = new TCanvas();
    c1->Divide(3,2);
    c1->cd(1);
    gPad->SetTopMargin(0.13);
    gPad->SetLeftMargin(0.28);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.13);
    hDPhi[0]->Draw("PE1");

    c1->cd(2);
    gPad->SetTopMargin(0.13);
    gPad->SetLeftMargin(0.28);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.13);
    hDPhi[1]->Draw("PE1");

    c1->cd(3);
    gPad->SetTopMargin(0.13);
    gPad->SetLeftMargin(0.28);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.13);
    hDPhi[2]->Draw("PE1");

    c1->cd(4);
    gPad->SetTopMargin(0.13);
    gPad->SetLeftMargin(0.28);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.13);
    hDPhi[3]->Draw("PE1");

    c1->cd(5);
    gPad->SetTopMargin(0.13);
    gPad->SetLeftMargin(0.28);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.13);
    hDPhi[4]->Draw("PE1");

    c1->cd(6);
    gPad->SetTopMargin(0.13);
    gPad->SetLeftMargin(0.28);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.13);
    hDPhi[5]->Draw("PE1");

    
}

void draw_correlation_plots(){
    plot_all("ccbar_Junctions_Hard_High/complete_root/DplusDminus.root");
    plot_all("ccbar_Junctions_Hard_High/complete_root/DplusDplus.root");
    plot_all("ccbar_Junctions_Hard_High/complete_root/LplusLminus.root");
    plot_all("ccbar_Junctions_Hard_High/complete_root/LplusLplus.root");
    plot_all("ccbar_Junctions_Hard_High/complete_root/LplusDminus.root");
    plot_all("ccbar_Junctions_Hard_High/complete_root/LplusDplus.root");

}