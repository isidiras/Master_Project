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

void inter_ratios(){
    TH1D *histSS[6];
    TH1D *histOS[6];
    TH1D* ratios[6];

    TString filenamesOS[6] = {"ccbar_Monash_Hard_Low/complete_root/LplusDminus.root","ccbar_Monash_Hard_High/complete_root/LplusDminus.root","ccbar_Monash_Soft/complete_root/LplusDminus.root","ccbar_Junctions_Hard_Low/complete_root/LplusDminus.root","ccbar_Junctions_Hard_High/complete_root/LplusDminus.root","ccbar_Junctions_Soft/complete_root/LplusDminus.root"};
    TString filenamesSS[6] = {"ccbar_Monash_Hard_Low/complete_root/LplusDplus.root","ccbar_Monash_Hard_High/complete_root/LplusDplus.root","ccbar_Monash_Soft/complete_root/LplusDplus.root","ccbar_Junctions_Hard_Low/complete_root/LplusDplus.root","ccbar_Junctions_Hard_High/complete_root/LplusDplus.root","ccbar_Junctions_Soft/complete_root/LplusDplus.root"};
    
    const char* Leg_entries[6] = {"MONASH HardQCD pTHatmin=1 GeV/c","MONASH HardQCD pTHatmin=10 GeV/c","MONASH SoftQCD","Junctions HardQCD pTHatmin=1 GeV/c","Junctions HardQCD pTHatmin=10 GeV/c","Junctions SoftQCD"};
    TLegend *leg = new TLegend();
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05);
    
    for(int i = 0; i < 6; i++){
        histSS[i] = Read_Hist(filenamesSS[i],"hDPhiII");
        histOS[i] = Read_Hist(filenamesOS[i],"hDPhiII");
        ratios[i] = Divide(histOS[i],histSS[i]);
        ratios[i]->SetTitle("");
        ratios[i]->GetYaxis()->SetTitle("#frac{O.S.}{S.S.}");
        ratios[i]->GetYaxis()->SetLabelSize(0.055);
        ratios[i]->GetYaxis()->SetTitleSize(0.055);
        ratios[i]->GetYaxis()->SetTitleOffset(1.1);
        ratios[i]->GetXaxis()->SetLabelSize(0.055);
        ratios[i]->GetXaxis()->SetTitleSize(0.055);
        ratios[i]->GetXaxis()->SetTitleOffset(1.0);
        ratios[i]->SetLineWidth(2);
        ratios[i]->SetLineColor(i+1);
        ratios[i]->SetStats(0);
        leg->AddEntry(ratios[i],Leg_entries[i],"L");

    }

    TLatex *text = new TLatex();
    text->SetTextFont(42);
    text->SetTextSize(0.07);

    ratios[0]->GetYaxis()->SetRangeUser(-0.6,8);
    TCanvas *c1 = new TCanvas();
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.1);
    c1->SetBottomMargin(0.12);
    ratios[0]->Draw("HIST");
    ratios[1]->Draw("HIST SAME");
    ratios[2]->Draw("HIST SAME");
    ratios[3]->Draw("HIST SAME");
    ratios[4]->Draw("HIST SAME");
    ratios[5]->Draw("HIST SAME");
    leg->Draw();
    text->DrawLatexNDC(0.42,0.8,"3 GeV/c < p_{T} < 8 GeV/c");


    
    
    }
