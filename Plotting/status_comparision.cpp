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

void plot_mothers_ID(TString filename, char* part){
	//This function returns a histogram were the mother's ID
	//of the decay produced particles is shown.

	//Reading Histograms
	TH1D* hist1 = Read_Hist(filename,"hTrMotherL");
	TH1D* hist2 = Read_Hist(filename,"hTrMotherI");
	TH1D* hist3 = Read_Hist(filename,"hTrMotherH");

	//Cosmetics
	hist1->SetLineColor(2);
	hist1->SetFillColor(2);
	hist2->SetLineColor(3);
	hist2->SetFillColor(3);
	hist3->SetLineColor(4);
	hist3->SetFillColor(4);

	//Creating Stack
	THStack* hStack = new THStack("hStack",Form(" ;ID;Counts",part));
	hStack->Add(hist1);
	hStack->Add(hist2);
	hStack->Add(hist3);

	//Creating Legend
	TLegend *leg = new TLegend(0.8,0.8,0.95,0.9);
	leg->SetBorderSize(0);
	leg->AddEntry(hist1,"1 GeV/c #leq #it{p_{T}} < 3 GeV/c","LF");
	leg->AddEntry(hist2,"3 GeV/c #leq #it{p_{T} < 8 GeV/c}","LF");
	leg->AddEntry(hist3,"8 GeV/c #leq #it{p_{T}}");

	//Drawing
	TCanvas *c1 = new TCanvas();
	hStack->Draw("nostackb");
	leg->Draw();
}



void status_comparision(){
    TString filenamesOS[6] = {"ccbar_Monash_Hard_Low/complete_root/LplusLminus.root","ccbar_Monash_Hard_High/complete_root/LplusLminus.root","ccbar_Monash_Soft/complete_root/LplusLminus.root","ccbar_Junctions_Hard_Low/complete_root/LplusLminus.root","ccbar_Junctions_Hard_High/complete_root/LplusLminus.root","ccbar_Junctions_Soft/complete_root/DplusDminus.root"};
    TString filenamesSS[6] = {"ccbar_Monash_Hard_Low/complete_root/LplusLplus.root","ccbar_Monash_Hard_High/complete_root/LplusLplus.root","ccbar_Monash_Soft/complete_root/LplusLplus.root","ccbar_Junctions_Hard_Low/complete_root/LplusLplus.root","ccbar_Junctions_Hard_High/complete_root/LplusLplus.root","ccbar_Junctions_Soft/complete_root/LplusLplus.root"};
    TString multiplicities[6] = {"/home/isidiras/university_staff/ccbar_MONASH_Hard_low/gen_hist.root","/home/isidiras/university_staff/ccbar_MONASH_Hard_High/general_histograms.root","/home/isidiras/university_staff/ccbar_MONASH_Soft/general_histograms.root","/home/isidiras/university_staff/ccbar_Junctions_Hard_Low/general_histograms.root","/home/isidiras/university_staff/ccbar_Junctions_Hard_High/general_histograms.root","/home/isidiras/university_staff/ccbar_Junctions_Soft/general_histograms.root"};
    const char* Leg_entries[6] = {"MONASH HardQCD pTHatmin=1 GeV/c","MONASH HardQCD pTHatmin=10 GeV/c","MONASH SoftQCD","Junctions HardQCD pTHatmin=1 GeV/c","Junctions HardQCD pTHatmin=10 GeV/c","Junctions SoftQCD"};

    TLegend *leg = new TLegend();
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05);

    TH1D *histSS[6];
    TH1D *histOS[6];
    TH1D *hSize[6];
    THStack *hs1 = new THStack();
    THStack *hs2 = new THStack();

    for(int i = 0; i < 6; i++){
        histSS[i] = Read_Hist(filenamesSS[i],"hStatusTr");
        histOS[i] = Read_Hist(filenamesSS[i],"hStatusAs");
        hSize[i] = Read_Hist(multiplicities[i],"hSize");
        histSS[i]->Scale(1/hSize[i]->GetEntries());
        histOS[i]->Scale(1/hSize[i]->GetEntries());
        histSS[i]->SetTitle("");
        histSS[i]->GetYaxis()->SetTitle("#frac{1}{N_{events}}#frac{dN_{tr.}}{dPr_{ID}}");
        histSS[i]->GetYaxis()->SetLabelSize(0.055);
        histSS[i]->GetYaxis()->SetTitleSize(0.055);
        histSS[i]->GetYaxis()->SetTitleOffset(1.3);
        histSS[i]->GetXaxis()->SetLabelSize(0.055);
        histSS[i]->GetXaxis()->SetTitleSize(0.055);
        histSS[i]->GetXaxis()->SetTitleOffset(1.0);
        histSS[i]->GetXaxis()->SetRangeUser(80,95);
        histSS[i]->SetLineWidth(2);
        histSS[i]->SetFillColor(i+1);
        histSS[i]->SetLineColor(i+1);
        histSS[i]->SetStats(0);
        hs1->Add(histSS[i]);
        histOS[i]->SetTitle("");
        histOS[i]->GetYaxis()->SetTitle("#frac{1}{N_{events}}#frac{dN_{as.}}{dPr_{ID}}");
        histOS[i]->GetYaxis()->SetLabelSize(0.055);
        histOS[i]->GetYaxis()->SetTitleSize(0.055);
        histOS[i]->GetYaxis()->SetTitleOffset(1.1);
        histOS[i]->GetXaxis()->SetLabelSize(0.055);
        histOS[i]->GetXaxis()->SetTitleSize(0.055);
        histOS[i]->GetXaxis()->SetTitleOffset(1.0);
        histOS[i]->GetXaxis()->SetRangeUser(80,95);
        histOS[i]->SetLineWidth(2);
        //histOS[i]->SetFillColor(i+1);
        histOS[i]->SetLineColor(i+1);
        histOS[i]->SetStats(0);
        leg->AddEntry(histSS[i],Leg_entries[i],"L");
    }


    
    histSS[4]->GetYaxis()->SetRangeUser(0.000001,0.7);
    TCanvas *c1 = new TCanvas();
    c1->SetLeftMargin(0.17);
    c1->SetRightMargin(0.07);
    c1->SetBottomMargin(0.12);
    c1->SetLogy();
    hs1->Draw("nostackb");
    leg->Draw();

    histOS[4]->GetYaxis()->SetRangeUser(0.000001,0.7);
    TCanvas *c2 = new TCanvas();
    c2->SetLeftMargin(0.17);
    c2->SetRightMargin(0.07);
    c2->SetBottomMargin(0.12);
    c2->SetLogy();
    histOS[4]->Draw("hist");
    histOS[5]->Draw("hist same");
    histOS[3]->Draw("hist same");
    histOS[1]->Draw("hist same");
    histOS[2]->Draw("hist same");
    histOS[0]->Draw("hist same");
    leg->Draw();

plot_mothers_ID("ccbar_Monash_Hard_High/complete_root/LplusLminus.root","sas");

}




