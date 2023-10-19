//With this small script I am reading the multiplicity plots and charm per event plot
//that are given directly from the simulation and I compare them by plotting the together


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
    input->Close();
    
}

void multiplicity_comparisson(){
    //Importing histograms
    TH1D *multiplicityLow = Read_Hist("/home/isidiras/university_staff/ccbar_MONASH_Hard_low/gen_hist.root","hSize");
    TH1D *multiplicityHigh = Read_Hist("/home/isidiras/university_staff/ccbar_MONASH_Hard_High/general_histograms.root","hSize");
    TH1D *multiplicitySoft = Read_Hist("/home/isidiras/university_staff/ccbar_MONASH_Soft/general_histograms.root","hSize");
    TH1D *multiplicityJunctionsLow = Read_Hist("/home/isidiras/university_staff/ccbar_Junctions_Hard_Low/general_histograms.root","hSize");
    TH1D *multiplicityJuctionsHigh = Read_Hist("/home/isidiras/university_staff/ccbar_Junctions_Hard_High/general_histograms.root","hSize");
    TH1D *multiplicityJunctionsSoft = Read_Hist("/home/isidiras/university_staff/ccbar_Junctions_Soft/general_histograms.root","hSize");
    TH1D *charmLow = Read_Hist("/home/isidiras/university_staff/ccbar_MONASH_Hard_low/gen_hist.root","hCharmPart");
    TH1D *charmHigh = Read_Hist("/home/isidiras/university_staff/ccbar_MONASH_Hard_High/general_histograms.root","hCharmPart");
    TH1D *charmSoft = Read_Hist("/home/isidiras/university_staff/ccbar_MONASH_Soft/general_histograms.root","hCharmPart");
    TH1D *charmJunctionsLow = Read_Hist("/home/isidiras/university_staff/ccbar_Junctions_Hard_Low/general_histograms.root","hCharmPart");
    TH1D *charmJuctionsHigh = Read_Hist("/home/isidiras/university_staff/ccbar_Junctions_Hard_High/general_histograms.root","hCharmPart");
    TH1D *charmJuctionsSoft = Read_Hist("/home/isidiras/university_staff/ccbar_Junctions_Soft/general_histograms.root","hCharmPart");

    //Normalizing to the number of events
    multiplicityLow->Scale(1/multiplicityLow->GetEntries());
    multiplicityHigh->Scale(1/multiplicityHigh->GetEntries());
    charmLow->Scale(1/charmLow->GetEntries());
    charmHigh->Scale(1/charmHigh->GetEntries());
    charmSoft->Scale(1/charmSoft->GetEntries());
    multiplicityJunctionsLow->Scale(1/multiplicityJunctionsLow->GetEntries());
    multiplicityJuctionsHigh->Scale(1/multiplicityJuctionsHigh->GetEntries());
    multiplicityJunctionsSoft->Scale(1/multiplicityJunctionsSoft->GetEntries());
    multiplicitySoft->Scale(1/multiplicitySoft->GetEntries());
    charmJunctionsLow->Scale(1/charmJunctionsLow->GetEntries());
    charmJuctionsHigh->Scale(1/charmJuctionsHigh->GetEntries());
    charmJuctionsSoft->Scale(1/charmJuctionsSoft->GetEntries());
    multiplicityLow->SetStats(0);
    multiplicityHigh->SetStats(0);
    charmLow->SetStats(0);
    charmHigh->SetStats(0);
    multiplicityJunctionsLow->SetStats(0);
    charmJunctionsLow->SetStats(0);
    charmJuctionsHigh->SetStats(0);
    multiplicityJuctionsHigh->SetStats(0);
    multiplicitySoft->SetStats(0);
    charmSoft->SetStats(0);
    multiplicityJunctionsSoft->SetStats(0);
    charmJuctionsSoft->SetStats(0);


    //Axis-Titles
     multiplicityLow->GetYaxis()->SetTitle("P(X=x)");
     multiplicityHigh->GetYaxis()->SetTitle("P(X=x)");
     charmLow->GetYaxis()->SetTitle("P(X=x)");
     charmHigh->GetYaxis()->SetTitle("P(X=x)");
     multiplicityJunctionsLow->GetYaxis()->SetTitle("P(X=x)");
     charmJunctionsLow->GetYaxis()->SetTitle("P(X=x)");
     multiplicityJuctionsHigh->GetYaxis()->SetTitle("P(X=x)");
     charmJuctionsHigh->GetYaxis()->SetTitle("P(X=x)");
     multiplicitySoft->GetYaxis()->SetTitle("P(X=x)");
     charmSoft->GetYaxis()->SetTitle("P(X=x)");
     multiplicityLow->GetYaxis()->SetLabelSize(0.055);
     multiplicityLow->GetYaxis()->SetTitleSize(0.055);
     multiplicityLow->GetYaxis()->SetTitleOffset(1.0);
    multiplicityLow->GetXaxis()->SetLabelSize(0.055);
    multiplicityLow->GetXaxis()->SetTitleSize(0.055);
     multiplicityLow->GetXaxis()->SetTitleOffset(1.0);
     multiplicityLow->GetXaxis()->SetTitle("#frac{Particles}{Event}");
     charmLow->GetYaxis()->SetLabelSize(0.055);
    charmLow->GetYaxis()->SetTitleSize(0.055);
    charmLow->GetYaxis()->SetTitleOffset(1.0);
    charmLow->GetXaxis()->SetLabelSize(0.055);
    charmLow->GetXaxis()->SetTitleSize(0.055);
    charmLow->GetXaxis()->SetTitleOffset(1.0);
     charmLow->GetXaxis()->SetTitle("#frac{Charm Particles}{Event}");
     //Line color
     multiplicityHigh->SetLineColor(2);
     charmHigh->SetLineColor(2);
     multiplicitySoft->SetLineColor(1);
     charmSoft->SetLineColor(1);
     multiplicityJunctionsLow->SetLineStyle(2);
     charmJunctionsLow->SetLineStyle(2);
     multiplicityJuctionsHigh->SetLineStyle(2);
     charmJuctionsHigh->SetLineStyle(2);
     multiplicityJuctionsHigh->SetLineColor(2);
     charmJuctionsHigh->SetLineColor(2);
     multiplicityJunctionsSoft->SetLineStyle(2);
     multiplicityJunctionsSoft->SetLineColor(1);
    charmJuctionsSoft->SetLineColor(1);
    charmJuctionsSoft->SetLineStyle(2);

     //Legends
     TLegend *legM = new TLegend(0.8,0.8,0.99,0.9);
     TLegend *legC = new TLegend(0.8,0.8,0.99,0.9);

     legM->AddEntry(multiplicityLow,"MONASH pTHatmin = 1 Gev/c","L");
     legM->AddEntry(multiplicityHigh,"MONASH pTHatmin = 10 GeV/c","L");
     legM->AddEntry(multiplicitySoft,"MONASH SoftQCD","L");
     legM->AddEntry( multiplicityJunctionsLow,"Juctions pTHatmin = 1 Gev/c","L");
     legM->AddEntry( multiplicityJuctionsHigh,"Juctions pTHatmin = 10 Gev/c","L");
     legM->AddEntry(multiplicityJunctionsSoft,"Juctions SoftQCD","L");
     legM->SetBorderSize(0);
     legC->AddEntry(charmLow,"MONASH pTHatmin = 1 Gev/c","L");
     legC->AddEntry(charmHigh,"MONASH pTHatmin = 10 Gev/c","L");
     legC->AddEntry(charmSoft,"MONASH SoftQCD","L");
     legC->AddEntry(charmJunctionsLow,"Juctions pTHatmin = 1 Gev/c","L");
     legC->AddEntry(charmJuctionsHigh,"Juctions pTHatmin = 10 Gev/c","L");
     legC->AddEntry(charmJuctionsSoft,"Juctions Soft","L");
     legC->SetBorderSize(0);

     TCanvas *c1 = new TCanvas();
     multiplicityLow->Draw("hist");
     multiplicityHigh->Draw("hist SAME");
     multiplicitySoft->Draw("hist SAME");
     multiplicityJunctionsLow->Draw("hist SAME");
     multiplicityJuctionsHigh->Draw("hist SAME");
     multiplicityJunctionsSoft->Draw("hist SAME");
     legM->Draw();

     TCanvas *c2 = new TCanvas();
     charmLow->Draw("hist");
     charmHigh->Draw("hist SAME");
     charmSoft->Draw("hist SAME");
     charmJunctionsLow->Draw("hist SAME");
     charmJuctionsHigh->Draw("hist SAME");
     charmJuctionsSoft->Draw("hist SAME");
     legC->Draw();




     



}