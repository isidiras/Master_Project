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

void draw_same_side_peak(){
     TString filenamesSS[3] = {"ccbar_Junctions_Hard_Low/complete_root/LplusDplus.root","ccbar_Junctions_Hard_High/complete_root/LplusDplus.root","ccbar_Junctions_Soft/complete_root/LplusDplus.root"};
    TH1D *plotsH[3];
    TH1D *plotsW[3];
    TH1D *pT[3];

    for(int i = 0; i < 3; i++){
        plotsH[i] = Read_Hist(filenamesSS[i],"hDPhiIIPr");
        plotsW[i] = Read_Hist(filenamesSS[i],"hDPhiIISc");
        pT[i] = Read_Hist(filenamesSS[i],"hTrPt");
        //Normalize
        normalize(plotsH[i],pT[i],3,8);
        normalize(plotsW[i],pT[i],3,8);
        //Cosmetics
        plotsH[i]->SetMarkerStyle(20);
        plotsH[i]->SetMarkerSize(1.5);
        plotsH[i]->SetMarkerColor(i+3);
        plotsH[i]->SetTitle("");
        plotsH[i]->GetYaxis()->SetLabelSize(0.055);
        plotsH[i]->GetYaxis()->SetTitleSize(0.055);
        plotsH[i]->GetYaxis()->SetTitleOffset(1.6);
        plotsH[i]->GetXaxis()->SetLabelSize(0.055);
        plotsH[i]->GetXaxis()->SetTitleSize(0.055);
        plotsH[i]->GetXaxis()->SetTitleOffset(1.0);
        plotsH[i]->SetStats(0);
        plotsW[i]->SetMarkerStyle(20);
        plotsW[i]->SetMarkerSize(1.5);
        plotsW[i]->SetMarkerColor(i+3);
        plotsW[i]->SetTitle("");
        plotsW[i]->GetYaxis()->SetLabelSize(0.055);
        plotsW[i]->GetYaxis()->SetTitleSize(0.055);
        plotsW[i]->GetYaxis()->SetTitleOffset(1.6);
        plotsW[i]->GetXaxis()->SetLabelSize(0.055);
        plotsW[i]->GetXaxis()->SetTitleSize(0.055);
        plotsW[i]->GetXaxis()->SetTitleOffset(1.0);
        plotsW[i]->SetStats(0);
    }
    
    TLegend *legJ = new TLegend(0.1,0.7,0.3,0.88);
     legJ->SetBorderSize(0);
     legJ->SetFillColor(0);
     legJ->SetTextFont(42);
     legJ->SetTextSize(0.045);
     legJ->SetHeader("Decay Produced","C");
     //legJ->AddEntry(plotsW[0],"JunctionspTHatMin = 1 GeV/c","P");
     legJ->AddEntry(plotsW[1]," Junctions pTHatMin = 10 GeV/c","P");
     //legJ->AddEntry(plotsW[2],"Junctions SoftQCD","P");

     TLegend *legM = new TLegend(0.1,0.7,0.3,0.88);
     legM->SetBorderSize(0);
     legM->SetFillColor(0);
     legM->SetTextFont(42);
     legM->SetTextSize(0.045);
     legM->SetHeader("Hadronisation","C");
     //legM->AddEntry(plotsH[0],"Junctions pTHatMin = 1 GeV/c","P");
     legM->AddEntry(plotsH[1],"Junctions pTHatMin = 10 GeV/c","P");
     //legM->AddEntry(plotsH[2],"Junctions SoftQCD","P");

     TLatex *text = new TLatex();
    text->SetTextFont(42);
    text->SetTextSize(0.07);
    

    TCanvas *c1 = new TCanvas();
    c1->SetLeftMargin(0.2);
    c1->SetRightMargin(0.1);
    c1->SetBottomMargin(0.15);
    plotsH[1]->GetYaxis()->SetRangeUser(0.00009,0.0002);
   // plotsH[0]->Draw("PE1");
    plotsH[1]->Draw("PE1");
   // plotsH[2]->Draw("PE1 SAME");
    legM->Draw();
    text->DrawLatexNDC(0.42,0.8,"3 GeV/c<p_{T}<8 GeV/c");

    TCanvas *c2 = new TCanvas();
    c2->SetLeftMargin(0.2);
    c2->SetRightMargin(0.1);
    c2->SetBottomMargin(0.15);
    plotsW[1]->GetYaxis()->SetRangeUser(0.00001,0.0003);
    //plotsW[0]->Draw("PE1");
    plotsW[1]->Draw("PE1");
    //plotsW[2]->Draw("PE1 SAME");
    legJ->Draw();
    text->DrawLatexNDC(0.42,0.8,"3 GeV/c<p_{T}<8 GeV/c");



}