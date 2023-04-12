//This script does the background reduction using the ZYAM method the functions are 
//given with the order of use. Furthermore, with this script the near and away side can be calculated.
//This can happen by changing the arguments in the yield function.


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
#include "TAttMarker.h"

#define PI 3.14159265
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
void trigger_normalization(TH1D *DPhi,TH1D *hTriggers){
	DPhi->Scale(1. / hTriggers->Integral(2,10));
	DPhi->GetYaxis()->SetTitle("#frac{1}{N_{tr}} #frac{dN_{assoc}}{d#Delta#phi}");
	DPhi->SetMarkerStyle(7);
	cout<<"The number of trigers is: "<<hTriggers->Integral(2,10)<<endl;
	
}

double* Zyam_and_Error(Double_t start, Double_t finish, TH1D *hist){
	//This function take the range and histogram as input and returns the Zyam
	//background and the error of it.
	int bin = hist->FindBin(start);
	int n = hist->FindBin(finish);
	double sum = 0;
	double sum_err = 0;
	for(int i = bin; i <= n; i++){
		
		sum = sum + hist->GetBinContent(i);
		sum_err = sum_err+Power(hist->GetBinError(i),2);
	}
	double mean = sum/(n-bin+1);//we add the p+1 because we use the <=
	double *mean_err = new double[2];
	*mean_err = mean;
	*(mean_err+1) = Sqrt(sum_err)/(n-bin);
	cout<<"The ZYAM background is: "<<mean_err[0]<<" +- "<<mean_err[1]<<" for "<<hist->GetName()<<endl;
	
	
	return mean_err;
}

TLine *background_line(Double_t ZYAM,TH1D *hist){

	Double_t x_1 = hist->GetBinCenter(1);
	Double_t x_2 = hist->GetBinCenter(hist->GetNbinsX());
	
	TLine *l = new TLine(x_1,ZYAM, x_2, ZYAM);
	l->SetLineColor(2);
	l->SetLineStyle(9);
	l->SetLineWidth(2);
	return l;
		
}

TH1D *signal(TH1D *hDPhi, Double_t ZYAM, Double_t Error){
//This function gives the signal
	TH1D *hist = (TH1D*)hDPhi->Clone();
	int bins = hist->GetNbinsX();
	for(int i = 1; i <= bins; i++){
		hist->SetBinError(i,Sqrt(Power(hist->GetBinError(i),2)+Power(Error,2)));
		hist->SetBinContent(i,hist->GetBinContent(i)-ZYAM);	
	}
	
	return hist;

}
TH1D *reduction(TH1D *hist,/*TString filename,TString hist_name,TString trig_name*/ Double_t start, Double_t finish){
//This function does the whole procces
	//TH1D *hist = Read_Hist(filename,hist_name);
	//TH1D *htrig = Read_Hist(filename,trig_name);
	//trigger_normalization(hist,htrig);
	double *p = Zyam_and_Error(start,finish,hist);
	TH1D *Signal = signal(hist,*p,*(p+1));
	cout<<"The ZYAM background is: "<<*p<<endl;

	delete p;
	return Signal; 
}

void automated_plot(TString filename, Double_t start[6], Double_t finish[6]){
	//This function plots everything across all momentum ranges
	//Dictionaries
	const char* DPhi_Names[6] = {"hDPhiLL","hDPhiIL","hDPhiII","hDPhiHL","hDPhiHI","hDPhiHH"};
	const char* Triggers_Names[3] = {"hTrigL","hTrigI","hTrigH"};
	
	//Allocation of arrays
	TH1D* hDPhi[6];
	TH1D* hTrig[3];
	double ZYAM[6];
	double ZYAM_err[6];
	TLine* line[6];
	TH1D* reduced[6];
	
	//Importing histograms
	for(int i = 0; i < 6; i++){
		hDPhi[i] = Read_Hist(filename, DPhi_Names[i]);
		if(i < 3){
			hTrig[i] = Read_Hist(filename,Triggers_Names[i]);
		}
	}
	//Nomralizing to the number of triggers
	for(int i = 0; i < 6; i++){
		if(i == 0) trigger_normalization(hDPhi[i],hTrig[0]);
		if(i == 1 || i == 2) trigger_normalization(hDPhi[i],hTrig[1]);
		if(i == 3 || i == 4 || i == 5) trigger_normalization(hDPhi[i],hTrig[2]);
	}
	//Calculation ZYAM and error and line creation
	for(int i = 0; i<6; i++){
		double *pZYAM = Zyam_and_Error(start[i],finish[i],hDPhi[i]);
		line[i] = background_line(*(pZYAM),hDPhi[i]);
		reduced[i] = reduction(hDPhi[i],start[i],finish[i]);
		delete pZYAM;
	}
	//Drawing everything
	TLegend *leg = new TLegend(0.89,0.7,0.99,0.9);
	leg->AddEntry(hDPhi[0], "D^{+}D{+} Correlation","P");
	leg->AddEntry(line[0],"ZYAM background","L");
	
	gStyle->SetPadLeftMargin(0.18); gStyle->SetPadRightMargin(0.01);
	TCanvas *c1 = new TCanvas();
	c1->Divide(3,2);
	TCanvas *c2 = new TCanvas();
	c2->Divide(3,2);
	
	for(int i = 0; i < 6;i++){
		c1->cd(i+1);
		if(i == 0) leg->Draw();
		hDPhi[i]->Draw("P");
		line[i]->Draw();
		c2->cd(i+1);
		reduced[i]->Draw();
	}
}
//untill here same as background_reduction.cpp
TGraphErrors *Yield_Low(TString filename,Double_t start[3], Double_t finish[3],Double_t low_limit, Double_t high_limit){
	//This functions takes the filename the start and finish for the background calculatiion and the limits for the yield calculation
	//and returns a graph where the per trigger yield is displayed when associating with low momentum particles.
	//Dictionaries used
	const char* DPhi_Names[6] = {"hDPhiLL","hDPhiIL","hDPhiHL"};
	const char* Triggers_Names[3] = {"hTrigL","hTrigI","hTrigH"};
	
	//Arrays used
	TH1D* hDPhi[3];
	TH1D* hTrig[3];
	Double_t x[3];
	Double_t y[3];
	Double_t xerr[3];
	Double_t yerr[3];
	TH1D* reduced[3];
	
	//Reading files
	for(int i = 0; i < 3; i++){
		hDPhi[i] = Read_Hist(filename, DPhi_Names[i]);
		hTrig[i] = Read_Hist(filename,Triggers_Names[i]);
	}
	
	TH1D *hTrigPt = Read_Hist(filename,"hPtTrig");
	
	//Normalizing and background reducing
	for(int i = 0; i < 3; i++){
		trigger_normalization(hDPhi[i],hTrig[i]);
		reduced[i] = reduction(hDPhi[i],start[i],finish[i]);
	}
	
	//X-axis input
	hTrigPt->GetXaxis()->SetRangeUser(1.,3.);
	x[0] = hTrigPt->GetMean();
	xerr[0] = hTrigPt->GetStdDev();
	
	hTrigPt->GetXaxis()->SetRangeUser(3.,8.);
	x[1] = hTrigPt->GetMean();
	xerr[1] = hTrigPt->GetStdDev();
	
	hTrigPt->GetXaxis()->SetRangeUser(8.,50.);
	x[2] = hTrigPt->GetMean();
	xerr[2] = hTrigPt->GetStdDev();
	
	//Input for y axis
	int low_bin = reduced[0]->FindBin(low_limit);
	int high_bin = reduced[0]->FindBin(high_limit);
	for(int i = 0; i<3; i++){
		y[i] = reduced[i]->IntegralAndError(low_bin,high_limit,yerr[i],"");
	}
	
	//Setting up and filling graph
	TGraphErrors *gr = new TGraphErrors();
	Int_t n = 0;
	cout<<"For low associate momentum we have:"<<endl;
	for(int i = 0; i<3; i++){
		cout<<"Filling x = "<<x[i]<<" +- "<<xerr[i]<<" and y = "<<y[i]<<" +- "<<yerr[i]<<endl;
		n = gr->GetN();
		gr->SetPoint(n,x[i],y[i]);
		gr->SetPointError(n,xerr[i],yerr[i]);
	}
	gr->SetMarkerStyle(7);
	gr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gr->GetYaxis()->SetTitle("Yield");
	return gr;
}


TGraphErrors *Yield_Intermediate(TString filename,Double_t start[2], Double_t finish[2],Double_t low_limit, Double_t high_limit){
	//This functions takes the filename the start and finish for the background calculatiion and the limits for the yield calculation
	//and returns a graph where the per trigger yield is displayed when associating with intermediate momentum particles.
	//Dictionaries used
	const char* DPhi_Names[2] = {"hDPhiII","hDPhiHI"};
	const char* Triggers_Names[2] = {"hTrigI","hTrigH"};
	
	//Arrays used
	TH1D* hDPhi[2];
	TH1D* hTrig[2];
	Double_t x[2];
	Double_t y[2];
	Double_t xerr[2];
	Double_t yerr[2];
	TH1D* reduced[2];
	
	//Reading files
	for(int i = 0; i < 2; i++){
		hDPhi[i] = Read_Hist(filename, DPhi_Names[i]);
		hTrig[i] = Read_Hist(filename,Triggers_Names[i]);
	}
	
	TH1D *hTrigPt = Read_Hist(filename,"hPtTrig");
	
	//Normalizing and background reducing
	for(int i = 0; i < 2; i++){
		trigger_normalization(hDPhi[i],hTrig[i]);
		reduced[i] = reduction(hDPhi[i],start[i],finish[i]);
	}
	
	//X-axis input
	
	
	hTrigPt->GetXaxis()->SetRangeUser(3.,8.);
	x[0] = hTrigPt->GetMean();
	xerr[0] = hTrigPt->GetStdDev();
	
	hTrigPt->GetXaxis()->SetRangeUser(8.,50.);
	x[1] = hTrigPt->GetMean();
	xerr[1] = hTrigPt->GetStdDev();
	
	//Input for y axis
	int low_bin = reduced[0]->FindBin(low_limit);
	int high_bin = reduced[0]->FindBin(high_limit);
	for(int i = 0; i<2; i++){
		y[i] = reduced[i]->IntegralAndError(low_bin,high_limit,yerr[i],"");
	}
	
	//Setting up and filling graph
	TGraphErrors *gr = new TGraphErrors();
	Int_t n = 0;
	cout<<"For intermediate associate momentum we have:"<<endl;
	for(int i = 0; i<2; i++){
		cout<<"Filling x = "<<x[i]<<" +- "<<xerr[i]<<" and y = "<<y[i]<<" +- "<<yerr[i]<<endl;
		n = gr->GetN();
		gr->SetPoint(n,x[i],y[i]);
		gr->SetPointError(n,xerr[i],yerr[i]);
	}
	gr->SetMarkerStyle(8);
	gr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gr->GetYaxis()->SetTitle("Yield");
	return gr;
}

TGraphErrors *Yield_High(TString filename,Double_t start, Double_t finish,Double_t low_limit, Double_t high_limit){
	//This functions takes the filename the start and finish for the background calculatiion and the limits for the yield calculation
	//and returns a graph where the per trigger yield is displayed when associating with intermediate momentum particles.
	
	

	//Reading files
	
	TH1D *hDPhi = Read_Hist(filename, "hDPhiHH");
	TH1D *hTrig = Read_Hist(filename,"hTrigH");
	
	TH1D *hTrigPt = Read_Hist(filename,"hPtTrig");
	
	//Normalizing and background reducing
	
	trigger_normalization(hDPhi,hTrig);
	TH1D *reduced = reduction(hDPhi,start,finish);
	
	
	//X-axis input
	
	
	Double_t x,xerr,y,yerr;
	
	hTrigPt->GetXaxis()->SetRangeUser(8.,50.);
	x = hTrigPt->GetMean();
	xerr = hTrigPt->GetStdDev();
	
	//Input for y axis
	int low_bin = reduced->FindBin(low_limit);
	int high_bin = reduced->FindBin(high_limit);
	y = reduced->IntegralAndError(low_bin,high_limit,yerr,"");
	
	
	//Setting up and filling graph
	TGraphErrors *gr = new TGraphErrors();
	Int_t n = 0;
	cout<<"For high associate momentum we have:"<<endl;

	cout<<"Filling x = "<<x<<" +- "<<xerr<<" and y = "<<y<<" +- "<<yerr<<endl;
	n = gr->GetN();
	gr->SetPoint(n,x,y);
	gr->SetPointError(n,xerr,yerr);
	
	gr->SetMarkerStyle(33);
	gr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gr->GetYaxis()->SetTitle("Yield");
	return gr;
}
void yield_calculation(){
//TH1D *signal = reduction("DD_correlation_MONASH.root","hDPhiHH","hTrigH",-1,1);

	Double_t start[6] = {-0.4,-0.3,-0.1,-0.8,-1,-1};
	Double_t finish[6] = {0.4,0.4,0.1,0.8,1,1};

	Double_t start_l[3] = {-0.4,-0.3,-0.8};
	Double_t finish_l[3] = {0.4,0.4,0.8};

	Double_t start_i[2] = {-0.1,-1.};
	Double_t finish_i[2] = {0.1,1.};
	
	Double_t start_h = -1;
	Double_t finish_h = 1;


	automated_plot("DD_correlations_Final.root",start,finish);
	TGraphErrors *gr = Yield_Low("DD_correlations_Final.root",start_l,finish_l,-PI/2,PI/2);
	TGraphErrors *gr_i = Yield_Intermediate("DD_correlations_Final.root",start_i,finish_i,-PI/2,PI/2);
	TGraphErrors *gr_h = Yield_High("DD_correlations_Final.root",start_h,finish_h,-PI/2,PI/2);
	TLegend *leg = new TLegend(0.89,0.7,0.99,0.9);
	leg->AddEntry(gr, "1-3 GeV/c","P");
	leg->AddEntry(gr_i,"3-8 GeV/c","P");
	leg->AddEntry(gr_h,"8- GeV/c","P");
	gr->SetTitle("Same Sign (D^{+}D^{+}) Near Side Yields for different associate particle momentum");
	TCanvas *c1 = new TCanvas();
	gr->Draw("AP");
	gr_i->Draw("SAME P");
	gr_h->Draw("SAME");
	leg->Draw();
}
