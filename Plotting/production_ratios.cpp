//This script produces 
//C++ libraries
#include <iostream>
#include <vector>
#include <fstream>
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

Double_t Covariance(int n ,Double_t x[],Double_t y[]){
    Double_t mean_1 = 0;
    Double_t mean_2 = 0;
    Double_t sum_1 = 0;
    Double_t sum_2 = 0;
    Double_t sum_3 = 0;
    Double_t cov[n];
    for(int i = 0; i < n; i++){
        sum_1 += x[i];
        sum_2 += y[i];
    }
    mean_1 = sum_1/n;
    mean_2 = sum_2/n;

    for(int i = 0; i < n; i++){
        cov[i] = (x[i]-mean_1)*(y[i]-mean_2); 
        sum_3 += cov[i];
    }
    

    return sum_3/n;
}

double_t give_mean(int n,double_t sample[]){
    double_t sum = 0.;
    for(int i = 0; i < n;i++){
        sum += sample[i];
    }


    return sum/(double_t)n;
}

double_t give_std(int n,double_t sample[]){
    double_t sum = 0.;
    double_t mean = give_mean(n,sample);
    for(int i = 0; i<n; i++){
        sum += pow(sample[i]-mean,2);
    }

    return sqrt(sum/((double_t)n-1));
}

double_t rho_test(int n,double_t x[], double_t y[]){
    return Covariance(n,x,y)/(give_std(n,x)*give_std(n,y));
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

TGraphErrors* plot_ratios(TString filename,const char* filepath,int i){
	//This function returns a plot where the ratio of string fragmentation pairs over the total pairs is shown 
	//as a function of p_T
	
	//Dictionaries used
	const char* DPhi_Names[6] = {"hDPhiLL","hDPhiIL","hDPhiII","hDPhiHL","hDPhiHI","hDPhiHH"};
	const char* DPhi_NamesPr[6] = {"hDPhiLLPr","hDPhiILPr","hDPhiIIPr","hDPhiHLPr","hDPhiHIPr","hDPhiHHPr"};
	//Arrays used
	Double_t x[3];
	Double_t entries_total[3];
	Double_t entries_string[3];
    Double_t entries_total_err[6];
    Double_t entries_string_err[6];
    Double_t y_err[0];
    Double_t y_string_err[3];
	TH1D* hDPhi[6];
	TH1D* hDPhiPr[6];
	TH1D* hTrigPt;
	
	//File creation
	fstream file;
	file.open(Form("%s/data%i.txt",filepath,i),ios::out);
	//Reading histograms
	for(int i = 0;i<6;i++){
		hDPhi[i] = Read_Hist(filename,DPhi_Names[i]);
		hDPhiPr[i] = Read_Hist(filename,DPhi_NamesPr[i]);
	}
	hTrigPt = Read_Hist(filename,"hTrPt");
	
	//x-axis
	hTrigPt->GetXaxis()->SetRangeUser(1.,3.);
	Double_t p_assoc_l = hTrigPt->GetMean();
	
	hTrigPt->GetXaxis()->SetRangeUser(1.,3.);
	x[0] = hTrigPt->GetMean();
	
	
	hTrigPt->GetXaxis()->SetRangeUser(3.,8.);
	x[1] = hTrigPt->GetMean();
	
	
	hTrigPt->GetXaxis()->SetRangeUser(8.,50.);
	x[2] = hTrigPt->GetMean();
	
	//y-axis
	entries_total[0] = hDPhi[0]->IntegralAndError(1,100,entries_total_err[0],"");
	entries_string[0] = hDPhiPr[0]->IntegralAndError(1,100,entries_string_err[0],"");
    y_err[0] = entries_total_err[0];
	y_string_err[0] = entries_string_err[0];

	entries_total[1] = hDPhi[1]->IntegralAndError(1,100,entries_total_err[1],"")+hDPhi[2]->IntegralAndError(1,100,entries_total_err[2],"");
	entries_string[1] = hDPhiPr[1]->IntegralAndError(1,100,entries_string_err[1],"")+hDPhiPr[1]->IntegralAndError(1,100,entries_string_err[2],"");
    y_err[1] = Sqrt(Power(entries_total_err[1],2)+Power(entries_total_err[2],2));
    y_string_err[1] =  Sqrt(Power(entries_string_err[1],2)+Power(entries_string_err[2],2));
	
	entries_total[2] = hDPhi[3]->IntegralAndError(1,100,entries_total_err[3],"")+hDPhi[4]->IntegralAndError(1,100,entries_total_err[4],"")+hDPhi[5]->IntegralAndError(1,100,entries_total_err[5],"");
	entries_string[2] = hDPhiPr[3]->IntegralAndError(1,100,entries_string_err[3],"")+hDPhiPr[4]->IntegralAndError(1,100,entries_string_err[4],"")+hDPhiPr[5]->IntegralAndError(1,100,entries_string_err[5],"");
	y_err[2] = Sqrt(Power(entries_total_err[3],2)+Power(entries_total_err[4],2)+Power(entries_total_err[5],2));
    y_string_err[2] =  Sqrt(Power(entries_string_err[3],2)+Power(entries_string_err[4],2)+Power(entries_string_err[5],2));
	//Allocating graph
	TGraphErrors *gr = new TGraphErrors();
	
	//Filling graph.
	Int_t n = 0;
	for(int i = 0; i<3; i++){
		n = gr->GetN();
		gr->SetPoint(n,x[i],entries_string[i]/entries_total[i]);
		gr->SetPointError(n,0,Sqrt(Abs(Power(y_err[i]/entries_total[i],2)+Power(y_string_err[i]/entries_string[i],2))));
		file<<x[i]<<" "<<entries_string[i]/entries_total[i]<<endl;
        cout<<(entries_string[i]/entries_total[i])*Sqrt(Abs(Power(y_err[i]/entries_total[i],2)+Power(y_string_err[i]/entries_string[i],2)-2*((rho_test(3,entries_total,entries_string)*give_std(3,entries_total)*give_std(3,entries_string))/(entries_string[i]*entries_total[i]))))<<endl;
        cout<<Covariance(3,entries_total,entries_string)<<endl;
        cout<<entries_string[i]*entries_total[i]<<endl;
        cout<<entries_total[i]<<" +- "<<entries_total_err[i]<<endl;
        cout<<rho_test(3,entries_total,entries_string)<<endl;
	}
	
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(.8);
	gr->SetMarkerColor(1);
	gr->GetXaxis()->SetTitle("#it{p_{T}^{tr.}} (GeV/c)");
	gr->GetYaxis()->SetTitle("Ratio");

	file.close();
	return gr;
	
	
}

void production_ratios(){
TGraphErrors *gr = plot_ratios("/home/isidiras/university_staff/ccbar_MONASH_Soft/complete_root/LplusLminus.root","/home/isidiras/university_staff/comparisions",2);
gr->Draw("AP");



}