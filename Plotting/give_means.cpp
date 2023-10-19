//This scripts calculates the pT average and STD of the trigger
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
    input->Close();
    
}

void give_means(){
    const char* filepath[6] = {"/home/isidiras/university_staff/ccbar_MONASH_Hard_low/complete_root/","/home/isidiras/university_staff/ccbar_MONASH_Hard_High/complete_root/","/home/isidiras/university_staff/ccbar_MONASH_Soft/complete_root/","/home/isidiras/university_staff/ccbar_Junctions_Hard_Low/complete_root/","/home/isidiras/university_staff/ccbar_Junctions_Hard_High/complete_root/","/home/isidiras/university_staff/ccbar_Junctions_Soft/complete_root/"};
    TH1D *LSpectra[6];
    TH1D *DSpectra[6];
    //Read_Histograms
    for(int i = 0; i<6; i++){
        LSpectra[i] = Read_Hist(Form("%sLplusLminus.root",filepath[i]),"hTrPt");
	    DSpectra[i] = Read_Hist(Form("%sDplusDminus.root",filepath[i]),"hTrPt");
        cout<<"_________________________"<<endl;
        cout<<"Tune: "<<i+1<<endl;
        cout<<"For Lambda:"<<endl;
        LSpectra[i]->GetXaxis()->SetRangeUser(1,3);
        cout<<"For the Low: "<<LSpectra[i]->GetMean()<<" +- "<<LSpectra[i]->GetStdDev()<<endl;
        LSpectra[i]->GetXaxis()->SetRangeUser(3,8);
        cout<<"For the Intermediate: "<<LSpectra[i]->GetMean()<<" +- "<<LSpectra[i]->GetStdDev()<<endl;
        LSpectra[i]->GetXaxis()->SetRangeUser(8,50);
        cout<<"For the High: "<<LSpectra[i]->GetMean()<<" +- "<<LSpectra[i]->GetStdDev()<<endl;
        cout<<"For D:"<<endl;
        DSpectra[i]->GetXaxis()->SetRangeUser(1,3);
        cout<<"For the Low: "<<DSpectra[i]->GetMean()<<" +- "<<DSpectra[i]->GetStdDev()<<endl;
        DSpectra[i]->GetXaxis()->SetRangeUser(3,8);
        cout<<"For the Intermediate: "<<DSpectra[i]->GetMean()<<" +- "<<DSpectra[i]->GetStdDev()<<endl;
        DSpectra[i]->GetXaxis()->SetRangeUser(8,50);
        cout<<"For the High: "<<DSpectra[i]->GetMean()<<" +- "<<DSpectra[i]->GetStdDev()<<endl;
    }


}