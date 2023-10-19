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
	hist->GetYaxis()->SetTitle("#frac{1}{N_{tr}} #frac{dN_{as.}}{d#Delta#phi}");
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

TH1D *reduction(TH1D *hist, Double_t start, Double_t finish){
//This function does the whole procces
	double *p = Zyam_and_Error(start,finish,hist);
	TH1D *Signal = signal(hist,*p,*(p+1));
	cout<<"The ZYAM background is: "<<*p<<endl;

	delete p;
	return Signal; 
}

void plot_all(const char* pair){
	//Dictionaries used
	const char* DPhi_Names[6] = {"hDPhiLL","hDPhiIL","hDPhiII","hDPhiHL","hDPhiHI","hDPhiHH"};
    const char* filepath[6] = {"/home/isidiras/university_staff/ccbar_MONASH_Hard_low/complete_root/","/home/isidiras/university_staff/ccbar_MONASH_Hard_High/complete_root/","/home/isidiras/university_staff/ccbar_MONASH_Soft/complete_root/","/home/isidiras/university_staff/ccbar_Junctions_Hard_Low/complete_root/","/home/isidiras/university_staff/ccbar_Junctions_Hard_High/complete_root/","/home/isidiras/university_staff/ccbar_Junctions_Soft/complete_root/"};
    //Histogram arrays
    TH1D* hDPhi[6][6];
    TH1D* reduced[6][6];
    TH1D* hTriggers[6];
    for(int i = 0; i<6;i++){
        hTriggers[i]=Read_Hist(Form("%s%s",filepath[i],pair),"hTrPt");
        for(int j = 0; j<6; j++){
            hDPhi[i][j] = Read_Hist(Form("%s%s",filepath[i],pair),DPhi_Names[j]);
        }
    }
    //Normalizing
    for(int j = 0; j<6; j++){
        for (int i = 0; i < 6; i++){
		    if(i == 0){
			    normalize(hDPhi[j][i],hTriggers[j],1.,3.);
            }
            if(i == 1 || i == 2){
			    normalize(hDPhi[j][i],hTriggers[j],3.,8.);
            }
            if(i == 3 || i == 4 || i == 5){
			    normalize(hDPhi[j][i],hTriggers[j],8.,50.);
            }

    
        }
    }

    //Plot cosmetics
    for(int j = 0; j<6; j++){
        for(int i=0;i<6;i++){
            reduced[j][i] = reduction(hDPhi[j][i],hDPhi[j][i]->GetBinCenter(hDPhi[j][i]->GetMinimumBin())-0.15,hDPhi[j][i]->GetBinCenter(hDPhi[j][i]->GetMinimumBin())+0.15);
            
            reduced[j][i]->SetMarkerStyle(20);
            reduced[j][i]->SetMarkerSize(1);
           // reduced[j][i]->SetMarkerColor(j+1);
            reduced[j][i]->GetYaxis()->SetLabelSize(0.055);
            reduced[j][i]->GetYaxis()->SetTitleSize(0.055);
            reduced[j][i]->GetYaxis()->SetTitleOffset(2.2);
            reduced[j][i]->GetXaxis()->SetLabelSize(0.055);
            reduced[j][i]->GetXaxis()->SetTitleSize(0.055);
            reduced[j][i]->GetXaxis()->SetTitleOffset(1.0);
            reduced[j][i]->SetStats(0);
            reduced[j][i]->SetTitle(" ");

        }
    }

        //Legends
    TLegend *legM = new TLegend(0.25,0.7,0.45,0.88);
     legM->SetBorderSize(0);
     legM->SetFillColor(0);
     legM->SetTextFont(42);
     legM->SetTextSize(0.055);
     legM->SetHeader("MONASH","C");
     legM->AddEntry(reduced[0][0],"pTHatMin = 1 GeV/c","P");
     legM->AddEntry(reduced[1][0],"pTHatMin = 10 GeV/c","P");
     legM->AddEntry(reduced[2][0],"SoftQCD","P");

     TLegend *legJ = new TLegend(0.6,0.7,0.8,0.88);
     legJ->SetBorderSize(0);
     legJ->SetFillColor(0);
     legJ->SetTextFont(42);
     legJ->SetTextSize(0.055);
     legJ->SetHeader("Junctions","C");
     legJ->AddEntry(reduced[3][0],"pTHatMin = 1 GeV/c","P");
     legJ->AddEntry(reduced[4][0],"pTHatMin = 10 GeV/c","P");
     legJ->AddEntry(reduced[5][0],"SoftQCD","P");

     //Drawing
    TCanvas *c1 = new TCanvas();
    c1->Divide(3,2);
    c1->cd(1);
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(0.28);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.13);
    //reduced[0][0]->GetYaxis()->SetRangeUser(0.0008,0.003);
   // reduced[0][0]->Draw("PE1");
    //reduced[1][0]->Draw("PE1 SAME");
    reduced[2][0]->Draw("PE1");
    //reduced[3][0]->Draw("PE1 SAME");
    //reduced[4][0]->Draw("PE1 SAME");
    //reduced[5][0]->Draw("PE1 SAME");
    TLatex *text = new TLatex();
    text->SetTextFont(42);
    text->SetTextSize(0.07);
    text->DrawLatexNDC(0.5,0.77,"#splitline{1<p^{tr.}_{T}<3 GeV/c}{1<p^{as.}_{T}<3 GeV/c}");
  
   

    c1->cd(2);
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(0.28);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.13);
    //reduced[0][1]->GetYaxis()->SetRangeUser(0.0007,0.003);
    //reduced[0][1]->Draw("PE1");
    //reduced[1][1]->Draw("PE1 SAME");
    reduced[2][1]->Draw("PE1");
    //reduced[3][1]->Draw("PE1 SAME");
    //reduced[4][1]->Draw("PE1 SAME");
   // reduced[5][1]->Draw("PE1 SAME");
    TLatex *text1 = new TLatex();
    text1->SetTextFont(42);
    text1->SetTextSize(0.07);
    text->DrawLatexNDC(0.5,0.77,"#splitline{3<p^{tr.}_{T}<8 GeV/c}{1<p^{as.}_{T}<3 GeV/c}");
    
  


    c1->cd(3);
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(0.28);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.13);
   // reduced[0][2]->GetYaxis()->SetRangeUser(0.0003,0.0038);
    //reduced[0][2]->Draw("PE1");
   // reduced[1][2]->Draw("PE1 SAME");
    reduced[2][2]->Draw("PE1");
    //reduced[3][2]->Draw("PE1 SAME");
    //reduced[4][2]->Draw("PE1 SAME");
    //reduced[5][2]->Draw("PE1 SAME");
    TLatex *text2 = new TLatex();
    text2->SetTextFont(42);
    text2->SetTextSize(0.07);
    text2->DrawLatexNDC(0.5,0.77,"#splitline{3<p^{tr.}_{T}<8 GeV/c}{3<p^{as.}_{T}<8 GeV/c}");
 

    c1->cd(4);
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(0.28);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.13);
    //reduced[0][3]->GetYaxis()->SetRangeUser(0.0006,0.0028);
    //reduced[0][3]->Draw("PE1");
    //reduced[1][3]->Draw("PE1 SAME");
    reduced[2][3]->Draw("PE1");
    //reduced[3][3]->Draw("PE1 SAME");
    //reduced[4][3]->Draw("PE1 SAME");
    //reduced[5][3]->Draw("PE1 SAME");
    TLatex *text3 = new TLatex();
    text3->SetTextFont(42);
    text3->SetTextSize(0.07);
    text3->DrawLatexNDC(0.32,0.8,"#splitline{p^{tr.}_{T}>8 GeV/c}{1<p^{as.}_{T}<3 GeV/c}");

    c1->cd(5);
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(0.28);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.13);
    //reduced[0][4]->GetYaxis()->SetRangeUser(0.0005,0.0035);
    //reduced[0][4]->Draw("PE1");
    //reduced[1][4]->Draw("PE1 SAME");
    reduced[2][4]->Draw("PE1 ");
   // reduced[3][4]->Draw("PE1 SAME");
    //reduced[4][4]->Draw("PE1 SAME");
    //reduced[5][4]->Draw("PE1 SAME");
    TLatex *text4 = new TLatex();
    text4->SetTextFont(42);
    text4->SetTextSize(0.07);
    text4->DrawLatexNDC(0.42,0.77,"#splitline{p^{tr.}_{T}>8 GeV/c}{3<p^{as.}_{T}<8 GeV/c}");


    c1->cd(6);
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(0.28);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.13);
    //reduced[0][5]->GetYaxis()->SetRangeUser(0.00008,0.0047);
    //reduced[0][5]->Draw("PE1");
   // reduced[1][5]->Draw("PE1 SAME");
    reduced[2][5]->Draw("PE1");
    //reduced[3][5]->Draw("PE1 SAME");
    //reduced[4][5]->Draw("PE1 SAME");
    //reduced[5][5]->Draw("PE1 SAME");
    TLatex *text5 = new TLatex();
    text5->SetTextFont(42);
    text5->SetTextSize(0.07);
    text5->DrawLatexNDC(0.55,0.77,"#splitline{p^{tr.}_{T}>8 GeV/c}{p^{as.}_{T}>8 GeV/c}");
   // legM->Draw();
   // legJ->Draw();

    TCanvas *c2 = new TCanvas();
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(0.28);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.13);
    reduced[0][0]->GetYaxis()->SetRangeUser(0.0008,0.003);
    reduced[0][0]->Draw("PE1");
    reduced[1][0]->Draw("PE1 SAME");
    reduced[2][0]->Draw("PE1 SAME");
    reduced[3][0]->Draw("PE1 SAME");
    reduced[4][0]->Draw("PE1 SAME");
    reduced[5][0]->Draw("PE1 SAME");
    text->DrawLatexNDC(0.5,0.77,"#splitline{1<p^{tr.}_{T}<3 GeV/c}{1<p^{as.}_{T}<3 GeV/c}");

    TCanvas *c3 = new TCanvas();
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(0.28);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.13);
    reduced[0][5]->GetYaxis()->SetRangeUser(0.00008,0.0047);
    reduced[0][5]->Draw("PE1");
    reduced[1][5]->Draw("PE1 SAME");
    reduced[2][5]->Draw("PE1 SAME");
    reduced[3][5]->Draw("PE1 SAME");
    reduced[4][5]->Draw("PE1 SAME");
    reduced[5][5]->Draw("PE1 SAME");
    text5->DrawLatexNDC(0.55,0.77,"#splitline{p^{tr.}_{T}>8 GeV/c}{p^{as.}_{T}>8 GeV/c}");
    legM->Draw();
    legJ->Draw();

    TLegend *leg3 = new TLegend(0.25,0.7,0.45,0.88);
     leg3->SetBorderSize(0);
     leg3->SetFillColor(0);
     leg3->SetTextFont(42);
     leg3->SetTextSize(0.055);
     leg3->SetHeader("pTHatmin=10 GeV/c","C");
     leg3->AddEntry(reduced[1][2],"MONASH","P");
     leg3->AddEntry(reduced[4][2],"Junctions","P");

    TCanvas *c4 = new TCanvas();
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(0.28);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.13);
    reduced[1][2]->SetMarkerColor(1);
    reduced[4][2]->SetMarkerColor(4);
    reduced[1][2]->Draw("PE1");
    reduced[4][2]->Draw("PE1 SAME");
    text2->DrawLatexNDC(0.5,0.77,"#splitline{3<p^{tr.}_{T}<8 GeV/c}{3<p^{as.}_{T}<8 GeV/c}");
    leg3->Draw();





}

void draw_corr_back(){
    plot_all("DplusDminus.root");
}