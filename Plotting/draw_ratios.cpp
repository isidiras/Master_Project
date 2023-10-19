#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
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
TGraphErrors *plot_text(int n, const char* filename,Int_t color){
    //Opening files
    fstream file[n];
    for(int i =  0; i < n; i++){
        file[i].open(Form("%sdata%i.txt",filename,i+1),ios::in);
        if(!file[i].is_open()) cout<<"Problem found in "<<i<<endl;
    }
    cout<<"Files opened succesfully!"<<endl;
    //Array allocation
    Double_t x[n][3];
    Double_t y[n][3];
    //Filling arrays
    for(int i = 0; i<n; i++){
        int z = 0;
        Double_t xt=0;
        Double_t yt=0;
        while(1){

            file[i]>>xt>>yt;
            if(file[i].eof()) break;
            x[i][z]=xt;
            y[i][z]=yt;
            cout<<x[i][z]<<" "<<y[i][z]<<endl;
            z++;
            
        }

    }
    //1d Arrays
    Double_t x1[10];
    Double_t x2[10];
    Double_t x3[10];
    Double_t y1[10];
    Double_t y2[10];
    Double_t y3[10];
    for(int i = 0; i<n; i++){
        x1[i]=x[i][0];
        x2[i]=x[i][1];
        x3[i]=x[i][2];
        y1[i]=y[i][0];
        y2[i]=y[i][1];
        y3[i]=y[i][2];
        file[i].close();
        cout<<"________________________"<<endl;
        cout<<x1[i]<<" "<<y1[i]<<endl;
        cout<<x2[i]<<" "<<y2[i]<<endl;
        cout<<x3[i]<<" "<<y3[i]<<endl;
        cout<<"________________________"<<endl;
    }
    //Final calculsation
    Double_t x_final[3];
    Double_t y_final[3];
    Double_t x_err[3];
    Double_t y_err[3];
    //Values
    x_final[0] = give_mean(n,x1);
    x_final[1] = give_mean(n,x2);
    x_final[2] = give_mean(n,x3);
    y_final[0] = give_mean(n,y1);
    y_final[1] = give_mean(n,y2);
    y_final[2] = give_mean(n,y3);
    //Errors
    x_err[0] = give_std(n,x1);
    x_err[1] = give_std(n,x2);
    x_err[2] = give_std(n,x3);
    y_err[0] = give_std(n,y1);
    y_err[1] = give_std(n,y2);
    y_err[2] = give_std(n,y3);

    for(int i = 0; i < 3; i++){
    cout<<"y= "<<y_final[i]<<" +- "<<y_err[i]<<" x= "<<x_final[i]<<" +- "<<x_err[i]<<endl;
}
    TGraphErrors *gr = new TGraphErrors(3,x_final,y_final,x_err,y_err);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.5);
    gr->SetMarkerColor(color);
    gr->SetTitle("");
    gr->GetYaxis()->SetTitle("Ratio");
    gr->GetXaxis()->SetTitle("p_{T}^{tr.} (GeV/c)");
    gr->GetYaxis()->SetLabelSize(0.055);
    gr->GetYaxis()->SetTitleSize(0.055);
    gr->GetYaxis()->SetTitleOffset(1.6);
    gr->GetXaxis()->SetLabelSize(0.055);
    gr->GetXaxis()->SetTitleSize(0.055);
    gr->GetXaxis()->SetTitleOffset(1.0);
    gr->GetXaxis()->SetLimits(0.15,16.1);
    return gr;
}



void draw_ratios(){
    TGraphErrors *gr1 = plot_text(10,"ccbar_Monash_Hard_Low/text_ratio/DplusDminus/",1);
    TGraphErrors *gr2 = plot_text(10,"ccbar_Monash_Hard_High/text_ratio/DplusDminus/",2);
    TGraphErrors *gr3 = plot_text(10,"ccbar_Monash_Soft/text_ratio/DplusDminus/",3);
    TGraphErrors *gr4 = plot_text(10,"ccbar_Junctions_Hard_Low/text_ratio/DplusDminus/",4);
    TGraphErrors *gr5 = plot_text(10,"ccbar_Junctions_Hard_High/text_ratio/DplusDminus/",5);
    TGraphErrors *gr6 = plot_text(10,"ccbar_Junctions_Soft/text_ratio/DplusDminus/",6);
    Double_t ymin = 0.4;
    Double_t ymax = 1.15;
    gr1->GetYaxis()->SetRangeUser(ymin,ymax);
    TLine *line_one_w = new TLine(3.,ymin,3.,ymax);
  	line_one_w->SetLineStyle(2);
  	line_one_w->SetLineWidth(1);
  	
  	TLine *line_two_w = new TLine(8.,ymin,8.,ymax);
  	line_two_w->SetLineStyle(2);
  	line_two_w->SetLineWidth(1);

    TLegend *legM = new TLegend(0.25,0.7,0.45,0.88);
     legM->SetBorderSize(0);
     legM->SetFillColor(0);
     legM->SetTextFont(42);
     legM->SetTextSize(0.035);
     legM->SetHeader("MONASH","C");
     legM->AddEntry(gr1,"pTHatMin = 1 GeV/c","P");
     legM->AddEntry(gr2,"pTHatMin = 10 GeV/c","P");
     legM->AddEntry(gr3,"SoftQCD","P");

     TLegend *legJ = new TLegend(0.6,0.7,0.8,0.88);
     legJ->SetBorderSize(0);
     legJ->SetFillColor(0);
     legJ->SetTextFont(42);
     legJ->SetTextSize(0.035);
     legJ->SetHeader("Junctions","C");
     legJ->AddEntry(gr4,"pTHatMin = 1 GeV/c","P");
     legJ->AddEntry(gr5,"pTHatMin = 10 GeV/c","P");
     legJ->AddEntry(gr6,"SoftQCD","P");

    TCanvas *c1 = new TCanvas("c1","Production Ratio",1000,700);
    c1->SetLeftMargin(0.18);
    c1->SetRightMargin(0.05);
    c1->SetBottomMargin(0.14);
    gr1->Draw("AP");
    gr2->Draw("SAME P");
    gr3->Draw("SAME P");
    gr4->Draw("SAME P");
    gr5->Draw("SAME P");
    gr6->Draw("SAME P");
    line_one_w->Draw();
    line_two_w->Draw();
    legM->Draw();
    legJ->Draw();

}