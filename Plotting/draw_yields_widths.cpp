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

TGraphErrors *plot_graph(int n,const char* filename,Int_t color,const char* title){
    fstream file[n][3];

    for(int i = 0; i < n; i++){
        for(int j = 0; j < 3; j++){
            if( j == 0 ) file[i][j].open(Form("%s/low%i.txt",filename,i+1),ios::in);
            if( j == 1 ) file[i][j].open(Form("%s/intermediate%i.txt",filename,i+1),ios::in);
            if( j == 2 ) file[i][j].open(Form("%s/high%i.txt",filename,i+1),ios::in);
            if(!file[i][j].is_open()) cout<<"Problem found in "<<i<<","<<j<<endl;
        }
    }
    cout<<"Files opened succesfully!"<<endl;

    double_t x_low[n][3];
    double_t x_int[n][2];
    double_t x_high[n];

    double_t y_low[n][3];
    double_t y_int[n][2];
    double_t y_high[n];

    for(int i = 0; i < n; i++){
        for(int j = 0; j < 3; j++){
            if( j == 0 ){
                int z = 0;
                while(z < 3){
                    double_t x, y;
                    file[i][j]>>x>>y;
                   // if(i%10 == 0) cout<<x<<" "<<y<<endl;
                    x_low[i][z] = x;
                    y_low[i][z] = y;
                    cout<<x_low[i][z]<<" "<<y_low[i][z]<<endl;
                    if(file[i][j].eof()) break;
                    z++;
                }
            }
            if( j == 1 ){
                int z = 0;
                while(z < 2){
                    double_t x, y;
                    file[i][j]>>x>>y;
                   // if(i%10 == 0) cout<<x<<" "<<y<<endl;
                    x_int[i][z] = x;
                    y_int[i][z] = y;
                  //  if(i%10 == 0) cout<<x_int[i][z]<<" "<<y_int[i][z]<<endl;
                    z++;
                    
                }
            }

            if( j == 2 ){
                while(!file[i][j].eof()){
                    double_t x, y;
                    file[i][j]>>x>>y;
                  //  if(i%10 == 0) cout<<x<<" "<<y<<endl;
                    x_high[i] = x;
                    y_high[i] = y;
                   // if(i%10 == 0) cout<<x_high[i]<<" "<<y_high[i]<<endl;
                    
                }
            }


        }
    }
    double_t x1[n];
    double_t x2[n];
    double_t x3[n];
    double_t x4[n];
    double_t x5[n];
    double_t y1[n];
    double_t y2[n];
    double_t y3[n];
    double_t y4[n];
    double_t y5[n];

    for(int i = 0; i < n; i++){
        x1[i] = x_low[i][0];
        x2[i] = x_low[i][1];
        x3[i] = x_low[i][2];
        x4[i] = x_int[i][0];
        x5[i] = x_int[i][1];
        y1[i] = y_low[i][0];
        y2[i] = y_low[i][1];
        y3[i] = y_low[i][2];
        y4[i] = y_int[i][0];
        y5[i] = y_int[i][1];

        cout<<y1[i]<<endl;

        for(int j = 0; j < 3; j++){
            file[i][j].close();
        }
    }
    cout<<"Files closed succesfully!"<<endl;

    double_t x[6];
    double_t y[6];

    double_t x_err[6];
    double_t y_err[6];

    x[0] = give_mean(n,x1);
    x[1] = give_mean(n,x2);
    x[2] = give_mean(n,x3);
    x[3] = give_mean(n,x4);
    x[4] = give_mean(n,x5);
    x[5] = give_mean(n,x_high);

    y[0] = give_mean(n,y1);
    y[1] = give_mean(n,y2);
    y[2] = give_mean(n,y3);
    y[3] = give_mean(n,y4);
    y[4] = give_mean(n,y5);
    y[5] = give_mean(n,y_high);

    x_err[0] = give_std(n,x1);
    x_err[1] = give_std(n,x2);
    x_err[2] = give_std(n,x3);
    x_err[3] = give_std(n,x4);
    x_err[4] = give_std(n,x5);
    x_err[5] = give_std(n,x_high);

    y_err[0] = give_std(n,y1);
    y_err[1] = give_std(n,y2);
    y_err[2] = give_std(n,y3);
    y_err[3] = give_std(n,y4);
    y_err[4] = give_std(n,y5);
    y_err[5] = give_std(n,y_high);


for(int i = 0; i<n;i++){
    cout<<i<<" "<<x1[i]<<", "<<y1[i]<<endl;
    cout<<i<<" "<<x2[i]<<", "<<y2[i]<<endl;
    cout<<i<<" "<<x3[i]<<", "<<y3[i]<<endl;
    cout<<i<<" "<<x4[i]<<", "<<y4[i]<<endl;
    cout<<i<<" "<<x5[i]<<", "<<y5[i]<<endl;
    cout<<i<<" "<<x_high[i]<<", "<<y_high[i]<<endl;
    cout<<"__________________"<<endl;

}

for(int i = 0; i < 6; i++){
    cout<<"y= "<<y[i]<<" +- "<<y_err[i]<<" x= "<<x[i]<<" +- "<<x_err[i]<<endl;
}
    TGraphErrors *gr = new TGraphErrors(6,x,y,x_err,y_err);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.5);
    gr->SetMarkerColor(color);
    gr->SetTitle("");
    gr->GetYaxis()->SetTitle(title);
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

void draw_yields_widths(){
    //Graph allocation
    //SameSideYield
    TGraphErrors *gr1 = plot_graph(10,"ccbar_Monash_Hard_Low/text_outputM/DplusDminus/SameSideYield/",1,"Yield");
    TGraphErrors *gr2 = plot_graph(10,"ccbar_Monash_Hard_High/text_outputM/DplusDminus/SameSideYield/",2,"Yield");
    TGraphErrors *gr3 = plot_graph(10,"ccbar_Monash_Soft/text_outputM/DplusDminus/SameSideYield/",1,"Yield");
    TGraphErrors *gr4 = plot_graph(10,"ccbar_Junctions_Hard_Low/text_outputM/DplusDminus/SameSideYield/",4,"Yield");
    TGraphErrors *gr5 = plot_graph(10,"ccbar_Junctions_Hard_High/text_outputM/DplusDminus/SameSideYield/",5,"Yield");
    TGraphErrors *gr6 = plot_graph(10,"ccbar_Junctions_Soft/text_outputM/DplusDminus/SameSideYield/",6,"Yield");
    gr1->GetYaxis()->SetRangeUser(0,0.003);

    //AwaySideYield
    TGraphErrors *grA1 = plot_graph(10,"ccbar_Monash_Hard_Low/text_outputM/DplusDplus/AwaySideYield/",1,"Yield");
    TGraphErrors *grA2 = plot_graph(10,"ccbar_Monash_Hard_High/text_outputM/DplusDplus/AwaySideYield/",2,"Yield");
    TGraphErrors *grA3 = plot_graph(10,"ccbar_Monash_Soft/text_outputM/DplusDminus/AwaySideYield/",1,"Yield");
    TGraphErrors *grA4 = plot_graph(10,"ccbar_Junctions_Hard_Low/text_outputM/DplusDplus/AwaySideYield/",4,"Yield");
    TGraphErrors *grA5 = plot_graph(10,"ccbar_Junctions_Hard_High/text_outputM/DplusDplus/AwaySideYield/",5,"Yield");
    TGraphErrors *grA6 = plot_graph(10,"ccbar_Junctions_Soft/text_outputM/DplusDplus/AwaySideYield/",6,"Yield");
    grA3->GetYaxis()->SetRangeUser(0,0.05);

    //SameSideWidth
    TGraphErrors *grW1 = plot_graph(5,"ccbar_Monash_Hard_Low/text_outputL/DplusDminus/SameSideWidth/",1,"Width (rad)");
    TGraphErrors *grW2 = plot_graph(10,"ccbar_Monash_Hard_High/text_outputM/DplusDminus/SameSideWidth/",2,"Width (rad)");
    TGraphErrors *grW3 = plot_graph(10,"ccbar_Monash_Soft/text_outputM/DplusDminus/SameSideWidth/",3,"Width (rad)");
    TGraphErrors *grW4 = plot_graph(10,"ccbar_Junctions_Hard_Low/text_outputM/DplusDminus/SameSideWidth/",4,"Width (rad)");
    TGraphErrors *grW5 = plot_graph(10,"ccbar_Junctions_Hard_High/text_outputM/DplusDminus/SameSideWidth/",5,"Width (rad)");
    TGraphErrors *grW6 = plot_graph(10,"ccbar_Junctions_Soft/text_outputM/DplusDminus/SameSideWidth/",6,"Width (rad)");
    grW1->GetYaxis()->SetRangeUser(0.0001,2);

    //AwaySideWidth
    TGraphErrors *grWA1 = plot_graph(5,"ccbar_Monash_Hard_Low/text_outputL/DplusDminus/AwaySideWidth/",1,"Width (rad)");
    TGraphErrors *grWA2 = plot_graph(5,"ccbar_Monash_Hard_High/text_outputL/DplusDminus/AwaySideWidth/",2,"Width (rad)");
    TGraphErrors *grWA3 = plot_graph(5,"ccbar_Monash_Soft/text_outputL/DplusDminus/AwaySideWidth/",3,"Width (rad)");
    TGraphErrors *grWA4 = plot_graph(5,"ccbar_Junctions_Hard_Low/text_outputL/DplusDminus/AwaySideWidth/",4,"Width (rad)");
    TGraphErrors *grWA5 = plot_graph(5,"ccbar_Junctions_Hard_High/text_outputL/DplusDminus/AwaySideWidth/",5,"Width (rad)");
    TGraphErrors *grWA6 = plot_graph(5,"ccbar_Junctions_Soft/text_outputL/DplusDminus/AwaySideWidth/",6,"Width (rad)");
    grWA1->GetYaxis()->SetRangeUser(0.0,2.5);


    //Plotting
	Double_t ymax = 0.003;//SameSideYield
	Double_t ymax_w = 2;//SameSideWidth
    Double_t ymaxA = 0.05;//AwaySideYield
    Double_t ymax_wA = 2.5; //AwaySideWidth
	//SameSide axis
	TGaxis *low = new TGaxis(8.,ymax,16.1,ymax,0.15,16.1,10,"-L");
	low->SetLabelColor(kRed+2);
 	low->SetLineColor(kRed+2);
  	low->SetTitleColor(kRed+2);
  	low->CenterTitle();
  	low->SetLabelSize(0.035);
  	low->SetTitleSize(0.03);
  	//low->SetTitle("#it{p_{T}^{tr}} (GeV/c)");
  	
  	TGaxis *inter = new TGaxis(3.,ymax,8.,ymax,0.15,8.,10,"-L");
	inter->SetLabelColor(kRed+2);
 	inter->SetLineColor(kRed+2);
  	inter->SetTitleColor(kRed+2);
  	inter->CenterTitle();
  	inter->SetLabelSize(0.035);
  	inter->SetTitleSize(0.035);
  	inter->SetTitle("#it{p_{T}^{as.}} (GeV/c)");
  	
  	TGaxis *high = new TGaxis(0.15,ymax,3.,ymax,0.15,3.1,5,"-L");
	high->SetLabelColor(kRed+2);
 	high->SetLineColor(kRed+2);
  	high->SetTitleColor(kRed+2);
  	high->CenterTitle();
  	high->SetLabelSize(0.035);
  	high->SetTitleSize(0.035);
  	//high->SetTitle("#it{p_{T}^{tr}} (GeV/c)");
  	
  	TLine *line_one = new TLine(3.,0.,3.,ymax);
  	line_one->SetLineStyle(2);
  	line_one->SetLineWidth(1);
  	
  	TLine *line_two = new TLine(8.,0.,8.,ymax);
  	line_two->SetLineStyle(2);
  	line_two->SetLineWidth(1);
    //Away Side Axis
    TGaxis *lowA = new TGaxis(8.,ymaxA,16.1,ymaxA,0.15,16.1,10,"-L");
	lowA->SetLabelColor(kRed+2);
 	lowA->SetLineColor(kRed+2);
  	lowA->SetTitleColor(kRed+2);
  	lowA->CenterTitle();
  	lowA->SetLabelSize(0.035);
  	lowA->SetTitleSize(0.035);
  	//low->SetTitle("#it{p_{T}^{tr}} (GeV/c)");
  	
  	TGaxis *interA = new TGaxis(3.,ymaxA,8.,ymaxA,0.15,8.,10,"-L");
	interA->SetLabelColor(kRed+2);
 	interA->SetLineColor(kRed+2);
  	interA->SetTitleColor(kRed+2);
  	interA->CenterTitle();
  	interA->SetLabelSize(0.035);
  	interA->SetTitleSize(0.035);
  	interA->SetTitle("#it{p_{T}^{as.}} (GeV/c)");
  	
  	TGaxis *highA = new TGaxis(0.15,ymaxA,3.,ymaxA,0.15,3.1,5,"-L");
	highA->SetLabelColor(kRed+2);
 	highA->SetLineColor(kRed+2);
  	highA->SetTitleColor(kRed+2);
  	highA->CenterTitle();
  	highA->SetLabelSize(0.035);
  	highA->SetTitleSize(0.035);
  	//high->SetTitle("#it{p_{T}^{tr}} (GeV/c)");
  	
  	TLine *line_oneA = new TLine(3.,0.,3.,ymaxA);
  	line_oneA->SetLineStyle(2);
  	line_oneA->SetLineWidth(1);
  	
  	TLine *line_twoA = new TLine(8.,0.,8.,ymaxA);
  	line_twoA->SetLineStyle(2);
  	line_twoA->SetLineWidth(1);
//Lines/axis for width
TGaxis *low_w = new TGaxis(0.15,ymax_w,3.,ymax_w,0.15,3.1,3,"-L");
	low_w->SetLabelColor(kRed+2);
 	low_w->SetLineColor(kRed+2);
  	low_w->SetTitleColor(kRed+2);
  	low_w->CenterTitle();
  	low_w->SetLabelSize(0.035);
  	low_w->SetTitleSize(0.035);
  	//low->SetTitle("#it{p_{T}^{tr}} (GeV/c)");
  	
	TGaxis *inter_w = new TGaxis(3.,ymax_w,8.,ymax_w,0.15,8.,10,"-L");
	inter_w->SetLabelColor(kRed+2);
 	inter_w->SetLineColor(kRed+2);
  	inter_w->SetTitleColor(kRed+2);
  	inter_w->CenterTitle();
  	inter_w->SetLabelSize(0.035);
  	inter_w->SetTitleSize(0.035);
  	inter_w->SetTitle("#it{p_{T}^{as.}} (GeV/c)");
  	
  	TGaxis *high_w = new TGaxis(8.,ymax_w,16.1,ymax_w,0.15,16.1,10,"-L");
	high_w->SetLabelColor(kRed+2);
 	high_w->SetLineColor(kRed+2);
  	high_w->SetTitleColor(kRed+2);
  	high_w->CenterTitle();
  	high_w->SetLabelSize(0.035);
  	high_w->SetTitleSize(0.035);
  	//high->SetTitle("#it{p_{T}^{tr}} (GeV/c)");
  	
  	TLine *line_one_w = new TLine(3.,0.,3.,ymax_w);
  	line_one_w->SetLineStyle(2);
  	line_one_w->SetLineWidth(1);
  	
  	TLine *line_two_w = new TLine(8.,0.,8.,ymax_w);
  	line_two_w->SetLineStyle(2);
  	line_two_w->SetLineWidth(1);
    //Axis Away Side Width
    TGaxis *low_wA = new TGaxis(0.15,ymax_wA,3.,ymax_wA,0.15,3.1,3,"-L");
	low_wA->SetLabelColor(kRed+2);
 	low_wA->SetLineColor(kRed+2);
  	low_wA->SetTitleColor(kRed+2);
  	low_wA->CenterTitle();
  	low_wA->SetLabelSize(0.035);
  	low_wA->SetTitleSize(0.035);
  	//low->SetTitle("#it{p_{T}^{tr}} (GeV/c)");
  	
	TGaxis *inter_wA = new TGaxis(3.,ymax_wA,8.,ymax_wA,0.15,8.,10,"-L");
	inter_wA->SetLabelColor(kRed+2);
 	inter_wA->SetLineColor(kRed+2);
  	inter_wA->SetTitleColor(kRed+2);
  	inter_wA->CenterTitle();
  	inter_wA->SetLabelSize(0.035);
  	inter_wA->SetTitleSize(0.035);
  	inter_wA->SetTitle("#it{p_{T}^{as.}} (GeV/c)");
  	
  	TGaxis *high_wA = new TGaxis(8.,ymax_wA,16.1,ymax_wA,0.15,16.1,10,"-L");
	high_wA->SetLabelColor(kRed+2);
 	high_wA->SetLineColor(kRed+2);
  	high_wA->SetTitleColor(kRed+2);
  	high_wA->CenterTitle();
  	high_wA->SetLabelSize(0.035);
  	high_wA->SetTitleSize(0.035);
  	//high->SetTitle("#it{p_{T}^{tr}} (GeV/c)");
  	
  	TLine *line_one_wA = new TLine(3.,0.,3.,ymax_wA);
  	line_one_wA->SetLineStyle(2);
  	line_one_wA->SetLineWidth(1);
  	
  	TLine *line_two_wA = new TLine(8.,0.,8.,ymax_wA);
  	line_two_wA->SetLineStyle(2);
  	line_two_wA->SetLineWidth(1);

    //Legends
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


    //Drawing
  	TCanvas *c1 = new TCanvas("c1","Same_Side_Yield",1200,900);
    c1->SetLeftMargin(0.18);
    c1->SetRightMargin(0.05);
    c1->SetBottomMargin(0.14);
    //gr1->Draw("AP");
    //gr2->Draw("SAME P");
    gr3->Draw("AP");
   // gr4->Draw("SAME P");
   // gr5->Draw("SAME P");
   // gr6->Draw("SAME P");
    low->Draw();
    inter->Draw();
    high->Draw();
    line_one->Draw();
    line_two->Draw();
    //legM->Draw();
    //legJ->Draw();
    
    TCanvas *c2 = new TCanvas("c2","Away_Side_Yield",1200,900);
    c2->SetLeftMargin(0.18);
    c2->SetRightMargin(0.05);
    c2->SetBottomMargin(0.14);
    //grA1->Draw("AP");
    //grA2->Draw("SAME P");
    grA3->Draw("AP");
    //grA4->Draw("SAME P");
    //grA5->Draw("SAME P");
   // grA6->Draw("SAME P");
    lowA->Draw();
    interA->Draw();
    highA->Draw();
    line_oneA->Draw();
    line_twoA->Draw();
    //legM->Draw();
    //legJ->Draw();

    TCanvas *c3 = new TCanvas("c3","Same_Side_Width",1200,900);
    c3->SetLeftMargin(0.18);
    c3->SetRightMargin(0.05);
    c3->SetBottomMargin(0.14);
    grW1->Draw("AP");
    grW2->Draw("SAME P");
    grW3->Draw("SAME P");
    grW4->Draw("SAME P");
    grW5->Draw("SAME P");
    grW6->Draw("SAME P");
    low_w->Draw();
    inter_w->Draw();
    high_w->Draw();
    line_one_w->Draw();
    line_two_w->Draw();
    legM->Draw();
    legJ->Draw();

    TCanvas *c4 = new TCanvas("c4","Away_Side_Width",1200,900);
    c4->SetLeftMargin(0.18);
    c4->SetRightMargin(0.05);
    c4->SetBottomMargin(0.14);
    grWA1->Draw("AP");
    grWA2->Draw("SAME P");
    grWA3->Draw("SAME P");
    grWA4->Draw("SAME P");
    grWA5->Draw("SAME P");
    grWA6->Draw("SAME P");
    low_wA->Draw();
    inter_wA->Draw();
    high_wA->Draw();
    line_one_wA->Draw();
    line_two_wA->Draw();
    legM->Draw();
    legJ->Draw();




}