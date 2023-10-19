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

double_t give_mean(double_t sample[100]){
    double_t sum = 0.;
    for(int i = 0; i < 100;i++){
        sum += sample[i];
    }


    return sum/100;
}

double_t give_std(double_t sample [100]){
    double_t sum = 0.;
    double_t mean = give_mean(sample);
    for(int i = 0; i<100; i++){
        sum += pow(sample[i]-mean,2);
    }

    return sqrt(sum/99);
}

void plot_graph(){
    fstream file[100][3];

    for(int i = 0; i < 100; i++){
        for(int j = 0; j < 3; j++){
            if( j == 0 ) file[i][j].open(Form("text_output/DplusDminus_yield_low%i.txt",i+1),ios::in);
            if( j == 1 ) file[i][j].open(Form("text_output/DplusDminus_yield_intermediate%i.txt",i+1),ios::in);
            if( j == 2 ) file[i][j].open(Form("text_output/DplusDminus_yield_High%i.txt",i+1),ios::in);
            if(!file[i][j].is_open()) cout<<"Problem found in "<<i<<","<<j<<endl;
        }
    }
    cout<<"Files opened succesfully!"<<endl;

    double_t x_low[100][3];
    double_t x_int[100][2];
    double_t x_high[100];

    double_t y_low[100][3];
    double_t y_int[100][2];
    double_t y_high[100];

    for(int i = 0; i < 100; i++){
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
    double_t x1[100];
    double_t x2[100];
    double_t x3[100];
    double_t x4[100];
    double_t x5[100];
    double_t y1[100];
    double_t y2[100];
    double_t y3[100];
    double_t y4[100];
    double_t y5[100];

    for(int i = 0; i < 100; i++){
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

    x[0] = give_mean(x1);
    x[1] = give_mean(x2);
    x[2] = give_mean(x3);
    x[3] = give_mean(x4);
    x[4] = give_mean(x5);
    x[5] = give_mean(x_high);

    y[0] = give_mean(y1);
    y[1] = give_mean(y2);
    y[2] = give_mean(y3);
    y[3] = give_mean(y4);
    y[4] = give_mean(y5);
    y[5] = give_mean(y_high);

    x_err[0] = give_std(x1);
    x_err[1] = give_std(x2);
    x_err[2] = give_std(x3);
    x_err[3] = give_std(x4);
    x_err[4] = give_std(x5);
    x_err[5] = give_std(x_high);

    y_err[0] = give_std(y1);
    y_err[1] = give_std(y2);
    y_err[2] = give_std(y3);
    y_err[3] = give_std(y4);
    y_err[4] = give_std(y5);
    y_err[5] = give_std(y_high);


cout<<y[0]<<" +- "<<y_err[0]<<endl;
    TGraphErrors *gr = new TGraphErrors(6,x,y,x_err,y_err);
    gr->SetTitle("D^{+}D^{-} Away Side yield");
    gr->SetMarkerStyle(20);
    gr->Draw("AP");
    




 
}