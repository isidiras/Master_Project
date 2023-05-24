//With this script I am doing a specific analysis concerning the azimuthial correlation plots of D^{+-} mesons from 
//production mechanism. This is layed out so it can be used for other particles as well by changing the IDs. The analysis
//is done across the hole momentum ranges and it also produces 3D histograms if one wants to analyze pseudorapidity correlations
//a s a function of the trigger p_T. 
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

Double_t DeltaPhi(Double_t phi1, Double_t phi2){
	//returns delta phi in range -pi/2 3pi/2
	return fmod(phi1-phi2+2.5*PI,2*PI)-0.5*PI;
	}
	
void status_file(Int_t id_trigger,Int_t id_associate, TString filename, const char* title){
	//This functions takes the trigger and associate id and creates a .root file
	// with the name filename. The title is a string so the particles will be seen
	//in histogram titles.
	
	//Define the TChain
	TChain *ch1 = new TChain("tree");
	TFile *output = new TFile(filename,"RECREATE");
	
	int ntrees = 4; //Number of trees we want to add to the TChain This can be changed by the user.
	
	//We put the trees to the chain
	for( int i = 1; i < ntrees+1;  i++){
	//One can change the file path accordingly to his/hers machine set up
	ch1->Add(Form("/home/isidiras/university_staff/ccbar_MONASH/output_MONASH_STATUS/Group%i/output.root",i));//File is for my local set up!
	}
	
	//Now we define vectors that carry the information at event level.
	vector<Int_t>* vID = 0;
	vector<Double_t>* vPt = 0;
	vector<Double_t>* vPhi = 0;
	vector<Double_t>* vStatus = 0;
	vector<Double_t>* vEta = 0;
	//Setting up chain branch addresses to the vectors defined above
	ch1->SetBranchAddress("ID",&vID);
	ch1->SetBranchAddress("PT",&vPt);
	ch1->SetBranchAddress("PHI",&vPhi);
	ch1->SetBranchAddress("ETA",&vEta);
	ch1->SetBranchAddress("STATUS",&vStatus);
	
	//Definition of variables I am going to use
	Int_t aID,pID;
	Double_t pPt,pPhi,pStatus,pEta;//For triger
	Double_t aPt,aPhi,aStatus,aEta;//For associate
	int nTrigger = 0;
	
	//These two can be changed from the user according to the status histogram
	Double_t primary_status = 83;
	Double_t secondary_status = 91;
	
	//Each vector is an event number of events analyzed is the total number of vectors
	int nEvents = ch1->GetEntries();
	
	cout<<"The number of events for this analysis is: "<<nEvents<<endl;
	
	//Definition of produced histograms
	//3D histograms
	TH3D *hPtrPaDEta = new TH3D("hDPtrPaDEta",Form("#Delta#phi p_{T} Triger p_{T} Associate for %s;p_{T} triger (Gev/c); p_{T} associate (GeV/c); #Delta#eta",title),100,0,50,100,0,50,80,-8,8);
	TH3D *hPtrPaDEtaPr = new TH3D("hDPtrPaDEtaPr",Form("#Delta#phi p_{T} Triger p_{T} Associate for %s from %.1f;p_{T} triger (Gev/c); p_{T} associate (GeV/c); #Delta#eta",title,primary_status),100,0,50,100,0,50,80,-8,8);
	TH3D *hPtrPaDEtaSc = new TH3D("hDPtrPaDEtaSc",Form("#Delta#phi p_{T} Triger p_{T} Associate for %s from %.1f;p_{T} triger (Gev/c); p_{T} associate (GeV/c); #Delta#eta",title,secondary_status),100,0,50,100,0,50,80,-8,8);
	
	//2D histograms
	TH2D *hDPhiDEta = new TH2D("hDPhiDEta",Form("%s #Delta#phi and #Delta#eta;#Delta#phi (rad);#Delta#eta",title),100,-PI/2,3*PI/2,80,-8,8);
	TH2D *hDPhiDEtaPr = new TH2D("hDPhiDEtaPr",Form("%s #Delta#phi and #Delta#eta from %.1f;#Delta#phi (rad);#Delta#eta",title,primary_status),100,-PI/2,3*PI/2,80,-8,8);
	TH2D *hDPhiDEtaSc = new TH2D("hDPhiDEtaSc",Form("%s #Delta#phi and #Delta#eta from %.1f;#Delta#phi (rad);#Delta#eta",title,secondary_status),100,-PI/2,3*PI/2,80,-8,8);
	
	//1D histogram accross the hole p_T
	TH1D *hTrPt = new TH1D("hTrPt",Form("Trigger Transverse Momentum for %s;p_{T} GeV/c;Counts",title),100,0,50);
	TH1D *hTrPtPr = new TH1D("hTrPtPr",Form("Trigger Transverse Momentum for %s from %.1f;p_{T} GeV/c;Counts",title,primary_status),100,0,50);
	TH1D *hTrPtSc = new TH1D("hTrPtSc",Form("Trigger Transverse Momentum for %s from %.1f;p_{T} GeV/c;Counts",title,secondary_status),100,0,50);
	TH1D *hStatusTr = new TH1D("hStatusTr",Form("Production mechanism of the trigger %s pair; Process ID; Counts",title),365,-182.5,182.5);
	TH1D *hStatusAs = new TH1D("hStatusAs",Form("Production mechanism of the associate %s pair; Process ID; Counts",title),365,-182.5,182.5);
	TH1D *hDPhi = new TH1D("hDPhi",Form("#Delta#phi for all processes for %s pair;#Delta#Phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D *hDPhiPr = new TH1D("hDPhiPr",Form("#Delta#phi for production mechanism %.1f for %s pair;#Delta#Phi (rad);Counts",primary_status,title),100,-PI/2,3*PI/2);
	TH1D *hDPhiSe = new TH1D("hDphiSe",Form("#Delta#phi for production mechanism %.1f for %s pair;#Delta#Phi (rad);Counts",secondary_status,title),100,-PI/2,3*PI/2);
	
	//no mechanism
	TH1D* hDPhiLL = new TH1D("hDPhiLL",Form("%s #Delta#phi correlation for low-low p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiIL = new TH1D("hDPhiIL",Form("%s #Delta#phi correlation for intermediate-low p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiII = new TH1D("hDPhiII",Form("%s #Delta#phi correlation for intermediate-intermediate p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiHL = new TH1D("hDPhiHL",Form("%s #Delta#phi correlation for high-low p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiHI = new TH1D("hDPhiHI",Form("%s #Delta#phi correlation for high-intermediate p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiHH = new TH1D("hDPhiHH",Form("%s #Delta#phi correlation for high-high p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	
	//primary mechanism
	TH1D* hDPhiLLPr = new TH1D("hDPhiLLPr",Form("%s #Delta#phi correlation for low-low p_{T} from %.1f;#Delta#phi (rad);Counts",title,primary_status),100,-PI/2,3*PI/2);
	TH1D* hDPhiILPr = new TH1D("hDPhiILPr",Form("%s #Delta#phi correlation for intermediate-low p_{T} from %.1f;#Delta#phi (rad);Counts",title, primary_status),100,-PI/2,3*PI/2);
	TH1D* hDPhiIIPr = new TH1D("hDPhiIIPr",Form("%s #Delta#phi correlation for intermediate-intermediate p_{T} from %.1f;#Delta#phi (rad);Counts",title, primary_status),100,-PI/2,3*PI/2);
	TH1D* hDPhiHLPr = new TH1D("hDPhiHLPr",Form("%s #Delta#phi correlation for high-low p_{T} from %.1f;#Delta#phi (rad);Counts",title, primary_status),100,-PI/2,3*PI/2);
	TH1D* hDPhiHIPr = new TH1D("hDPhiHIPr",Form("%s #Delta#phi correlation for high-intermediate p_{T} from %.1f;#Delta#phi (rad);Counts",title, primary_status),100,-PI/2,3*PI/2);
	TH1D* hDPhiHHPr = new TH1D("hDPhiHHPr",Form("%s #Delta#phi correlation for high-high p_{T} from %.1f;#Delta#phi (rad);Counts",title, primary_status),100,-PI/2,3*PI/2);
	
	//Secondary mechanism
	TH1D* hDPhiLLSc = new TH1D("hDPhiLLSc",Form("%s #Delta#phi correlation for low-low p_{T} from %.1f;#Delta#phi (rad);Counts",title, secondary_status),100,-PI/2,3*PI/2);
	TH1D* hDPhiILSc = new TH1D("hDPhiILSc",Form("%s #Delta#phi correlation for intermediate-low p_{T} from %.1f;#Delta#phi (rad);Counts",title, secondary_status),100,-PI/2,3*PI/2);
	TH1D* hDPhiIISc = new TH1D("hDPhiIISc",Form("%s #Delta#phi correlation for intermediate-intermediate p_{T} from %.1f;#Delta#phi (rad);Counts",title, secondary_status),100,-PI/2,3*PI/2);
	TH1D* hDPhiHLSc = new TH1D("hDPhiHLSc",Form("%s #Delta#phi correlation for high-low p_{T} from %.1f;#Delta#phi (rad);Counts",title, secondary_status),100,-PI/2,3*PI/2);
	TH1D* hDPhiHISc = new TH1D("hDPhiHISc",Form("%s #Delta#phi correlation for high-intermediate p_{T} from %.1f;#Delta#phi (rad);Counts",title, secondary_status),100,-PI/2,3*PI/2);
	TH1D* hDPhiHHSc = new TH1D("hDPhiHHSc",Form("%s #Delta#phi correlation for high-high p_{T} from %.1f;#Delta#phi (rad);Counts",title, secondary_status),100,-PI/2,3*PI/2);	
	//Event Loop
	for(int iEvent = 0; iEvent < nEvents; iEvent++){
		ch1->GetEntry(iEvent);
		int nparticles = vID->size();
		for(int ipart = 0; ipart < nparticles; ipart++){
			pID = (*vID)[ipart];
			pPhi = (*vPhi)[ipart];
			pPt = (*vPt)[ipart];
			pStatus = (*vStatus)[ipart];
			pEta =(*vEta)[ipart];
			if(pID == id_trigger){
				nTrigger++;
				hTrPt->Fill(pPt);
				hStatusTr->Fill(pStatus);
				
				for(int jpart = 0; jpart < nparticles; jpart++){
					if(jpart == ipart) continue;//Do not correlate with it self
						aID = (*vID)[jpart];
						aPhi = (*vPhi)[jpart];
						aPt = (*vPt)[jpart];
						aStatus = (*vStatus)[jpart];
						aEta = (*vEta)[jpart];
						if(aID == id_associate){
						hStatusAs->Fill(aStatus);
						hPtrPaDEta->Fill(pPt,aPt,pEta-aEta);
						hDPhiDEta->Fill(DeltaPhi(pPhi,aPhi),pEta-aEta);
						hDPhi->Fill(DeltaPhi(pPhi,aPhi));
						if( 81 <= pStatus && pStatus <= 89 &&  81 <= aStatus && aStatus <= 89 ){//Hadroniazation processes
							hTrPtPr->Fill(pPt);
							hDPhiPr->Fill(DeltaPhi(pPhi,aPhi));
							hPtrPaDEtaPr->Fill(pPt,aPt,pEta-aEta);
							hDPhiDEtaPr->Fill(DeltaPhi(pPhi,aPhi),pEta-aEta);
						
						}
						if(  91 <= pStatus && pStatus <= 99 &&  91 <= aStatus && aStatus <= 99 ){//Decay products
							hTrPtSc->Fill(pPt);
							hDPhiSe->Fill(DeltaPhi(pPhi,aPhi));
							hPtrPaDEtaSc->Fill(pPt,aPt,pEta-aEta);
							hDPhiDEtaSc->Fill(DeltaPhi(pPhi,aPhi),pEta-aEta);
						
						}
						//Filling triger momentum range correlations
							if(pPt >= 1. && pPt < 3. && aPt > 1. && aPt < 3.){
								hDPhiLL->Fill(DeltaPhi(pPhi,aPhi));
								if( 81 <= pStatus && pStatus <= 89 &&  81 <= aStatus && aStatus <= 89 ) hDPhiLLPr->Fill(DeltaPhi(pPhi,aPhi));
								if(  91 <= pStatus && pStatus <= 99 &&  91 <= aStatus && aStatus <= 99 ) hDPhiLLSc->Fill(DeltaPhi(pPhi,aPhi));
							}//End of low transverse momentum
							if(pPt >= 3. && pPt < 8. && aPt > 1. && aPt < 3.){
								hDPhiIL->Fill(DeltaPhi(pPhi,aPhi));
								if( 81 <= pStatus && pStatus <= 89 &&  81 <= aStatus && aStatus <= 89 ) hDPhiILPr->Fill(DeltaPhi(pPhi,aPhi));
								if(  91 <= pStatus && pStatus <= 99 &&  91 <= aStatus && aStatus <= 99 ) hDPhiILSc->Fill(DeltaPhi(pPhi,aPhi));
																
							}//End of intermidiate-low transverse momentum
							
							if(pPt >= 3. && pPt < 8. && aPt >= 3. && aPt < 8.){
								hDPhiII->Fill(DeltaPhi(pPhi,aPhi));
								if( 81 <= pStatus && pStatus <= 89 &&  81 <= aStatus && aStatus <= 89 ) hDPhiIIPr->Fill(DeltaPhi(pPhi,aPhi));
								if(  91 <= pStatus && pStatus <= 99 &&  91 <= aStatus && aStatus <= 99 ) hDPhiIISc->Fill(DeltaPhi(pPhi,aPhi));
								
							}
							
							if(pPt >= 8. && aPt > 1. && aPt < 3.){
								hDPhiHL->Fill(DeltaPhi(pPhi,aPhi));
								if( 81 <= pStatus && pStatus <= 89 &&  81 <= aStatus && aStatus <= 89 ) hDPhiHLPr->Fill(DeltaPhi(pPhi,aPhi));
								if(  91 <= pStatus && pStatus <= 99 &&  91 <= aStatus && aStatus <= 99 ) hDPhiHLSc->Fill(DeltaPhi(pPhi,aPhi));
								
							}//End of high-low transverse momentum
							
							if(pPt >= 8. && aPt >=3. && aPt < 8.){
								hDPhiHI->Fill(DeltaPhi(pPhi,aPhi));
								if( 81 <= pStatus && pStatus <= 89 &&  81 <= aStatus && aStatus <= 89 ) hDPhiHIPr->Fill(DeltaPhi(pPhi,aPhi));
								if(  91 <= pStatus && pStatus <= 99 &&  91 <= aStatus && aStatus <= 99 ) hDPhiHISc->Fill(DeltaPhi(pPhi,aPhi));
								
							}//End of high-intermediate transverse momentum
							
							if(pPt >= 8. && aPt >= 8.){
								hDPhiHH->Fill(DeltaPhi(pPhi,aPhi));
								if( 81 <= pStatus && pStatus <= 89 &&  81 <= aStatus && aStatus <= 89 ) hDPhiHHPr->Fill(DeltaPhi(pPhi,aPhi));
								if(  91 <= pStatus && pStatus <= 99 &&  91 <= aStatus && aStatus <= 99 ) hDPhiHHSc->Fill(DeltaPhi(pPhi,aPhi));
							
							}//End of high-high transverse momentum
						}//Associate Codition	
				}//Associate Loop
			}//Trigger Codition
		}//Trigger Loop
	}//End of event loop
	
	if(nTrigger == 0){
		cout<<"Have not found any trigger particle with id: "<<id_trigger<<endl;
		output->Close();
		return 0;
	}
	output->Write();
	output->Close();
	cout<<"The total number of triggers is: "<<nTrigger<<endl;
	
	cout<<"File: "<<filename<<" has been created!"<<endl;
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

void normalize(TH1D *hist, TH1D *trig, Double_t xmin, Double_t xmax){
	Int_t bin_start = trig->FindBin(xmin);
	Int_t bin_finish = trig->FindBin(xmax);
	Double_t integral = trig->Integral(bin_start, bin_finish);
	hist->Scale(1./integral);
	hist->GetYaxis()->SetTitle("#frac{1}{N_{tr}} #frac{dN}{d#Delta#phi}");
}

void plot_all(TString filename){
	//Dictionaries used
	const char* DPhi_Names[6] = {"hDPhiLL","hDPhiIL","hDPhiII","hDPhiHL","hDPhiHI","hDPhiHH"};
	const char* DPhi_NamesPr[6] = {"hDPhiLLPr","hDPhiILPr","hDPhiIIPr","hDPhiHLPr","hDPhiHIPr","hDPhiHHPr"};
	const char* DPhi_NamesSc[6] = {"hDPhiLLSc","hDPhiILSc","hDPhiIISc","hDPhiHLSc","hDPhiHISc","hDPhiHHSc"};
	const char* Triggers_Names[3] = {"hTrPt","hTrPtPr","hTrPtSc"};
	const char* DPhiDEta_Names[3] = {"hDPhiDEta","hDPhiDEtaPr","hDPhiDEtaSc"};
	
	//Histogram arrays
	TH2D* hDPhiDEta[3];
	TH1D* hTriggers[3];
	TH1D* hDPhi[6];
	TH1D* hDPhiPr[6];
	TH1D* hDPhiSc[6];
	
	//Reading histograms
	for(int i = 0; i<6;i++){
		hDPhi[i] = Read_Hist(filename,DPhi_Names[i]);
		hDPhiPr[i] = Read_Hist(filename,DPhi_NamesPr[i]);
		hDPhiSc[i] = Read_Hist(filename,DPhi_NamesSc[i]);
		if( i < 3 ){
			hDPhiDEta[i] =(TH2D*) Read_Hist(filename,DPhiDEta_Names[i]);// Castind needed because function returns TH1D
			hTriggers[i] = Read_Hist(filename,Triggers_Names[i]);
		}
	}
	
	//Normalizing 1D azimuthial plots
	for (int i = 0; i < 6; i++){
		if(i == 0){
			normalize(hDPhi[i],hTriggers[0],1.,3.);
			normalize(hDPhiPr[i],hTriggers[0],1.,3.);
			normalize(hDPhiSc[i],hTriggers[0],1.,3.);
		}
		if(i == 1 || i == 2){
			normalize(hDPhi[i],hTriggers[0],3.,8.);
			normalize(hDPhiPr[i],hTriggers[0],3.,8.);
			normalize(hDPhiSc[i],hTriggers[0],3.,8.);
			
		}
		if(i == 3 || i == 4 || i == 5){
			normalize(hDPhi[i],hTriggers[0],8.,50.);
			normalize(hDPhiPr[i],hTriggers[0],8.,50.);
			normalize(hDPhiSc[i],hTriggers[0],8.,50.);
		}
	}
	
	//Normalizing 2D azimuthial plot
	for(int i =0;i<3;i++){
		hDPhiDEta[i]->Scale(1./(hTriggers[0]->Integral(hTriggers[0]->FindBin(0.15),hTriggers[0]->FindBin(50.))));
		hDPhiDEta[i]->GetZaxis()->SetTitle("#frac{1}{N_{tr}} #frac{dN}{#Delta#phi#Delta#eta}");
	}
	//Drawing
	//Setting up canvases and legends
	TCanvas *c2D[3];
	TCanvas *c1D[6];
	
	//Drawing 2D histograms
	gStyle->SetPalette(kRainBow);
	for(int i = 0;i<3;i++){
	hDPhiDEta[i]->GetZaxis()->SetTitleOffset(1.1);
	hDPhiDEta[i]->SetStats(0);
	hDPhiDEta[i]->SetContour(80);
	c2D[i] = new TCanvas();
	hDPhiDEta[i]->Draw("SURF1");
	}
	//Drawing 1D histograms
	gStyle->SetPadLeftMargin(2); 
	gStyle->SetPadRightMargin(0.01);
	for(int i = 0; i<6;i++){
		c1D[i] = new TCanvas();
		c1D[i]->Divide(3,1);
		hDPhi[i]->SetMarkerStyle(20);
		hDPhi[i]->SetMarkerColor(2);
		hDPhiPr[i]->SetMarkerStyle(20);
		hDPhiPr[i]->SetMarkerColor(3);
		hDPhiSc[i]->SetMarkerStyle(20);
		hDPhiSc[i]->SetMarkerColor(4);
		hDPhi[i]->SetStats(0);
		hDPhiPr[i]->SetStats(0);
		hDPhiSc[i]->SetStats(0);
		c1D[i]->cd(1);
		hDPhi[i]->Draw("PE1");
		c1D[i]->cd(2);
		hDPhiPr[i]->Draw("PE1");
		c1D[i]->cd(3);
		hDPhiSc[i]->Draw("PE1");
	}
}

TGraphErrors* plot_ratios(TString filename){
	//This function returns a plot where the ratio of string fragmentation pairs over the total pairs is shown 
	//as a function of p_T
	
	//Dictionaries used
	const char* DPhi_Names[6] = {"hDPhiLL","hDPhiIL","hDPhiII","hDPhiHL","hDPhiHI","hDPhiHH"};
	const char* DPhi_NamesPr[6] = {"hDPhiLLPr","hDPhiILPr","hDPhiIIPr","hDPhiHLPr","hDPhiHIPr","hDPhiHHPr"};
	//Arrays used
	Double_t x[3];
	Double_t entries_total[3];
	Double_t entries_string[3];
	TH1D* hDPhi[6];
	TH1D* hDPhiPr[6];
	TH1D* hTrigPt;
	
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
	entries_total[0] = hDPhi[0]->GetEntries();
	entries_string[0] = hDPhiPr[0]->GetEntries();
	
	entries_total[1] = hDPhi[1]->GetEntries()+hDPhi[2]->GetEntries();
	entries_string[1] = hDPhiPr[1]->GetEntries()+hDPhiPr[1]->GetEntries();
	
	entries_total[2] = hDPhi[3]->GetEntries()+hDPhi[4]->GetEntries()+hDPhi[5]->GetEntries();
	entries_string[2] = hDPhiPr[3]->GetEntries()+hDPhiPr[4]->GetEntries()+hDPhiPr[5]->GetEntries();
	
	//Allocating graph
	TGraphErrors *gr = new TGraphErrors();
	
	//Filling graph.
	Int_t n = 0;
	for(int i = 0; i<3; i++){
		n = gr->GetN();
		gr->SetPoint(n,x[i],entries_string[i]/entries_total[i]);
		gr->SetPointError(n,0,0);
	}
	
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(2.);
	gr->SetMarkerColor(1);
	gr->GetXaxis()->SetTitle("#it{p_{T}^{tr.}} (GeV/c)");
	gr->GetYaxis()->SetTitle("Ratio");
	return gr;
	
	
}

//+-+-+-++-+-+-+-+-+-+-+-+++-++-+-+-+-+-+-+-+-+-+-+-+-+-
//Here I include functions that are used for background reduction
//Yield and width calculation.
//ATTENTION !!!!!
//The code is not the same as in yield_width_calculation.cpp things needed to be optimized
//for the production mechanism related analysis.

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

TH1D *reduction(TH1D *hist, Double_t start, Double_t finish){
//This function does the whole procces
	double *p = Zyam_and_Error(start,finish,hist);
	TH1D *Signal = signal(hist,*p,*(p+1));
	cout<<"The ZYAM background is: "<<*p<<endl;

	delete p;
	return Signal; 
}

void background_plot(TString filename, Double_t start[6], Double_t finish[6]){
	//This function plots everything across all momentum ranges
	//Dictionaries
	const char* DPhi_Names[6] = {"hDPhiLL","hDPhiIL","hDPhiII","hDPhiHL","hDPhiHI","hDPhiHH"};
	
	
	//Allocation of arrays
	TH1D* hDPhi[6];
	TLine* line[6];
	TH1D* reduced[6];
	TH1D* hTrigPt = Read_Hist(filename,"hTrPt");
	
	//Importing histograms
	for(int i = 0; i < 6; i++){
		hDPhi[i] = Read_Hist(filename, DPhi_Names[i]);
	}
	//Nomralizing to the number of triggers
	for(int i = 0; i < 6; i++){
		if(i == 0) normalize(hDPhi[i],hTrigPt,1.,3.);
		if(i == 1 || i == 2) normalize(hDPhi[i],hTrigPt,3.,8.);
		if(i == 3 || i == 4 || i == 5) normalize(hDPhi[i],hTrigPt,8.,50.);
		hDPhi[i]->SetStats(0);
	}
	//Calculation ZYAM and error and line creation
	for(int i = 0; i<6; i++){
		double *pZYAM = Zyam_and_Error(start[i],finish[i],hDPhi[i]);
		line[i] = background_line(*(pZYAM),hDPhi[i]);
		reduced[i] = reduction(hDPhi[i],start[i],finish[i]);
		delete pZYAM;
	}
	//Drawing everything
	TLegend *leg = new TLegend();
	leg->AddEntry(hDPhi[0], "D^{+}D^{+} Correlation","P");
	leg->AddEntry(line[0],"ZYAM background","L");
	
	gStyle->SetPadLeftMargin(0.18); gStyle->SetPadRightMargin(0.01);
	TCanvas *c1 = new TCanvas();
	c1->Divide(3,2);
	TCanvas *c2 = new TCanvas();
	c2->Divide(3,2);
	
	for(int i = 0; i < 6;i++){
		c1->cd(i+1);
		hDPhi[i]->Draw("P");
		line[i]->Draw();
		leg->Draw();
		c2->cd(i+1);
		reduced[i]->Draw();
	}
}
//At this point I am introducing the yield and width calculation in all functions
//There will be the int mechanism argument that is implemented to the user
//can choose if he or she wants to do the computation for all the correlations
//or specific to correlations created from string fragmentation or decay 
//specifically. This can be done by commenting in and out the DPhiNames 
//arrays.A develpment to tune this better is on the way!


TGraphErrors *Yield_Low(TString filename,Double_t start[3], Double_t finish[3],Double_t low_limit, Double_t high_limit){
	//This function takes the filename the start and finish for the background calculatiion and the limits for the yield calculation
	//and returns a graph where the per trigger yield is displayed when associating with low momentum particles.
	//Dictionaries used
	
	
	const char* DPhi_Names[3] = {"hDPhiLL","hDPhiIL","hDPhiHL"};
	//const char* DPhi_Names[3]= {"hDPhiLLPr","hDPhiILPr","hDPhiHLPr"};
	//const char* DPhi_Names[3] = {"hDPhiLLSc","hDPhiILSc","hDPhiHLSc"};
	
	
	
	//Arrays used
	TH1D* hDPhi[3];
	Double_t x[3];
	Double_t y[3];
	Double_t xerr[3];
	Double_t yerr[3];
	TH1D* reduced[3];
	
	//Reading files
	for(int i = 0; i < 3; i++){
		hDPhi[i] = Read_Hist(filename, DPhi_Names[i]);
	}
	
	TH1D *hTrigPt = Read_Hist(filename,"hTrPt");
	
	//Normalizing and background reducing
	for(int i = 0; i < 3; i++){
		if(i == 0) normalize(hDPhi[i],hTrigPt,1.,3.);
		if(i == 1) normalize(hDPhi[i],hTrigPt,3.,8.);
		if(i == 2) normalize(hDPhi[i],hTrigPt,8.,50.);
		hDPhi[i]->SetStats(0);
		reduced[i] = reduction(hDPhi[i],start[i],finish[i]);
	}
	
	//X-axis input
	//Rescaling is needed in order to 
	//be corretly plotted with the associate p_t axis
	hTrigPt->GetXaxis()->SetRangeUser(1.,3.);
	Double_t p_assoc_l = hTrigPt->GetMean();
	
	hTrigPt->GetXaxis()->SetRangeUser(1.,3.);
	x[0] = hTrigPt->GetMean();
	xerr[0] = hTrigPt->GetStdDev();
	
	hTrigPt->GetXaxis()->SetRangeUser(3.,8.);
	x[1] = hTrigPt->GetMean();
	xerr[1] = hTrigPt->GetStdDev();
	x[1] = ((x[1]-3.)/2.2)*p_assoc_l+3.;
	
	hTrigPt->GetXaxis()->SetRangeUser(8.,50.);
	x[2] = hTrigPt->GetMean();
	xerr[2] = hTrigPt->GetStdDev();
	x[2] = ((x[2]-8.)/6.)*p_assoc_l+8.;
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
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.5);
	gr->SetMarkerColor(1);
	gr->GetXaxis()->SetTitle("#it{p_{T}^{tr.}} (GeV/c)");
	gr->GetYaxis()->SetTitle("Yield");
	gr->GetXaxis()->SetLimits(0.15,16.1);
	return gr;
}

TGraphErrors *Yield_Intermediate(TString filename,Double_t start[2], Double_t finish[2],Double_t low_limit, Double_t high_limit){
	//This function takes the filename the start and finish for the background calculatiion and the limits for the yield calculation
	//and returns a graph where the per trigger yield is displayed when associating with low momentum particles.
	//Dictionaries used
	
	
	const char* DPhi_Names[2] = {"hDPhiII","hDPhiHI"};
	//const char* DPhi_Names[2]= {"hDPhiIIPr","hDPhiHIPr"};
	//const char* DPhi_Names[2] = {"hDPhiIISc","hDPhiHISc"};
	
	
	
	//Arrays used
	TH1D* hDPhi[2];
	Double_t x[2];
	Double_t y[2];
	Double_t xerr[2];
	Double_t yerr[2];
	TH1D* reduced[2];
	
	//Reading files
	for(int i = 0; i < 2; i++){
		hDPhi[i] = Read_Hist(filename, DPhi_Names[i]);
	}
	
	TH1D *hTrigPt = Read_Hist(filename,"hTrPt");
	
	//Normalizing and background reducing
	for(int i = 0; i < 2; i++){
		if(i == 0) normalize(hDPhi[i],hTrigPt,3.,8.);
		if(i == 1) normalize(hDPhi[i],hTrigPt,8.,50.);
		hDPhi[i]->SetStats(0);
		reduced[i] = reduction(hDPhi[i],start[i],finish[i]);
	}
	
	//X-axis input
	//Rescaling is needed in order to 
	//be corretly plotted with the associate p_t axis
	hTrigPt->GetXaxis()->SetRangeUser(3.,8.);
	Double_t p_assoc_l = hTrigPt->GetMean();
	
	hTrigPt->GetXaxis()->SetRangeUser(3.,8.);
	x[0] = hTrigPt->GetMean();
	xerr[0] = hTrigPt->GetStdDev();
	x[0] = ((x[0]-3.)/2.2)*p_assoc_l+3.;
	
	hTrigPt->GetXaxis()->SetRangeUser(8.,50.);
	x[1] = hTrigPt->GetMean();
	xerr[1] = hTrigPt->GetStdDev();
	x[1] = ((x[1]-8.)/6.)*p_assoc_l+8.;
	//Input for y axis
	int low_bin = reduced[0]->FindBin(low_limit);
	int high_bin = reduced[0]->FindBin(high_limit);
	for(int i = 0; i<2; i++){
		y[i] = reduced[i]->IntegralAndError(low_bin,high_limit,yerr[i],"");
	}
	
	//Setting up and filling graph
	TGraphErrors *gr = new TGraphErrors();
	Int_t n = 0;
	cout<<"For low associate momentum we have:"<<endl;
	for(int i = 0; i<2; i++){
		cout<<"Filling x = "<<x[i]<<" +- "<<xerr[i]<<" and y = "<<y[i]<<" +- "<<yerr[i]<<endl;
		n = gr->GetN();
		gr->SetPoint(n,x[i],y[i]);
		gr->SetPointError(n,xerr[i],yerr[i]);
	}
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.5);
	gr->SetMarkerColor(1);
	gr->GetXaxis()->SetTitle("#it{p_{T}^{tr.}} (GeV/c)");
	gr->GetYaxis()->SetTitle("Yield");
	gr->GetXaxis()->SetLimits(0.15,16.1);
	return gr;
}

TGraphErrors *Yield_High(TString filename,Double_t start, Double_t finish,Double_t low_limit, Double_t high_limit){
	//This functions takes the filename the start and finish for the background calculatiion and the limits for the yield calculation
	//and returns a graph where the per trigger yield is displayed when associating with intermediate momentum particles.
	
	

	//Reading files
	
	TH1D *hDPhi = Read_Hist(filename, "hDPhiHH");
	//TH1D *hDPhi = Read_Hist(filename, "hDPhiHHPr");
	//TH1D *hDPhi = Read_Hist(filename, "hDPhiHHSc");
	TH1D *hTrigPt = Read_Hist(filename,"hTrPt");
	
	//Normalizing and background reducing
	
	normalize(hDPhi,hTrigPt,8.,50.);
	TH1D *reduced = reduction(hDPhi,start,finish);
	
	
	//X-axis input
	
	
	Double_t x,xerr,y,yerr;
	
	hTrigPt->GetXaxis()->SetRangeUser(8.,50.);
	x = hTrigPt->GetMean();
	xerr = hTrigPt->GetStdDev();
	x = ((x-8.)/6.)*x+8.;
	
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
	
	gr->SetMarkerStyle(20.);
	gr->SetMarkerSize(1.5);
	gr->SetMarkerColor(1);
	gr->GetXaxis()->SetTitle("#it{p_{T}^{tr.}} (GeV/c)");
	gr->GetYaxis()->SetTitle("Yield");
	gr->GetXaxis()->SetLimits(0.15,16.1);
	return gr;
}

//End of yield calculations the width calculations are following

TF1 *fitting(TH1D* hist,Double_t start,Double_t finish){
	// This finction does the fitting of the histogram
	//hist for the given range start finish
	TF1 *fit = new TF1("fit","gaus",start,finish);
	fit->SetLineWidth(1);
	fit->SetLineColor(2);
	fit->SetLineStyle(2);
	hist->Fit("fit","R");
	cout<< "The width for the range "<<start<<" to "<<finish<<" is: "<<fit->GetParameter(2)<< " +- "<<fit->GetParError(2)<<endl;
	cout<<"The chi^2/d.o.f is: "<<fit->GetChisquare()/((hist->FindBin(finish)-hist->FindBin(start))-3)<<endl;
	return fit;
}


void fit_and_plot_all(TString filename,Double_t start[6],Double_t finish[6],Double_t fit_start, Double_t fit_finish){
	//With this function one can plot the fitted histograms all together
	//Dictionaries used
	const char* DPhi_Names[6] = {"hDPhiLL","hDPhiIL","hDPhiII","hDPhiHL","hDPhiHI","hDPhiHH"};
	
	//Allocation of arrays
	TH1D* hDPhi[6];
	TH1D*red[6];
	TF1* fits[6];
	TLatex* text[6];
	TH1D *hTrigPt = Read_Hist(filename,"hTrPt");
	double fits_chi[6];
	double fits_width[6];
	double fits_width_err[6];
	
	//Importing histograms
	for(int i = 0; i < 6; i++){
		hDPhi[i] = Read_Hist(filename, DPhi_Names[i]);
	}
	
	//Normalizing
	for(int i = 0; i < 6; i++){
		if(i == 0) normalize(hDPhi[i],hTrigPt,1.,3.);
		if(i == 1 || i == 2) normalize(hDPhi[i],hTrigPt,3.,8.);
		if(i == 3 || i == 4 || i == 5) normalize(hDPhi[i],hTrigPt,8.,50.);
		hDPhi[i]->SetStats(0);
	}
	
	//Reducing background and fitting
	for(int i = 0; i<6; i++){
		double *pZYAM = Zyam_and_Error(start[i],finish[i],hDPhi[i]);//Dynamic allocation inside this function
		red[i] = reduction(hDPhi[i],start[i],finish[i]);
		red[i]->SetStats(0);
		fits[i] = fitting(red[i],fit_start,fit_finish);
		//Storing fit parameters
		fits_chi[i] = fits[i]->GetChisquare()/((red[i]->FindBin(fit_finish)-red[i]->FindBin(fit_start))-3);
		fits_width[i] = fits[i]->GetParameter(2);
		fits_width_err[i] = fits[i]->GetParError(2);
		//Showing results in plot
		text[i] = new TLatex(-1.,0.000045,Form("#scale[0.6]{#splitline{Width = %f #pm %f}{#chi^{2}/d.o.f. = %f}}",fits_width[i],fits_width_err[i],fits_chi[i]));
		delete pZYAM; //to avoid meemory leaks
	}
	TLegend* leg = new TLegend(0.89,0.8,0.99,0.99);
	leg->AddEntry(red[2],"ZYAM reduced correlation plot","P");
	leg->AddEntry(fits[2],"Fitting","L");
	leg->SetBorderSize(0);
	//Plotting
	TCanvas *mycanvas = new TCanvas();
	mycanvas->Divide(3,2);
	for(int i = 0; i < 6; i++){
		mycanvas->cd(i+1);
		red[i]->Draw();
		if(i == 2) leg->Draw();
		text[i]->Draw();
	}	
}


//For width plotting
TGraphErrors *Width_Low(TString filename,Double_t start[3], Double_t finish[3],Double_t low_limit, Double_t high_limit){
const char* DPhi_Names[3] = {"hDPhiLL","hDPhiIL","hDPhiHL"};
	//const char* DPhi_Names[3]= {"hDPhiLLPr","hDPhiILPr","hDPhiHLPr"};
	//const char* DPhi_Names[3] = {"hDPhiLLSc","hDPhiILSc","hDPhiHLSc"};
	
	
	
	//Arrays used
	TH1D* hDPhi[3];
	TF1* fits[3];
	Double_t x[3];
	Double_t y[3];
	Double_t xerr[3];
	Double_t yerr[3];
	TH1D* reduced[3];
	
	//Reading files
	for(int i = 0; i < 3; i++){
		hDPhi[i] = Read_Hist(filename, DPhi_Names[i]);
	}
	
	TH1D *hTrigPt = Read_Hist(filename,"hTrPt");
	
	//Normalizing and background reducing
	for(int i = 0; i < 3; i++){
		if(i == 0) normalize(hDPhi[i],hTrigPt,1.,3.);
		if(i == 1) normalize(hDPhi[i],hTrigPt,3.,8.);
		if(i == 2) normalize(hDPhi[i],hTrigPt,8.,50.);
		hDPhi[i]->SetStats(0);
		reduced[i] = reduction(hDPhi[i],start[i],finish[i]);
	}
	
	//X-axis input
	hTrigPt->GetXaxis()->SetRangeUser(1.,3.);
	Double_t p_assoc_l = hTrigPt->GetMean();
	
	hTrigPt->GetXaxis()->SetRangeUser(1.,3.);
	x[0] = hTrigPt->GetMean();
	xerr[0] = hTrigPt->GetStdDev();
	
	hTrigPt->GetXaxis()->SetRangeUser(3.,8.);
	x[1] = hTrigPt->GetMean();
	xerr[1] = hTrigPt->GetStdDev();
	x[1] = ((x[1]-3.)/2.2)*p_assoc_l+3.;
	
	hTrigPt->GetXaxis()->SetRangeUser(8.,50.);
	x[2] = hTrigPt->GetMean();
	xerr[2] = hTrigPt->GetStdDev();
	x[2] = ((x[2]-8.)/6.)*p_assoc_l+8.;
	//Input for y axis
	for(int i = 0; i<3; i++){
		fits[i] = fitting(reduced[i],low_limit,high_limit);
		y[i] = fits[i]->GetParameter(2);
		yerr[i] = fits[i]->GetParError(2);
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
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.5);
	gr->SetMarkerColor(1);
	gr->GetXaxis()->SetTitle("#it{p_{T}^{tr.}} (GeV/c)");
	gr->GetYaxis()->SetTitle("Width");
	gr->GetXaxis()->SetLimits(0.15,16.1);
	return gr;
}

TGraphErrors *Width_Intermediate(TString filename,Double_t start[2], Double_t finish[2],Double_t low_limit, Double_t high_limit){
	//This function take the filename the start and finish for the background calculatiion and the limits for the yield calculation
	//and returns a graph where the per trigger yield is displayed when associating with intermediate momentum particles.
	//Dictionaries used
	const char* DPhi_Names[2] = {"hDPhiII","hDPhiHI"};
	//const char* DPhi_Names[2]= {"hDPhiIIPr","hDPhiHIPr"};
	//const char* DPhi_Names[2] = {"hDPhiIISc","hDPhiHISc"};
	
	
	
	//Arrays used
	TH1D* hDPhi[2];
	TF1* fits[2];
	Double_t x[2];
	Double_t y[2];
	Double_t xerr[2];
	Double_t yerr[2];
	TH1D* reduced[2];
	
	//Reading files
	for(int i = 0; i < 2; i++){
		hDPhi[i] = Read_Hist(filename, DPhi_Names[i]);
	}
	
	TH1D *hTrigPt = Read_Hist(filename,"hTrPt");
	
	//Normalizing and background reducing
	for(int i = 0; i < 2; i++){
		if(i == 0) normalize(hDPhi[i],hTrigPt,3.,8.);
		if(i == 1) normalize(hDPhi[i],hTrigPt,8.,50.);
		hDPhi[i]->SetStats(0);
		reduced[i] = reduction(hDPhi[i],start[i],finish[i]);
	}
	
	//X-axis input
	hTrigPt->GetXaxis()->SetRangeUser(3.,8.);
	Double_t p_assoc_i = hTrigPt->GetMean();
	
	hTrigPt->GetXaxis()->SetRangeUser(3.,8.);
	x[0] = hTrigPt->GetMean();
	xerr[0] = hTrigPt->GetStdDev();
	x[0] = ((x[0]-3.)/2.2)*p_assoc_i+3.;
	
	hTrigPt->GetXaxis()->SetRangeUser(8.,50.);
	x[1] = hTrigPt->GetMean();
	xerr[1] = hTrigPt->GetStdDev();
	x[1] = ((x[1]-8.)/6.)*p_assoc_i+8.;
	
	
	//Input for y axis
	for(int i = 0; i<2; i++){
		fits[i] = fitting(reduced[i],low_limit,high_limit);
		y[i] = fits[i]->GetParameter(2);
		yerr[i] = fits[i]->GetParError(2);
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
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.5);
	gr->SetMarkerColor(1);
	gr->GetXaxis()->SetTitle("#it{p_{T}^{tr.}} (GeV/c)");
	gr->GetYaxis()->SetTitle("Width");
	gr->GetXaxis()->SetLimits(0.15,16.1);
	return gr;
}

TGraphErrors *Width_High(TString filename,Double_t start, Double_t finish,Double_t low_limit, Double_t high_limit){
	//This functions takes the filename the start and finish for the background calculatiion and the limits for the yield calculation
	//and returns a graph where the per trigger yield is displayed when associating with intermediate momentum particles.
	
	

	//Reading files
	
	//Reading files
	
	TH1D *hDPhi = Read_Hist(filename, "hDPhiHH");
	//TH1D *hDPhi = Read_Hist(filename, "hDPhiHHPr");
	//TH1D *hDPhi = Read_Hist(filename, "hDPhiHHSc");
	TH1D *hTrigPt = Read_Hist(filename,"hTrPt");
	
	//Normalizing and background reducing
	
	normalize(hDPhi,hTrigPt,8.,50.);
	TH1D *reduced = reduction(hDPhi,start,finish);
	
	
	//X-axis input
	
	
	Double_t x,xerr,y,yerr;
	
	hTrigPt->GetXaxis()->SetRangeUser(8.,50.);
	x = hTrigPt->GetMean();
	xerr = hTrigPt->GetStdDev();
	x = ((x-8.)/6.)*x+8.;
	
	//Input for y axis
	TF1* fit = fitting(reduced,low_limit,high_limit);
	y = fit->GetParameter(2);
	yerr = fit->GetParError(2);
	
	
	//Setting up and filling graph
	TGraphErrors *gr = new TGraphErrors();
	Int_t n = 0;
	cout<<"For high associate momentum we have:"<<endl;

	cout<<"Filling x = "<<x<<" +- "<<xerr<<" and y = "<<y<<" +- "<<yerr<<endl;
	n = gr->GetN();
	gr->SetPoint(n,x,y);
	gr->SetPointError(n,xerr,yerr);
	
	gr->SetMarkerStyle(20.);
	gr->SetMarkerSize(1.5);
	gr->SetMarkerColor(1);
	gr->GetXaxis()->SetTitle("#it{p_{T}^{tr.}} (GeV/c)");
	gr->GetYaxis()->SetTitle("Width");
	gr->GetXaxis()->SetLimits(0.15,16.1);
	return gr;
}

void status_analysis_v2(){
	//status_file(411,411,"DD_correlation_status.root","D^{+}D^{+}");
	/*
	plot_all("complete_root/DplusDminus_complete.root");
	TGraphErrors* gr1 = plot_ratios("complete_root/DplusDplus_complete.root");
	TGraphErrors* gr2 = plot_ratios("complete_root/DplusDminus_complete.root");
	TLegend *leg = new TLegend(0.8,0.7,0.9,0.9);
	leg->AddEntry(gr1,"D^{+}D^{+}","P");
	leg->AddEntry(gr2,"D^{+}D^{-}","P");
	leg->SetBorderSize(0);
	gr2->SetMarkerColor(3);
	TCanvas *c1 = new TCanvas();
	gr1->Draw("AP");
	gr2->Draw("SAME P");
	leg->Draw();
	*/
	
	//Arrays with calculation limits these can be changed from the user.
	Double_t start[6] = {-0.4,-0.3,-0.1,-0.8,-1,-1};
	Double_t finish[6] = {0.4,0.4,0.1,0.8,1,1};
	
	Double_t start_l[3] = {-0.4,-0.3,-0.8};
	Double_t finish_l[3] = {0.4,0.4,0.8};

	Double_t start_i[2] = {-0.1,-1.};
	Double_t finish_i[2] = {0.1,1.};
	
	Double_t start_h = -1;
	Double_t finish_h = 1;

	
	
	//background_plot("complete_root/DplusDplus_complete.root",start,finish);
	//fit_and_plot_all("complete_root/DplusDplus_complete.root",start,finish,PI/2,3*PI/2);
	TGraphErrors *gr = Yield_Low("complete_root/DplusDplus_complete.root",start_l,finish_l,PI/2,3*PI/2);
	TGraphErrors *gr_i = Yield_Intermediate("complete_root/DplusDplus_complete.root",start_i,finish_i,PI/2,3*PI/2);
	TGraphErrors *gr_h = Yield_High("complete_root/DplusDplus_complete.root",start_h,finish_h,PI/2,3*PI/2);
	TGraphErrors *grw = Width_Low("complete_root/DplusDplus_complete.root",start_l,finish_l,PI/2,3*PI/2);
	TGraphErrors *gr_iw = Width_Intermediate("complete_root/DplusDplus_complete.root",start_i,finish_i,PI/2,3*PI/2);
	TGraphErrors *gr_hw = Width_High("complete_root/DplusDplus_complete.root",start_h,finish_h,PI/2,3*PI/2);
	
	//Plotting
	TGaxis *low = new TGaxis(8.,0.0051,16.1,0.0051,0.15,16.1,10,"-L");
	low->SetLabelColor(kRed+2);
 	low->SetLineColor(kRed+2);
  	low->SetTitleColor(kRed+2);
  	low->CenterTitle();
  	low->SetLabelSize(0.03);
  	low->SetTitleSize(0.03);
  	//low->SetTitle("#it{p_{T}^{tr}} (GeV/c)");
  	
  	TGaxis *inter = new TGaxis(3.,0.0051,8.,0.0051,0.15,8.,10,"-L");
	inter->SetLabelColor(kRed+2);
 	inter->SetLineColor(kRed+2);
  	inter->SetTitleColor(kRed+2);
  	inter->CenterTitle();
  	inter->SetLabelSize(0.03);
  	inter->SetTitleSize(0.03);
  	inter->SetTitle("#it{p_{T}^{as.}} (GeV/c)");
  	
  	TGaxis *high = new TGaxis(0.15,0.0051,3.,0.0051,0.15,3.1,10,"-L");
	high->SetLabelColor(kRed+2);
 	high->SetLineColor(kRed+2);
  	high->SetTitleColor(kRed+2);
  	high->CenterTitle();
  	high->SetLabelSize(0.03);
  	high->SetTitleSize(0.03);
  	//high->SetTitle("#it{p_{T}^{tr}} (GeV/c)");
  	
  	TLine *line_one = new TLine(3.,0.001,3.,0.0051);
  	line_one->SetLineStyle(2);
  	line_one->SetLineWidth(1);
  	
  	TLine *line_two = new TLine(8.,0.001,8.,0.0051);
  	line_two->SetLineStyle(2);
  	line_two->SetLineWidth(1);
  	
  	TCanvas *c1 = new TCanvas();
	gr->GetYaxis()->SetRangeUser(0.001,0.0051);
	gr->GetXaxis()->SetLimits(0.15,16.1);
	gr->Draw("A P");
	high->Draw();
	inter->Draw();
	low->Draw();
	gr_i->Draw("SAME P");
	gr_h->Draw("SAME P");
	line_one->Draw();
	line_two->Draw();
	
	return 0;
}



