//C++ libraries
#include <iostream>
#include <vector>
#include <cstring>
//Root Libraries
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"

#define PI 3.14159265
using namespace std;
using namespace TMath;

Double_t DeltaPhi(Double_t phi1, Double_t phi2){
	//returns delta phi in range -pi/2 3pi/2
	return fmod(phi1-phi2+2.5*PI,2*PI)-0.5*PI;
	}

void correlation_plot(Int_t id_trigger,Int_t id_associate, TString filename, const char* title){
	//This function takes the pdg ids of the triger and associate particle as inputs and produces all kinds of
	//correlation plots concerning the DeltaPhi and DeltaEta as a function of pT.
	//filename: filename of root file that the histograms are going to be stored
	//title: Part of the histograms title concerning the particle for example it can be something like "D^{+}D^{+}"
	
	//Define the TChain
	TChain *ch1 = new TChain("tree");
	TFile *output = new TFile(filename,"RECREATE");
	
	int ntrees = 3; //Number of trees we want to add to the TChain
	
	//We put the trees to the chain
	for( int i = 1; i < ntrees+1;  i++){
	//One can change the file path accordingly to his/hers machine set up
	ch1->Add(Form("/home/isidiras/university_staff/ccbar_MONASH/MONASH_BATCH1/Group%i/output.root",i));
	}
	
	//Now we define the vectors and set branch addresses remember each vector corresponds to an event
	vector<Int_t>* vID = 0;
	vector<Double_t>* vPt = 0;//Carefull of null pointers.
	vector<Double_t>* vEta = 0;
	vector<Double_t>* vPhi = 0;
	vector<Double_t>* vCharge = 0;
	//Setting up chain addresses
	ch1->SetBranchAddress("ID",&vID);
	ch1->SetBranchAddress("PT",&vPt);
	ch1->SetBranchAddress("ETA",&vEta);
	ch1->SetBranchAddress("PHI",&vPhi);
	ch1->SetBranchAddress("CHARGE",&vCharge);
	
	//Here we define the variables for each particle
	Int_t pid,aid;
	int nTrigger = 0;
	Double_t pPt,pEta,pPhi,pCharge;//For triger
	Double_t aPt,aEta,aPhi,aCharge;//For associate
	Int_t trigL,trigI,trigH,pairLL,pairIL, pairII, pairHL, pairHI, pairHH;//Trigger particle counter variable
	
	int nEvents = ch1->GetEntries();
	
	cout<<"The number of events for this analysis is: "<<nEvents<<endl; 
	
	//Setting up histograms
	//General simulation variables distribution
	TH1D *hPt = new TH1D("hPt","Transverse momentum distribution of all the charm partile produced;p_{T} (Gev/c);Counts",100,0,50);
	TH1D *hEta = new TH1D("hEta","#eta distribution of all the charm particle  produced;#eta;Counts",40,-4,4);
	TH1D *hPhi = new TH1D("hPhi","Azimuthial angle distribution of all charm particle produced;#phi (rad);Counts",60,-PI,PI);
	TH1D *hID = new TH1D("hID","ID distribution of all the charm particle produced;ID;Counts",12001,-6000.5,6000.5);
	TH1D *hCharge = new TH1D("hCharge","Charge distribution of charm particles;q (e);Counts",100,-4.5,4.5);
	//Multidimensional correlation plots
	TH3D *hPtrPaDPhi = new TH3D("hDPtrPaDPhi",Form("#Delta#phi p_{T} Triger p_{T} Associate for %s;p_{T} triger (Gev/c); p_{T} associate (GeV/c); #Delta#phi (rad)",title),100,0,50,100,0,50,100,-PI/2,3*PI/2);
	TH2D *hDeltaPhiPt = new TH2D("hDeltaPhiPt",Form("%s #Delta#phi and trigger p_{T};p_{T} (GeV/c); #Delta#phi (rad)",title),100,0,50,100,-PI/2,3*PI/2);
	TH2D *hDEtaPt = new TH2D("hDEtaPt",Form("%s #Delta#eta and trigger p_{T};p_{T} (GeV/c);#Delta#eta",title),100,0,50,80,-8,8);
	TH2D *hDPhiDEta = new TH2D("hDPhiDEta",Form("%s #Delta#phi and #Delta#eta;#Delta#phi (rad);#Delta#eta",title),100,-PI/2,3*PI/2,80,-8,8);
	TH2D *htrIDasID = new TH2D("htrIDasID",Form("%s Triger and Associate ID; trigger ID;associate ID",title),12001,-6000.5,6000.5,12001,-6000.5,6000.5);
	//DeltaPhi correlation plots
	TH1D* hDeltaPhi = new TH1D("hDeltaPhi",Form("%s correlations;#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiLL = new TH1D("hDPhiLL",Form("%s #Delta#phi correlation for low-low p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiIL = new TH1D("hDPhiIL",Form("%s #Delta#phi correlation for intermediate-low p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiII = new TH1D("hDPhiII",Form("%s #Delta#phi correlation for intermediate-intermediate p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiHL = new TH1D("hDPhiHL",Form("%s #Delta#phi correlation for high-low p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiHI = new TH1D("hDPhiHI",Form("%s #Delta#phi correlation for high-intermediate p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	TH1D* hDPhiHH = new TH1D("hDPhiHH",Form("%s #Delta#phi correlation for high-high p_{T};#Delta#phi (rad);Counts",title),100,-PI/2,3*PI/2);
	//DeltaEta correlation plots
	TH1D* hDEta = new TH1D("hDEta",Form("%s #Delta#eta correlation;#Delta#eta;Counts",title),80,-8,8);
	TH1D* hDEtaLL = new TH1D("hDEtaLL",Form("%s #Delta#eta correlation for low-low p_{T};#Delta#eta;Counts",title),80,-8,8);
	TH1D* hDEtaIL = new TH1D("hDEtaIL",Form("%s #Delta#eta correlation for intermediate-low p_{T};#Delta#eta;Counts",title),80,-8,8);
	TH1D* hDEtaII = new TH1D("hDEtaII",Form("%s #Delta#eta correlation for intermediate-intermediate p_{T};#Delta#eta;Counts",title),80,-8,8);
	TH1D* hDEtaHL = new TH1D("hDEtaHL",Form("%s #Delta#eta correlation for high-low p_{T};#Delta#eta;Counts",title),80,-8,8);
	TH1D* hDEtaHI = new TH1D("hDEtaHI",Form("%s #Delta#eta correlation for high-intermediate p_{T};#Delta#eta;Counts",title),80,-8,8);
	TH1D* hDEtaHH = new TH1D("hDEtaHH",Form("%s #Delta#eta correlation for high-high p_{T};#Delta#eta;Counts",title),80,-8,8);
	//Trigger counter plots
	TH1D* hTrigL = new TH1D("hTrigL",Form("%s Triggers for low p_{T};#frac{Triggers}{Event};Counts",title),11,-0.5,10.5);
	TH1D* hTrigI = new TH1D("hTrigI",Form("%s Triggers for intermediate p_{T};#frac{Triggers}{Event};Counts",title),11,-0.5,10.5);
	TH1D* hTrigH = new TH1D("hTrigH",Form("%s Triggers for high p_{T};#frac{Triggers}{Event};Counts",title),11,-0.5,10.5);
	//Pair counter
	TH1D* hPairLL = new TH1D("hPairLL",Form("%s Number of pairs for Low-Low correlations;#frac{Pairs}{Event};Events",title),11,-0.5,10.5);
	TH1D* hPairIL = new TH1D("hPairIL",Form("%s Number of pairs for Intermediate-Low correlations;#frac{Pairs}{Event};Events",title),11,-0.5,10.5);
	TH1D* hPairII = new TH1D("hPairII",Form("%s Number of pairs for Intermediate-Intermediate correlations;#frac{Pairs}{Event};Event",title),11,-0.5,10.5);
	TH1D* hPairHL = new TH1D("hPairHL",Form("%s Number of pairs for High-Low correlations;#frac{Pairs}{Event};Event;",title),11,-0.5,10.5);
	TH1D* hPairHI = new TH1D("hPairHI",Form("%s Number of pairs for High-Intermediate correlations;#frac{Pairs}{Event};Event",title),11,-0.5,10.5);
	TH1D* hPairHH = new TH1D("hPairHH",Form("%s Number of pairs for High-High correlations;#frac{Pairs}{Events};Event",title),11,-0.5,10.5);
	
	//Begining of event loop
	for(int iEvent = 0; iEvent < nEvents; iEvent++){
		ch1->GetEntry(iEvent);
		int nparticles = vID->size();
		//Initializing trigger counters
		trigL = 0;
		trigI = 0;
		trigH = 0;
		pairLL = 0;
		pairIL = 0;
		pairII = 0;
		pairHL = 0;
		pairHI = 0;
		pairHH = 0;

		//Starting particle loop
		for(int iparticle = 0; iparticle < nparticles; iparticle++){
			pid = (*vID)[iparticle];
			pPhi = (*vPhi)[iparticle];
			pPt = (*vPt)[iparticle];
			pEta = (*vEta)[iparticle];
			pCharge = (*vCharge)[iparticle];
			//Filling general histograms
			hID->Fill((Double_t) pid);
			hPhi->Fill(pPhi);
			hPt->Fill(pPt);
			hEta->Fill(pEta);
			hCharge->Fill(pCharge);
			//Trigger loop
			if(pid == id_trigger){
				if(pPt >= 1.) nTrigger++;
				if(pPt >= 1. && pPt < 3.) trigL++;
				if(pPt >= 3. && pPt < 8.) trigI++;
				if(pPt >= 8.) trigH++;
				for(int jparticle = 0; jparticle < nparticles; jparticle++){
					if (iparticle == jparticle) continue;//So we do not correlate with itself.
						aid = (*vID)[jparticle];
						aPhi = (*vPhi)[jparticle];
						aEta = (*vEta)[jparticle];
						aPt = (*vPt)[jparticle];
						if(aid == id_associate){
							//Filling general correlationplots
							hPtrPaDPhi->Fill(pPt,aPt,DeltaPhi(pPhi,aPhi));
							hDeltaPhi->Fill(DeltaPhi(pPhi,aPhi));
							hDeltaPhiPt->Fill(pPt,DeltaPhi(pPhi,aPhi));
							hDPhiDEta->Fill(DeltaPhi(pPhi,aPhi),pEta-aEta);
							hDEtaPt->Fill(pPt,pEta-aEta);
							hDEta->Fill(pEta-aEta);
							htrIDasID->Fill(pid,aid);
							//Filling triger momentum range correlations
							if(pPt >= 1. && pPt < 3. && aPt > 1. && aPt < 3.){
								hDPhiLL->Fill(DeltaPhi(pPhi,aPhi));
								hDEtaLL->Fill(pEta-aEta);
								pairLL++;
								
								
							}//End of low transverse momentum
							else if(pPt >= 3. && pPt < 8. && aPt < 3.){
								hDPhiIL->Fill(DeltaPhi(pPhi,aPhi));
								hDEtaIL->Fill(pEta-aEta);
								pairIL++;
								//htrIDasID->Fill(pid,aid);
																
							}//End of intermidiate-low transverse momentum
							
							else if(pPt >= 3. && pPt < 8. && aPt >= 3. && aPt < 8.){
								hDPhiII->Fill(DeltaPhi(pPhi,aPhi));
								hDEtaII->Fill(pEta-aEta);
								pairII++;
								//htrIDasID->Fill(pid,aid);
								
							}
							
							else if(pPt >= 8. && aPt < 3.){
								hDPhiHL->Fill(DeltaPhi(pPhi,aPhi));
								hDEtaHL->Fill(pEta-aEta);
								pairHL++;
								//htrIDasID->Fill(pid,aid);
								
							}//End of high-low transverse momentum
							
							else if(pPt >= 8. && aPt >=3. && aPt < 8.){
								hDPhiHI->Fill(DeltaPhi(pPhi,aPhi));
								hDEtaHI->Fill(pEta-aEta);
								pairHI++;
								//htrIDasID->Fill(pid,aid);
								
							}//End of high-intermediate transverse momentum
							
							else if(pPt >= 8. && aPt >= 8.){
								hDPhiHH->Fill(DeltaPhi(pPhi,aPhi));
								hDEtaHH->Fill(pEta-aEta);
								pairHH++;
								//htrIDasID->Fill(pid,aid);
							
							}//End of high-high transverse momentum
					}//End of DD correlation
				}//End of associate particle loop D+D+
			}//D+ triger
			
		}//End of  particle loop
		hTrigL->Fill((Double_t)trigL);
		hTrigI->Fill((Double_t)trigI);
		hTrigH->Fill((Double_t)trigH);
		hPairLL->Fill((Double_t)pairLL);
		hPairIL->Fill((Double_t)pairIL);
		hPairII->Fill((Double_t)pairII);
		hPairHL->Fill((Double_t)pairHL);
		hPairHI->Fill((Double_t)pairHI);
		hPairHH->Fill((Double_t)pairHH);
			
	}//End of event loop
	if(nTrigger == 0){
		cout<<"Have not found any trigger particle with id: "<<id_trigger<<endl;
		output->Close();
		return;
	}
	output->Write();
	output->Close();
	cout<<"The total number of triggers is: "<<nTrigger<<endl;
	
	cout<<"File: "<<filename<<" has been created!"<<endl;
	
}//End of correlation plot function

void plot_all_combinations(TString filename1, TString filename2, TString filename3){
//This functions takes the three filenames as input reads them takes the correlation and trigger
//histograms normalizes to the number of triggers and plot the correlations of each kind for a
//specific momentum range
//always put same sign correlations in filename2 and filename3
	//Allocations of arrays
	TH1D *DPhi_Hist1[6];
	TH1D *DPhi_Hist2[6];
	TH1D *DPhi_Hist3[6];
	TH1D *hTriggers1[3];
	TH1D *hTriggers2[3];
	TH1D *hTriggers3[3];
	TLegend *leg[6];
	TCanvas *c1[6];
	
	//Dictionaries used
	const char* DPhi_Names[6] = {"hDPhiLL","hDPhiIL","hDPhiII","hDPhiHL","hDPhiHI","hDPhiHH"};
	const char* Triggers_Names[3] = {"hTrigL","hTrigI","hTrigH"};
	const char* titles[6] = {"Low-Low p_{T} correlations","Intermediate-Low p_{T} correlations","Intermediate-Intermediate  p_{T} correlations","High-Low p_{T} correlations","High-Intermediate p_{T} correlations","High-High p_{T} correlations"};
	
	//Reading the files
	TFile *input1 = TFile::Open(filename1);
	if(!input1->IsOpen() || !input1){
		cout<<"File: "<<filename1<<" is not found!"<<endl;
		cout<<"Terminating!"<<endl;
		return 0;
	}
	
	TFile *input2 = TFile::Open(filename2);
	if(!input2->IsOpen() || !input2){
		cout<<"File: "<<filename2<<" is not found!"<<endl;
		cout<<"Terminating!"<<endl;
		return 0;
	}
	
	TFile *input3 = TFile::Open(filename3);
	if(!input3->IsOpen() || !input3){
		cout<<"File: "<<filename3<<" is not found!"<<endl;
		cout<<"Terminating!"<<endl;
		return 0;
	}
	
	//Importing histograms
	for(int i = 0; i<6; i++){
		DPhi_Hist1[i] = (TH1D*)input1->Get(DPhi_Names[i]);
		if(!DPhi_Hist1[i]){
			cout<<DPhi_Names[i]<<" not found in "<<filename1<<endl;
			cout<<"Terminating"<<endl;
		}
		DPhi_Hist2[i] = (TH1D*)input2->Get(DPhi_Names[i]);
		if(!DPhi_Hist2[i]){
			cout<<DPhi_Names[i]<<" not found in "<<filename2<<endl;
			cout<<"Terminating"<<endl;
		}
		DPhi_Hist3[i] = (TH1D*)input3->Get(DPhi_Names[i]);
		if(!DPhi_Hist3[i]){
			cout<<DPhi_Names[i]<<" not found in "<<filename3<<endl;
			cout<<"Terminating"<<endl;
		}
		if(i<3){
			hTriggers1[i] = (TH1D*)input1->Get(Triggers_Names[i]);
			if(!hTriggers1[i]){
				cout<<Triggers_Names[i]<<" not found in "<<filename1<<endl;
				cout<<"Terminating!"<<endl;
				return 0;
			}
			hTriggers2[i] = (TH1D*)input2->Get(Triggers_Names[i]);
			if(!hTriggers2[i]){
				cout<<Triggers_Names[i]<<" not found in "<<filename2<<endl;
				cout<<"Terminating!"<<endl;
				return 0;
			}
			hTriggers3[i] = (TH1D*)input3->Get(Triggers_Names[i]);
			if(!hTriggers3[i]){
				cout<<Triggers_Names[i]<<" not found in "<<filename3<<endl;
				cout<<"Terminating!"<<endl;
				return 0;
			}			
		
		}//End of Importing triggers
	}
	
	//Normalizing
	for(int i = 0; i<6; i++){
               //Normlizing to number of particles
		 if(i == 0){
			DPhi_Hist1[i]->Scale(1. /hTriggers1[0]->Integral(2,10));
			DPhi_Hist2[i]->Scale(1. /hTriggers2[0]->Integral(2,10));
			DPhi_Hist3[i]->Scale(1. /hTriggers3[0]->Integral(2,10));
			
		}
		if(i == 1 || i == 2){
			DPhi_Hist1[i]->Scale(1. /hTriggers1[1]->Integral(2,10));
			DPhi_Hist2[i]->Scale(1. /hTriggers2[1]->Integral(2,10));
			DPhi_Hist3[i]->Scale(1. /hTriggers3[1]->Integral(2,10));
		}
		if(i == 3 || i == 4 || i == 5){
			DPhi_Hist1[i]->Scale(1. /hTriggers1[2]->Integral(2,10));
			DPhi_Hist2[i]->Scale(1. /hTriggers2[2]->Integral(2,10));
			DPhi_Hist3[i]->Scale(1. /hTriggers3[2]->Integral(2,10));
		}
	
		//Cosmetics
		DPhi_Hist1[i]->GetYaxis()->SetTitle("#frac{1}{N_{tr}} #frac{dN_{assoc}}{d#Delta#phi}");
		DPhi_Hist2[i]->GetYaxis()->SetTitle("#frac{1}{N_{D^{+}}} #frac{dN_{assoc}}{d#Delta#phi}");
		DPhi_Hist3[i]->GetYaxis()->SetTitle("#frac{1}{N_{D^{-}}} #frac{dN_{assoc}}{d#Delta#phi}");
		
		DPhi_Hist1[i]->SetMarkerStyle(7);
		DPhi_Hist2[i]->SetMarkerStyle(7);
		DPhi_Hist3[i]->SetMarkerStyle(7);
		DPhi_Hist1[i]->SetMarkerColor(2);
		DPhi_Hist2[i]->SetMarkerColor(3);
		DPhi_Hist3[i]->SetMarkerColor(4);
		DPhi_Hist1[i]->SetTitle(titles[i]);
		DPhi_Hist2[i]->SetTitle(titles[i]);
		DPhi_Hist1[i]->SetStats(0);
		DPhi_Hist2[i]->SetStats(0);
		DPhi_Hist3[i]->SetStats(0);
	}
	
	//Legends
	for(int i = 0; i<6; i++){
		leg[i]=new TLegend(0.89,0.7,0.99,0.9);
		leg[i]->AddEntry(DPhi_Hist1[i],"D^{+}D^{-}","P");
		leg[i]->AddEntry(DPhi_Hist2[i],"D^{+}D^{+}","P");
		leg[i]->AddEntry(DPhi_Hist3[i],"D^{-}D^{-}","P");
		leg[i]->SetBorderSize(0);
	}
	TCanvas *c2 = new TCanvas();
	c2->Divide(3,2);
	
	//Plotting
	
	for(int i = 0; i<6; i++){
		gStyle->SetPadLeftMargin(0.18); gStyle->SetPadRightMargin(0.01);
		//c1[i] = new TCanvas(DPhi_Hist1[i]->GetName(),DPhi_Hist1[i]->GetName(),1920,1080);
		//c1[i]->Divide(2,1);
		//c1[i]->cd(1);
		//DPhi_Hist2[i]->GetYaxis()->SetRangeUser(0.00001,0.0025);
		//DPhi_Hist2[i]->Draw("PE1");
		//DPhi_Hist3[i]->Draw("PE1 same");
		//DPhi_Hist1[i]->Draw("PE1 same");
		//c1[i]->cd(2);
		//DPhi_Hist2[i]->Draw("PE1");
		//DPhi_Hist3[i]->Draw("PE1 same");
		//leg[i]->Draw();
		c2->cd(i+1);
		DPhi_Hist1[i]->Draw("PE1");
		DPhi_Hist3[i]->Draw("PE1 same");
		DPhi_Hist2[i]->Draw("PE1 same");
		leg[i]->Draw();
	}
	
	for(int i = 0;i<3;i++){
		cout<<"The integral of trigers in is: "<<hTriggers2[i]->Integral(2,10)<<endl;
		cout<<"The integral of at least two particles is: "<<hTriggers2[i]->Integral(3,10)<<endl;
		cout<<"Printing bin content:"<<endl;
		for(int j = 0; j < 10; j++){
			cout<<" Bin content for histogram "<<hTriggers1[i]->GetTitle()<<" Bin: "<<j<<" Bin Center: "<<hTriggers1[i]->GetBinCenter(j)<<" Content: "<<hTriggers1[i]->GetBinContent(j)<<" events. "<<endl;
		}
	}
	
	for (int i = 0; i<6; i++){
		cout<<"The number of entries for SS is: "<< DPhi_Hist2[i]->GetEntries()<<" for OS: "<<DPhi_Hist1[i]->GetEntries()<<endl;
	
	}

}//End of plot_all_combinations functions

void plot_particle_number(TString filename){
	TH1D *hTriggers1[3];
	const char* Triggers_Names[3] = {"hTrigL","hTrigI","hTrigH"};
	//Opening file
	TFile *input1 = TFile::Open(filename);
	if(!input1->IsOpen() || !input1){
		cout<<"File: "<<filename<<" is not found!"<<endl;
		cout<<"Terminating!"<<endl;
		return 0;
	}
	//Importing histograms
	for(int i = 0; i < 3; i++){
		hTriggers1[i] = (TH1D*)input1->Get(Triggers_Names[i]);
				if(!hTriggers1[i]){
					cout<<Triggers_Names[i]<<" not found in "<<filename<<endl;
					cout<<"Terminating!"<<endl;
					return 0;
				}
		
	}
	
	hTriggers1[0]->SetTitle("D^{+} production in 1-3 GeV/c");
	hTriggers1[0]->GetXaxis()->SetTitle("#frac{Paricles}{Event}");
	hTriggers1[1]->SetTitle("D^{+} production in 3-8 GeV/c");
	hTriggers1[1]->GetXaxis()->SetTitle("#frac{Paricles}{Event}");
	hTriggers1[2]->SetTitle("D^{+} production in p_{T}>8 GeV/c");
	hTriggers1[2]->GetXaxis()->SetTitle("#frac{Paricles}{Event}");
	
	TCanvas *c2 = new TCanvas();
	c2->Divide(3,1);
	for(int i = 0; i < 3; i++){
		c2->cd(i+1)->SetLogy();
		hTriggers1[i]->Draw();
	}
	
}

void triger_analysis(){
	//correlation_plot(411,-411,"DDbar_correlations.root","D^{+}D^{-}");
	//correlation_plot(411,411,"DD_correlations.root","D^{+}D^{+}");
	//correlation_plot(-411,-411,"DbarDbar_correlations.root","D^{-}D^{-}");
	
	plot_all_combinations("DDbar_correlation_MONASH.root","DD_correlation_MONASH.root","DbarDbar_correlation_MONASH.root");
	//plot_all_combinations("DDbar_correlations.root","DD_correlations.root","DbarDbar_correlations.root");
//	plot_particle_number("DD_correlation_MONASH.root");
}

