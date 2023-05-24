//This script will provide trees with the pT phi eta PDG codes so they can be used to create more correlatiion plots. 
//In principle every vector is an event and every element of the vector is a particle produced at that event.  It also 
//creates a D+D+ correlation histogram so we can compare it with histogram created by the trees to check if everything works
//alright. Finally it produces a multiplicity histogram and charm particles per event histogram.

//C++ libraries we are using
#include <iostream>
#include <cmath>
#include <cstring>
#include <chrono>
#include <vector>
//Root and pythia libraries
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TTree.h"

#define PI 3.14159265
using namespace std;
using namespace Pythia8;

//Here we define some functions we are going to use
bool IsCharm( Int_t particlepdg){
	//This function will test if a particle has any quark content
	Int_t pdg = abs(particlepdg);
	pdg /= 10;//Last digid does not have to do with quark content
	if(pdg % 10 == 4) return true; //3rd quark
	pdg /= 10;
	if(pdg % 10 == 4) return true; //2nd quark
	pdg /= 10;
	if(pdg % 10 == 4) return true; //1st quark
	return false;
	}

Double_t DeltaPhi(Double_t phi1, Double_t phi2){
	//returns delta phi in range -pi/2 3pi/2
	return fmod(phi1-phi2+2.5*PI,2*PI)-0.5*PI;
	}


int main(int argc, char** argv){

	if(argc != 2){//main sould have two arguments one will be the ./CharmDeltaPhi and one other will be the filename or path	
		cout<<"Error in the number of arguments provided"<<endl;
		cout<<"Provide only the filepath/name."<<endl;
		cout<<"Terminating program"<<endl;
		return 0;
	}

	//Here we start keeping track of time
	auto start = chrono::high_resolution_clock::now();
	
	//Create output file
	TFile* output = new TFile(argv[1],"CREATE");
	if(!output->IsOpen()){
		cout<<"Error: File "<<argv[1]<<"already exists terminating program!"<<endl;
		return 1;
	}
	
	//Here I define the tree that the data will be stored.
	TTree *tree = new TTree("tree","ccbar correlations");

	
	//here we define the variables we are going to need.
	Double_t pT, eta, phi, pTTriger, pTAssociate, charge, DeltaPhiDD,status;
	Int_t  id,idCharm,size,nEvents,charmness;
	
	//Here I define the vectors
	vector<Int_t> vID;
	vector<Double_t> vPt;
	vector<Double_t> vEta;
	vector<Double_t> vPhi;
	vector<Double_t> vCharge;
	vector<Double_t> vStatus;

	
	//Setting up tree Branches and histograms
	tree->Branch("ID",&vID);
	tree->Branch("PT",&vPt);
	tree->Branch("ETA",&vEta);
	tree->Branch("PHI",&vPhi);
	tree->Branch("CHARGE",&vCharge);
	tree->Branch("STATUS",&vStatus);
	tree->Branch("MULTIPLICITY",&size,"x/I");
	TH1D* hSize = new TH1D("hSize","Multiplicity",301,-0.5,300.5);//Corrected it was 300!
	TH1D* hidCharm = new TH1D("hidCharm","PDG Codes for Charm hadrons",12000,-6000,6000);
	TH1D* hPtTriger = new TH1D("hPtTriger","p_{T} for triger D^{+} ",50,0,10);
	TH1D* hPtAssociate = new TH1D("hPtAssociate", "p_{T} for associate D^{+}",50,0,10);
	TH1D* hDeltaPhiDD = new TH1D("hDeltaPhiDD","D^{0}#barD^{0} correlations",100,-PI/2,3*PI/2);
	TH1D* hCharmPart = new TH1D("hCharmPart", "Charm Particles Per Event",200,-0.5,200.5);

	
	//Kinematics constraints
	const Double_t pTmin = 0.15;//minimum pT
	const Double_t etamax = 4.;// maximum eta
	
	//The simulation is an object so we define it like
	Pythia pythia;
	
	//Simulation settings from a file
	pythia.readFile("pythiasettings.cmnd");
	nEvents = pythia.mode("Main:numberOfEvents");
	
	//Here we create a radnom seed using the runtime so eachtime the simulation runs the outcome will be trully radnom.
	Int_t proccessid = getpid();
	string seedstr = "Random:seed = "+std::to_string((time(0)+proccessid)%900000000);
	pythia.readString("Random:setSeed = on");
	pythia.readString(seedstr);
	
	//Initializing simulation.
	pythia.init();//important must be at every simulation
	
	cout<<"Generating "<<nEvents<<" events!"<<endl;
	
	//Here we start the event loop
	for(int iEvent = 0; iEvent<nEvents; iEvent++){
		if(!pythia.next()) continue; //generate the next event if there is an error in the current. This command must exist at every pythia script.
	
		int nPart = pythia.event.size();//particles produced in this event
		size = 0;//initialiazing for multiplicity plot.
		charmness = 0;//intialiazing for charm production plot.
		//Initializing vectors
		vID.clear();
		vPt.clear();
		vEta.clear();
		vPhi.clear();
		vCharge.clear();
		vStatus.clear();
		//Particle loop
		for(int iPart = 0; iPart<nPart; iPart++){
			const Particle &particle = pythia.event[iPart];
			if(!particle.isFinal()) continue; //skip if the particle is not at its final state.
			
			id = particle.id();
			pT = particle.pT();
			eta = particle.eta();
			phi = particle.phi();
			charge = particle.charge();
			status = static_cast<Double_t> (particle.status());
			
			
			//kinematics check
			if(pT < pTmin || abs(eta) > etamax ) continue;
			
			if(charge != 0) size++;
			
			if(!IsCharm(id)) continue;
				idCharm = id;
				hidCharm->Fill((Double_t) id);
				charmness++;
				//Filling vectors
				vID.push_back(id);
				vPt.push_back(pT);
				vEta.push_back(eta);
				vPhi.push_back(phi);
				vCharge.push_back(charge);
				vStatus.push_back(status);
				
			//Creating correlation plots
				if(id==421 && pT >2.){//Dzero meson triger
				
					for(int jPart =0; jPart<nPart; jPart++){
						const Particle particleD = pythia.event[jPart];
						if(jPart == iPart) continue;//So we do not correlate a particle with itself.
						Int_t associate_id = particleD.id();
						Double_t associate_pT = particleD.pT();
						Double_t associate_eta = particleD.eta();
						if(associate_pT < pTmin  || abs(associate_eta) > etamax) continue;
						if( associate_id == -421){//DD correlation
							pTTriger = pT;
							DeltaPhiDD = DeltaPhi(phi,particleD.phi());
							hPtTriger->Fill(pT);
							hPtAssociate->Fill(associate_pT);
							hDeltaPhiDD->Fill(DeltaPhiDD);					
						}	
						else continue;	
					}// End of D meson triger particle loop.
				}//End of D meson triger
		}//1st particle loop
		hSize->Fill((Double_t) size);
		hCharmPart->Fill((Double_t) charmness);
		//In order not fill trees with empty vectors.
		if(vID.empty() || vPt.empty() || vEta.empty() || vPhi.empty() || vCharge.empty() || vStatus.empty() ) continue; 
		tree->Fill();
	}//End of event loop
	
	//write output and close it
	
	output->Write();
	cout<<"File has been created and its name is: "<<output->GetName()<<endl;
	output->Close();
	
	//Stop keeping time
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::minutes>(end - start);
        cout << "This script took " << duration.count() << " minutes to run." << endl;
	
	return 0;
}




