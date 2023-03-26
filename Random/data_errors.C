Int_t one_event(TH1D *hist){
//This function returns the bin with one events
	Int_t c = hist->GetBinCenter(0);
	int i = 0;
	while(c < 2){
		c = hist->GetBinCenter(i);
		if(c == 1){
			break ;
			}
		i++;
		}
	cout<<"The scalling is for bin "<<i<<" and center "<<c<<endl;
	return i;
	}
void scale(TH1D *hist, Int_t bin){
//This function normalizes histograms to the number of events that produced at least one particle.
	hist->Scale(1/hist->Integral(bin,hist->GetNbinsX()));
	}



void data_errors(){
using namespace std;
using namespace TMath;

//Importing the file
TFile *input = TFile::Open("merged_hf2.root");
if(!input->IsOpen()|| !input){
	cout<<"Input not found"<<endl;
	return;
	}

//Importing histograms
TH1D *deltaphi = (TH1D*)input->Get("hLambdaDPhi");
if(!deltaphi){
	cout<<"Histogram not found"<<endl;
	return;
	}

TH1D *multiplicity = (TH1D*)input->Get("hNparticles");
if(!multiplicity){
	cout<<"Histogram not found"<<endl;	
	}
//Normalizing to number of events that generated at least one particle.
scale(deltaphi,one_event(multiplicity));

//Defining graph for error to data comparisson
TGraph *gr = new TGraph();
Double_t x,y,err;

//Filling the Graph
for(int i = 1; i < deltaphi->GetNbinsX()+1;i++){
	x = deltaphi->GetBinCenter(i);
	y = deltaphi->GetBinContent(i);
	err = deltaphi->GetBinError(i);
	gr->SetPoint(gr->GetN(),x,err/y);	
	}//end of filling graph loop
	
//Histogram and Graph cosmetics
deltaphi->GetYaxis()->SetTitle("#frac{1}{N_{events}} #frac{dN}{d#varphi}");
deltaphi->SetTitle("#Delta#phi between #Lambda_{c} and every other charm particle");
gr->GetXaxis()->SetTitle("#Delta#phi (Rad)");
gr->GetYaxis()->SetTitle("Error/Value");
gr->SetMarkerStyle(7);


//Drawing everything
TCanvas *c1 = new TCanvas();
c1->Divide(1,2,0,0);

c1->cd(1);
deltaphi->Draw("E1");

c1->cd(2);
gr->Draw("APC");



}
