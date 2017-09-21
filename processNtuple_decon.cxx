//compile independently with: g++ -std=c++11 -o processNtuple_decon processNtuple_decon.cxx `root-config --cflags --glibs` -lMinuit
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
using namespace std;

#include "TROOT.h"
#include "TMath.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile2D.h"
#include "TString.h"
#include "TCanvas.h"
#include <TGraphErrors.h>
#include "TMultiGraph.h"
#include "TImage.h"
#include "TSystem.h"
#include <TVirtualFFT.h>

#include "FitFeElecResponse.hxx"

using namespace std;

//global TApplication object declared here for simplicity
TApplication *theApp;

class Analyze {
	public:
	Analyze(std::string inputFileName);
	void SetupNominalResponse();
	void SetupFilters();
	void DrawResponse();
	void doAnalysis();
	void getResponseWaveform();
	void analyzeChannel();
	void simChannel();
	void deconChannel();

	TFile* inputFile;
	TFile *statusfile;
	TFile *gOut;

	//ROI tr_rawdata variables
	TTree *tr_rawdata;
	std::vector<short> *adc_v = 0;
	unsigned short run,subrun, event, chan;

	//Constants
	const int numAsic = 512; //35t
	const int numChan = 8256;// 35t
	const int maxNumBin = 5000; //35t
	const Int_t maxTicks = 9594;

	TGraphErrors *gTest;
	TCanvas* c0;
	TGraph *gCh;
	int eventNum;

	TH2D *pCorrWaveVsCh;
	TH1D *respWf;
        TH1F *hWf;
	TGraph *gResp;
	TGraph *gFFT_mag;
	TGraph *gFFT_phase;
	TGraph *gFFT_nominal_mag;
	TGraph *gFFT_nominal_phase;
	TGraph *gCh_corr;
	TGraph *gNominal;

	FeElecResponse *elecRespSig;

	TF1 *fFilterU;
  	TF1 *fFilterV;
  	TF1 *fFilterW;

	//histograms
	TH2F *hSampVsChan;
	TProfile *pSampVsChan;
	TH2F *hMeanVsChan;
	TProfile *pMeanVsChan;
	TH2F *hRmsVsChan;
	TProfile *pRmsVsChan;

	TH2F *hCorrVsChan;
};

Analyze::Analyze(std::string inputFileName){

	//get input file
	if( inputFileName.empty() ){
		std::cout << "Error invalid file name" << std::endl;
		return;
	}

	inputFile = new TFile(inputFileName.c_str());
	if (inputFile->IsZombie()) {
		std::cout << "Error opening input file" << std::endl;
		return;
	}

	if( !inputFile ){
		std::cout << "Error opening input file" << std::endl;
		return;
	}

	tr_rawdata = (TTree*) inputFile->Get("GetWF/fTree");
  	if( !tr_rawdata ){
		std::cout << "Error opening input file tree" << std::endl;
		return;
  	}
	tr_rawdata->SetBranchAddress("_run", &run);
  	tr_rawdata->SetBranchAddress("_chan", &chan);
  	tr_rawdata->SetBranchAddress("_event", &event);
  	tr_rawdata->SetBranchAddress("adc_v",&adc_v);

	//open response file
        TFile *respFile = new TFile("responseWaveforms.root");
	if (respFile->IsZombie()) {
		std::cout << "Error opening response file" << std::endl;
		return;
	}

	if( !respFile ){
		std::cout << "Error opening response file" << std::endl;
		return;
	}
	pCorrWaveVsCh = (TH2D*) respFile->Get("pCorrWaveVsCh");

	//make output file
  	std::string outputFileName = "output_processNtuple_decon_" + inputFileName;
  	gOut = new TFile(outputFileName.c_str() , "RECREATE");

  	//initialize canvas
  	c0 = new TCanvas("c0", "c0",1400,800);

  	//initialize graphs + profiles
  	gTest = new TGraphErrors();
  	gCh = new TGraph();
	hWf = new TH1F("hWf","",maxTicks,-0.5,maxTicks-0.5);
	gResp = new TGraph();
	gFFT_mag = new TGraph();
	gFFT_phase = new TGraph();
	gFFT_nominal_mag = new TGraph();
	gFFT_nominal_phase = new TGraph();
	gCh_corr = new TGraph();
	gNominal = new TGraph();

	//electronics response object
	elecRespSig = new FeElecResponse();

  	//histograms + grraphs
  	hSampVsChan = new TH2F("hSampVsChan","",numChan,0-0.5,numChan-0.5,4096,-0.5,4096-0.5);
 	pSampVsChan = new TProfile("pSampVsChan","",numChan,0-0.5,numChan-0.5);
  	hMeanVsChan = new TH2F("hMeanVsChan","",numChan,0-0.5,numChan-0.5,4096,-0.5,4096-0.5);
	pMeanVsChan = new TProfile("pMeanVsChan","",numChan,0-0.5,numChan-0.5);
  	hRmsVsChan = new TH2F("hRmsVsChan","",numChan,0-0.5,numChan-0.5,300,0,300.);
  	pRmsVsChan = new TProfile("pRmsVsChan","",numChan,0-0.5,numChan-0.5);

	hCorrVsChan = new TH2F("hCorrVsChan","",numChan,0-0.5,numChan-0.5,maxTicks,0-0.5,maxTicks-0.5);
}

void Analyze::doAnalysis(){
  	//loop over tr_rawdata entries
  	Long64_t nEntries(tr_rawdata->GetEntries());

	SetupFilters();
	SetupNominalResponse();

	tr_rawdata->GetEntry(0);
	int currEvent = event;
	int numEvent = 0;
	//loop over pulse waveform
	for(Long64_t entry(0); entry<nEntries; ++entry) { 
		tr_rawdata->GetEntry(entry);

		//try to get response waveform
		getResponseWaveform();
		if( respWf->GetEntries() == 0 )
			continue;



		std::cout << entry << "\t" << nEntries << std::endl;

		//analyze current entry
    		analyzeChannel();
		//simChannel();

		//deconvolution
		deconChannel();
  	}//entries

    	gOut->Cd("");
  	hSampVsChan->Write();
	pSampVsChan->Write();
  	hMeanVsChan->Write();
	pMeanVsChan->Write();
  	hRmsVsChan->Write();
  	pRmsVsChan->Write();
	hCorrVsChan->Write();
  	gOut->Close();
}

void Analyze::DrawResponse()
{
	c0->Clear();
	c0->Divide(1,2);
	c0->cd(1);
	gResp->GetXaxis()->SetRangeUser(0,35);
	gResp->GetYaxis()->SetRangeUser(-0.001,0.010);
	gResp->Draw("ALP");
	c0->cd(2);
	gNominal->GetXaxis()->SetRangeUser(0,35);
	gNominal->GetYaxis()->SetRangeUser(-0.001,0.010);
	gNominal->SetLineColor(kRed);
	gNominal->Draw("ALP");
	c0->Update();
	char cr;
	std::cin >> cr;
}


void Analyze::SetupFilters()
{
  TF1 *filter_u = new TF1("filter_u","(x>0.0)*gaus*exp(-0.5*pow(x/[3],[4]))");
  double par[5] = {1.73/0.959301,1.69,1.55,0.19,3.75};
  filter_u->SetParameters(par);
  fFilterU = (TF1*)filter_u->Clone("fFilterU");
  
  TF1 *filter_v = new TF1("filter_v","(x>0.0)*gaus*exp(-0.5*pow(x/[3],[4]))");
  double par1[5] = {1.74/0.941034,1.46,1.33,0.23,4.89};
  filter_v->SetParameters(par1);
  fFilterV = (TF1*)filter_v->Clone("fFilterU");
  
  TF1 *filter_w = new TF1("filter_w","(x>0.0)*[0]*exp(-0.5*pow(pow((x-[1])/[2],2),[3]))");
  double par2[4] = {1.03/0.995635,0.08,0.15,2.17};
  filter_w->SetParameters(par2);
  fFilterW = (TF1*)filter_w->Clone("fFilterU");

  delete filter_u;
  delete filter_v;
  delete filter_w;

  return;
}

void Analyze::SetupNominalResponse()
{
	//TH1F *currentHist = new TH1F("currentHist","",maxTicks,-0.5,maxTicks-0.5);
	TH1F *currentHist = new TH1F("currentHist","",maxNumBin,-0.5,maxNumBin-0.5);

	gNominal->Set(0);
	for(int i = 0 ; i < currentHist->GetNbinsX() ; i++ ){
		double time = i*0.050; //50ns steps, like pulser response waveforms
		double base = 0;
		//double pulseStart = 0; //us
		double pulseStart = 3.6;
		double shapeTime = 2.2; //nominal shaping time
		double amp = 1./191.6/0.0988165; //amp factor is not pulse height
		
		double simVal = 0;
		elecRespSig->getSignalValue(time, base, pulseStart, shapeTime, amp, simVal);

		//gNominal->SetPoint(gNominal->GetN(),7.2+i*0.1,simVal);
		gNominal->SetPoint(gNominal->GetN(),i*0.1,simVal);
		currentHist->SetBinContent(i+1, simVal );
	}

	TH1 *hm = 0;
        hm = currentHist->FFT(0,"MAG");
  	TH1 *hp = 0;
        hp = currentHist->FFT(0,"PH");

	//nominal waveform period = 50ns
	int numBins = hm->GetNbinsX();
	for(int bin = 0 ; bin < numBins ; bin++ ){
		//hm->GetBinContent();
    		Double_t freq;
		freq = ((Double_t) bin)/(1.0*numBins)*20.0;
		gFFT_nominal_mag->SetPoint(gFFT_nominal_mag->GetN(),freq,hm->GetBinContent(bin+1));
		gFFT_nominal_phase->SetPoint(gFFT_nominal_phase->GetN(),freq,hp->GetBinContent(bin+1));
	}

	delete hm;
	delete hp;
	delete currentHist;
}

void Analyze::getResponseWaveform(){

	char name[200];
	memset(name,0,sizeof(char)*100 );
        sprintf(name,"hAvgSignal_%.5i",chan);
	respWf = pCorrWaveVsCh->ProjectionY(name,chan+1,chan+1);

	gResp->Set(0);
	gFFT_mag->Set(0);
	gFFT_phase->Set(0);
	if( respWf->GetEntries() == 0 )
		return;

	//TH1F *currentHist = new TH1F("currentHist","",maxTicks,-0.5,maxTicks-0.5);
	TH1F *currentHist = new TH1F("currentHist","",maxNumBin,-0.5,maxNumBin-0.5);
	for( int bin = 0 ; bin < respWf->GetNbinsX() ; bin++ )
		currentHist->SetBinContent(bin+1, respWf->GetBinContent(bin+1)/10000. );
	for( int bin = 0 ; bin < currentHist->GetNbinsX() ; bin++ )
		gResp->SetPoint( gResp->GetN() , bin*0.1, currentHist->GetBinContent(bin+1) );

	TH1 *hm = 0;
        hm = currentHist->FFT(0,"MAG");
  	TH1 *hp = 0;
        hp = currentHist->FFT(0,"PH");

	//response waveform period = 50ns
	int numBins = hm->GetNbinsX();
	for(int bin = 0 ; bin < numBins ; bin++ ){
		//hm->GetBinContent();
    		Double_t freq;
		freq = ((Double_t) bin)/(1.0*numBins)*20.0;
		gFFT_mag->SetPoint(gFFT_mag->GetN(),freq,hm->GetBinContent(bin+1));
		gFFT_phase->SetPoint(gFFT_phase->GetN(),freq,hp->GetBinContent(bin+1));
	}

	delete hm;
	delete hp;
	delete currentHist;
}

void Analyze::analyzeChannel(){
	if( adc_v->size() == 0 )
		return;

	//calculate mean
	double mean = 0;
	int count = 0;
	for( int s = 0 ; s < adc_v->size() ; s++ ){
		double value = adc_v->at(s);
		mean += value;
		count++;
	}
	if( count > 0 )
		mean = mean / (double) count;

	//calculate rms
	double rms = 0;
	count = 0;
	for( int s = 0 ; s < adc_v->size() ; s++ ){
		double value = adc_v->at(s);
		rms += (value-mean)*(value-mean);
		count++;
	}	
	if( count > 1 )
		rms = TMath::Sqrt( rms / (double)( count - 1 ) );

	hMeanVsChan->Fill( chan, mean );
	pMeanVsChan->Fill( chan, mean );
	hRmsVsChan->Fill(chan, rms);
	pRmsVsChan->Fill(chan, rms);

	gCh->Set(0);
	for( int s = 0 ; s < adc_v->size() ; s++ )
		gCh->SetPoint(gCh->GetN() , s , adc_v->at(s) );
}

void Analyze::simChannel(){
	if( gResp->GetN() == 0 )
		return;

	double pulseSample = 1001; //us
	double pulseQ = 10000.; //e-
	double base = 432;

	gCh->Set(0);
	for( int s = 0 ; s < maxTicks ; s++ ){
		double time = s*0.5;
		double sampleVal = base;

		if( s - pulseSample > -10 && s - pulseSample < 140 ){
			sampleVal += pulseQ*gResp->Eval( s - pulseSample  );
			//std::cout << s << "\t" <<  pulseQ*gResp->Eval( s - pulseSample ) << std::endl;
		}
		gCh->SetPoint(gCh->GetN() , s , sampleVal );	
	}

	/*
	c0->Clear();
	c0->Divide(1,2);
	c0->cd(1);
	gResp->Draw("ALP");
	c0->cd(2);
	gCh->GetXaxis()->SetRangeUser(900,1050);
	gCh->Draw("ALP");
	c0->Update();
	char ct;
	std::cin >> ct;
	*/
}

void Analyze::deconChannel(){
	if( gCh->GetN() == 0 )
		return;
	if( gFFT_mag->GetN() == 0 )
		return;
	if( gFFT_phase->GetN() == 0 )
		return;

	Int_t numBins = gCh->GetN();
  	TH1F *currentHist = new TH1F("currentHist","",numBins,-0.5,numBins-0.5);
	for( int s = 0 ; s < numBins ; s++ ){
		double dataX,dataY;
		gCh->GetPoint(s,dataX,dataY);
		currentHist->SetBinContent(s+1, dataY );
	}

	TH1 *hm = 0;
  	hm = currentHist->FFT(0,"MAG");
  	hm->SetName("tempFFT_mag");

  	TH1 *hp = 0;
  	hp = currentHist->FFT(0,"PH");
  	hp->SetName("tempFFT_phase");

	Double_t *value_re = new Double_t[maxTicks];
  	Double_t *value_im = new Double_t[maxTicks];

	TF1 *filter;
  	if(chan < 2400)
    		filter = fFilterU;
  	if( chan >= 2400 && chan < 4800 )
		filter = fFilterV;
  	if( chan >= 4800 )
		filter = fFilterW;

  	for(Int_t i = 0; i < numBins; i++)
  	{
    		Double_t freq;
		freq = ((Double_t) i)/(1.0*numBins)*2.0;
		
		/*
		//this is only required if using filter
    		if(i < numBins/2.0)
    		{
      			freq = ((Double_t) i)/(1.0*numBins)*2.0;
    		}
    		else
    		{
      			freq = ((Double_t) numBins-i)/(1.0*numBins)*2.0;
    		}
		*/
		
		double respMag = gFFT_mag->Eval(freq);
		double respPhase = gFFT_phase->Eval(freq);
		double nominalMag = gFFT_nominal_mag->Eval(freq);
		double nominalPhase = gFFT_nominal_phase->Eval(freq);
		double filterVal = filter->Eval(freq);
		if( respMag <= 0 ){ 
			std::cout << "HERE" << std::endl;
			continue; //bad, should not happen
		}

		//do nothing
    		//Double_t rho = hm->GetBinContent(i+1)*filterVal;
		//Double_t rho = hm->GetBinContent(i+1);
    		//Double_t phi = hp->GetBinContent(i+1);		

		//nominal only
    		//Double_t rho = hm->GetBinContent(i+1)/nominalMag*filterVal;
    		//Double_t phi = hp->GetBinContent(i+1)-nominalPhase;

		//elec response
		//Double_t rho = hm->GetBinContent(i+1)/respMag*filterVal;
    		//Double_t phi = hp->GetBinContent(i+1)-respPhase;

		//deconvolution and reconvolution
    		//Double_t rho = hm->GetBinContent(i+1)/respMag*filterVal*nominalMag;
    		Double_t rho = hm->GetBinContent(i+1)/respMag*nominalMag;
    		Double_t phi = hp->GetBinContent(i+1)-respPhase+nominalPhase;
      
    		if(i == 0) 
      			rho = 0.0;

    		value_re[i] = rho*cos(phi)/((Double_t) numBins);
    		value_im[i] = rho*sin(phi)/((Double_t) numBins);
  	} 

  	Int_t nFreqBins = numBins;
  	TVirtualFFT *invCurrentFFTObject = TVirtualFFT::FFT(1,&nFreqBins,"C2R M K");
  	invCurrentFFTObject->SetPointsComplex(value_re,value_im);
   	invCurrentFFTObject->Transform();

  	TH1F *newHist = new TH1F("","",numBins,-0.5,numBins-0.5);
  	newHist = (TH1F*)TH1::TransformHisto(invCurrentFFTObject,0,"Re");

	//output corrected waveform
	gCh_corr->Set(0);
  	for(Int_t i = 0; i < numBins; i++)
  	{
    		//=fHist->SetBinContent(i+1,newHist->GetBinContent(i+1)/(14.0*4096.0/2000.0));
		gCh_corr->SetPoint(gCh_corr->GetN(),i,newHist->GetBinContent(i+1) );
		hCorrVsChan->SetBinContent(chan+1,i+1,newHist->GetBinContent(i+1));
  	}

	//draw result
	if(0){
	char name[200];
	memset(name,0,sizeof(char)*100 );
        sprintf(name,"Waveform Channel # %.5i",chan);

	c0->Clear();
	c0->Divide(1,3);
	c0->cd(1);
	gResp->SetTitle("Electronics Response Waveform");
	gResp->GetXaxis()->SetRangeUser(-10,60);
	gResp->Draw("ALP");
	gNominal->SetLineColor(kRed);
	gNominal->Draw("L");
	c0->cd(2);
	//gCh->GetXaxis()->SetRangeUser(1475,1685);
	//gCh->GetXaxis()->SetRangeUser(1000,1500);
	gCh->GetXaxis()->SetRangeUser(0,500);
	gCh->SetTitle(name);
	gCh->Draw("ALP");
	c0->cd(3);
	gCh_corr->SetTitle("Corrected Waveform");
	//gCh_corr->GetXaxis()->SetRangeUser(1000,1500);
	gCh_corr->GetXaxis()->SetRangeUser(0,500);
	gCh_corr->Draw("ALP");
	c0->Update();

	char ct;
	std::cin >> ct;
	}

	delete currentHist;
	delete hm;
	delete hp;
	delete invCurrentFFTObject;
	delete newHist;
        delete value_re;
	delete value_im;
}

void processNtuple(std::string inputFileName) {

  Analyze ana(inputFileName);
  ana.doAnalysis();

  return;
}

int main(int argc, char *argv[]){
  if(argc!=2){
    cout<<"Usage: processNtuple [inputFilename]"<<endl;
    return 0;
  }

  std::string inputFileName = argv[1];
  std::cout << "inputFileName " << inputFileName << std::endl;

  //define ROOT application object
  theApp = new TApplication("App", &argc, argv);
  processNtuple(inputFileName); 

  //return 1;
  gSystem->Exit(0);
}
