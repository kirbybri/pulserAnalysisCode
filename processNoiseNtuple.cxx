//compile independently with: g++ -std=c++11 -o processNoiseNtuple processNoiseNtuple.cxx `root-config --cflags --glibs`
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
#include "TProfile2D.h"
#include "TString.h"
#include "TCanvas.h"
#include <TGraphErrors.h>
#include "TMultiGraph.h"
#include "TImage.h"
#include "TSystem.h"
#include "TVirtualFFT.h"

using namespace std;

//global TApplication object declared here for simplicity
TApplication *theApp;

class Analyze {
	public:
	Analyze(std::string inputFileName);
	void doAnalysis();
	void analyzeChannel();
	void analyzeEvent();
	void summarizeNoise();
	void corrWave();
	void correctGroup( int s, int baseCh, std::vector<int> ch_v);

	//file pointers
	TFile* inputFile;
	TFile *gOut;

	//ROI tr_rawdata variables
	TTree *tr_rawdata;
	std::vector<short> *adc_v = 0;
	unsigned short run,subrun, event, chan;

	//Constants
	const int numChan = 2048;// 35t
	const int maxNumBin = 1000;
	const double SAMP_PERIOD = 0.5;//us

	TGraph *gTest;
	TCanvas* c0;
	TGraph *gChannel;
	int eventNum;

	//histograms
	TH2F *hRawSampVsChan;
	TH2F *hWaveNoStuckSampVsChan;
	TH2F *hCorrWaveSampVsChan;
	TH2F *hRawVsChan;
	TH2F *hWaveNoStuckVsChan;
	TH2F *hCorrWaveVsChan;
	TH2F *hMeanVsChan;
	TProfile *pMeanVsChan;
	TH2F *hRmsVsChan;
	TProfile *pRmsVsChan;
	TProfile2D *pFFTVsChan;
	TH1F *hEvMeanVsChan;

	TProfile *pRawRMSVsChan;
	TProfile *pWaveNoStuckRMSVsChan;
	TProfile *pCorrWaveRMSVsChan;
};

Analyze::Analyze(std::string inputFileName){

	//get input file
	if( inputFileName.empty() ){
		std::cout << "Error invalid file name" << std::endl;
		gSystem->Exit(0);
	}

	inputFile = new TFile(inputFileName.c_str());
	if (inputFile->IsZombie()) {
		std::cout << "Error opening input file" << std::endl;
		gSystem->Exit(0);
	}

	if( !inputFile ){
		std::cout << "Error opening input file" << std::endl;
		gSystem->Exit(0);
	}

	//initialize tr_rawdata branches
  	tr_rawdata = (TTree*) inputFile->Get("GetWF/_wf_tree");
  	if( !tr_rawdata ){
		std::cout << "Error opening input file tree" << std::endl;
		gSystem->Exit(0);
  	}
	tr_rawdata->SetBranchAddress("_run", &run);
  	tr_rawdata->SetBranchAddress("_chan", &chan);
  	tr_rawdata->SetBranchAddress("_event", &event);
  	tr_rawdata->SetBranchAddress("adc_v",&adc_v);

	//make output file
  	std::string outputFileName = "output_processNoiseNtuple_" + inputFileName;
  	gOut = new TFile(outputFileName.c_str() , "RECREATE");

  	//initialize canvas
  	c0 = new TCanvas("c0", "c0",1400,800);

  	//initialize graphs + profiles
  	gTest = new TGraph();
  	gChannel = new TGraph();

  	//histograms + grraphs
	hRawSampVsChan = new TH2F("hRawSampVsChan","",numChan,0-0.5,numChan-0.5,4096,-0.5,4096-0.5);
	hWaveNoStuckSampVsChan = new TH2F("hWaveNoStuckSampVsChan","",numChan,0-0.5,numChan-0.5,4100,-2050-0.5,2050-0.5);
	hCorrWaveSampVsChan = new TH2F("hCorrWaveSampVsChan","",numChan,0-0.5,numChan-0.5,4100,-2050-0.5,2050-0.5);
  	hRawVsChan = new TH2F("hRawVsChan","",numChan,0-0.5,numChan-0.5,maxNumBin,0-0.5,maxNumBin-0.5);
	hWaveNoStuckVsChan = new TH2F("hWaveNoStuckVsChan","",numChan,0-0.5,numChan-0.5,maxNumBin,0-0.5,maxNumBin-0.5);
	hCorrWaveVsChan = new TH2F("hCorrWaveVsChan","",numChan,0-0.5,numChan-0.5,maxNumBin,0-0.5,maxNumBin-0.5);
  	hMeanVsChan = new TH2F("hMeanVsChan","",numChan,0-0.5,numChan-0.5,4096,-0.5,4096-0.5);
	pMeanVsChan = new TProfile("pMeanVsChan","",numChan,0-0.5,numChan-0.5);
  	hRmsVsChan = new TH2F("hRmsVsChan","",numChan,0-0.5,numChan-0.5,300,0,300.);
  	pRmsVsChan = new TProfile("pRmsVsChan","",numChan,0-0.5,numChan-0.5);
	pFFTVsChan = new TProfile2D("pFFTVsChan","",numChan,-0.5,numChan-0.5,500,0,1);
	hEvMeanVsChan = new TH1F("hEvMeanVsChan","",numChan,0-0.5,numChan-0.5);

	pRawRMSVsChan = new TProfile("pRawRMSVsChan","",numChan,0-0.5,numChan-0.5);
	pWaveNoStuckRMSVsChan = new TProfile("pWaveNoStuckRMSVsChan","",numChan,0-0.5,numChan-0.5);
	pCorrWaveRMSVsChan = new TProfile("pCorrWaveRMSVsChan","",numChan,0-0.5,numChan-0.5);
}

void Analyze::doAnalysis(){
  	//loop over tr_rawdata entries
  	Long64_t nEntries(tr_rawdata->GetEntries());

	tr_rawdata->GetEntry(0);
	int currEvent = event;
	int numEvent = 0;
	//loop over puls	e waveform
	for(Long64_t entry(0); entry<nEntries; ++entry) { 
		tr_rawdata->GetEntry(entry);
		//if( event != 11 && event != 42 && event != 56 && event != 60 )
		//	continue;
		//detect new event
		if( event != currEvent ){
			//if( numEvent > 1 )			
			//	break;
			//analyze previous event
			analyzeEvent();
	
			currEvent = event;
			numEvent++;
			std::cout << "Event " << currEvent << std::endl;
			
			hRawVsChan->Reset();
			hWaveNoStuckVsChan->Reset();
			hCorrWaveVsChan->Reset();
			hEvMeanVsChan->Reset();
		}
		//analyze individual channel
    		analyzeChannel();
  	}//entries
	//analyze final event
	analyzeEvent();

	//summarize noise
	summarizeNoise();

    	gOut->Cd("");
	hRawSampVsChan->Write();
	hWaveNoStuckSampVsChan->Write();
	hCorrWaveSampVsChan->Write();
  	hRawVsChan->Write();
	hWaveNoStuckVsChan->Write();
	hCorrWaveVsChan->Write();
  	hMeanVsChan->Write();
	pMeanVsChan->Write();
  	hRmsVsChan->Write();
  	pRmsVsChan->Write();
	pFFTVsChan->Write();
	pRawRMSVsChan->Write();
	pWaveNoStuckRMSVsChan->Write();
	pCorrWaveRMSVsChan->Write();
  	gOut->Close();
}

void Analyze::analyzeChannel(){

	//require at least certain number of samples
	if( adc_v->size() <  maxNumBin )
		return;

	//calculate mean
	double mean = 0;
	int count = 0;
	for( int s = 0 ; s < adc_v->size() ; s++ ){
		short value = adc_v->at(s);
		if( (value & 0x3F ) == 0x0 || (value & 0x3F ) == 0x3F || (value & 0x3F ) == 0x1 )
			continue;
		if( value < 10 ) continue;
		mean += value;
		count++;
	}
	if( count > 0 )
		mean = mean / (double) count;

	//calculate rms
	double rms = 0;
	count = 0;
	for( int s = 0 ; s < adc_v->size() ; s++ ){
		short value = adc_v->at(s);
		if( (value & 0x3F ) == 0x0 || (value & 0x3F ) == 0x3F || (value & 0x3F ) == 0x1 )
			continue;
		if( value < 10 ) continue;
		rms += (value-mean)*(value-mean);
		count++;
	}	
	if( count > 1 )
		rms = TMath::Sqrt( rms / (double)( count - 1 ) );

	//do channel selection here
	/*
	if( mean > 850 || mean < 650 )
		return;
	if( rms > 50 )
		return;
	if( chan >= 1528 && chan <= 1534 )
		return;
	if( chan >= 1912 && chan <= 1918 )
		return;
	if( chan == 172 ) return;
	*/

	//fill channel waveform hists
	for( int s = 0 ; s < adc_v->size() ; s++ ){
		short value = adc_v->at(s);

		hRawSampVsChan->Fill( chan, value);
		hRawVsChan->SetBinContent( chan+1, s+1 , value);
		if( (value & 0x3F ) == 0x0 || (value & 0x3F ) == 0x3F || (value & 0x3F ) == 0x1 )
			continue;
		if( value < 10 ) continue;
		hWaveNoStuckSampVsChan->Fill( chan, value-mean);
		hWaveNoStuckVsChan->SetBinContent( chan+1, s+1, value-mean);
	}

	//do FFT
	std::string histName = "hData_" + std::to_string(chan);
	int numBins = adc_v->size();
    	TH1F *hData = new TH1F(histName.c_str(),"",numBins,0,numBins*0.5 );
    	for(int i = 0 ; i < adc_v->size() ; i++ ){
		if( adc_v->at(i) > 10 )
			hData->SetBinContent(i+1, adc_v->at(i) - mean );
		else
			hData->SetBinContent(i+1, 0 );
    	}

	TH1F *hFftData = new TH1F("hFftData","",numBins,0,numBins );
        hData->FFT(hFftData,"MAG");
    	for(int i = 1 ; i < hFftData->GetNbinsX() ; i++ ){
		double freq = 2.* i / (double) hFftData->GetNbinsX() ;
		pFFTVsChan->Fill( chan, freq, hFftData->GetBinContent(i+1) );
    	}

	delete hData;
	delete hFftData;

	//channel specific
	hMeanVsChan->Fill( chan, mean );
	pMeanVsChan->Fill( chan, mean );
	hRmsVsChan->Fill(chan, rms);
	pRmsVsChan->Fill(chan, rms);
	
	//event specific
	hEvMeanVsChan->Fill( chan, mean);
}

void Analyze::analyzeEvent(){
	corrWave();

	if( 0 ){
		c0->Clear();	
		hCorrWaveVsChan->Draw("COLZ");
		c0->Update();
		char ct;
		std::cin >> ct;
	}
}

void Analyze::corrWave(){

	hCorrWaveVsChan->Reset();

	//define set of induciton and collection channels in each regulator group
	std::vector<int> indCh1;
	std::vector<int> indCh2;
	std::vector<int> indCh3;
	std::vector<int> indCh4;
	std::vector<int> colCh1;
	std::vector<int> colCh2;
	std::vector<int> colCh3;
	std::vector<int> colCh4;
	for(int i = 0 ; i < 64 ; i++){
		if( i < 32 || i == 32 || i == 47 || i == 48 || i == 63 ){
			if( i < 8 )
				indCh1.push_back(i);
			if( i >= 8 && i < 16 )
				indCh2.push_back(i);
			if( i >= 16 && i < 24 )
				indCh3.push_back(i);
			if( i >= 24 )
				indCh4.push_back(i);
		}
		else{
			if( i < 40 )
				colCh1.push_back(i);
			if( i >= 40 && i < 48 )
				colCh2.push_back(i);
			if( i >= 48 && i < 56 )
				colCh3.push_back(i);
			if( i >= 56  )
				colCh4.push_back(i);
		}
	}
	
	//loop through time slices
	for(int s = 0 ; s < hRawVsChan->GetNbinsY() ; s++ ){
		//loop through regulator groups
		for(int g = 0 ; g < 16*2 ; g++){
			int baseCh = g*64;
			//get induction plane correction
			correctGroup(s,g*64, indCh1);
			correctGroup(s,g*64, indCh2);
			correctGroup(s,g*64, indCh3);
			correctGroup(s,g*64, indCh4);
			//get collection plane correction
			correctGroup(s,g*64, colCh1);
			correctGroup(s,g*64, colCh2);
			correctGroup(s,g*64, colCh3);
			correctGroup(s,g*64, colCh4);
		}//end loop over regulator groups
	}//end loop over samples

	return;
}

//derive correction factors - require raw adc waveform and pedestal for each channel
void Analyze::correctGroup( int s, int baseCh, std::vector<int> ch_v){
	std::vector<Double_t> corrVals;
	corrVals.clear();
	for(int c = 0 ; c < ch_v.size() ; c++){
		int ch = baseCh + ch_v.at(c);
		int adc = hRawVsChan->GetBinContent(ch+1,s+1);
		if( (adc & 0x3F ) == 0x0 || (adc & 0x3F ) == 0x3F || (adc & 0x3F ) == 0x1 )
			continue;
		if( adc < 10  ) //skip "sample dropping" problem
			continue;
		double mean = hEvMeanVsChan->GetBinContent(ch+1);
		if( mean < 10 )
			continue;
		corrVals.push_back(adc - mean);
	}
			
	//derive correction
	int corrValSize = corrVals.size();
    	sort(corrVals.begin(),corrVals.end());
	double correction = 0;
    	if(corrValSize < 2)
      		correction = 0.0;
    	else if((corrValSize % 2) == 0)
      		correction = (corrVals[corrValSize/2] + corrVals[(corrValSize/2)-1])/2.0;
   	else
      		correction = corrVals[(corrValSize-1)/2];
	
	for(int c = 0 ; c < ch_v.size() ; c++){
		int ch = baseCh + ch_v.at(c);
		int adc = hRawVsChan->GetBinContent(ch+1,s+1);
		hCorrWaveVsChan->SetBinContent(ch+1,s+1, 0);
		if( (adc & 0x3F ) == 0x0 || (adc & 0x3F ) == 0x3F || (adc & 0x3F ) == 0x1 )
			continue;
		if( adc < 10  ) //skip "sample dropping" problem
			continue;
		double mean = hEvMeanVsChan->GetBinContent(ch+1);
		if( mean < 10 )
			continue;
		double newAdc = adc - correction;
		hCorrWaveSampVsChan->Fill(ch, newAdc - mean);
		hCorrWaveVsChan->SetBinContent(ch+1,s+1, newAdc);
		//hCorrWaveVsChan->SetBinContent(ch+1,s+1, adc);
	}
}

void Analyze::summarizeNoise(){

	for(int ch = 0 ; ch < hRawSampVsChan->GetNbinsX() ; ch++ ){
		char name[200];
		memset(name,0,sizeof(char)*100 );
	        sprintf(name,"hCh_%.3i",ch);
		TH1D *hChan = hRawSampVsChan->ProjectionY(name,ch+1,ch+1);
		hChan->GetXaxis()->SetRangeUser(0,4095);
		double rms = hChan->GetRMS();
		pRawRMSVsChan->Fill(ch,rms);
	}
	for(int ch = 0 ; ch < hWaveNoStuckSampVsChan->GetNbinsX() ; ch++ ){
		char name[200];
		memset(name,0,sizeof(char)*100 );
	        sprintf(name,"hCh_%.3i",ch);
		TH1D *hChan = hWaveNoStuckSampVsChan->ProjectionY(name,ch+1,ch+1);
		hChan->GetXaxis()->SetRangeUser(-2099,2099);
		double rms = hChan->GetRMS();
		pWaveNoStuckRMSVsChan->Fill(ch,rms);
	}
	for(int ch = 0 ; ch < hCorrWaveSampVsChan->GetNbinsX() ; ch++ ){
		char name[200];
		memset(name,0,sizeof(char)*100 );
	        sprintf(name,"hCh_%.3i",ch);
		TH1D *hChan = hCorrWaveSampVsChan->ProjectionY(name,ch+1,ch+1);
		hChan->GetXaxis()->SetRangeUser(-2099,2099);
		double rms = hChan->GetRMS();
		pCorrWaveRMSVsChan->Fill(ch,rms);
	}
}

void processNoiseNtuple(std::string inputFileName) {

  Analyze ana(inputFileName);
  ana.doAnalysis();

  return;
}

int main(int argc, char *argv[]){
  if(argc!=2){
    cout<<"Usage: processNoiseNtuple [inputFilename]"<<endl;
    return 0;
  }

  std::string inputFileName = argv[1];
  std::cout << "inputFileName " << inputFileName << std::endl;

  //define ROOT application object
  theApp = new TApplication("App", &argc, argv);
  processNoiseNtuple(inputFileName); 

  //return 1;
  gSystem->Exit(0);
}
