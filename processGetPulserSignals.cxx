//compile independently with: g++ -std=c++11 -o processGetPulserSignals processGetPulserSignals.cxx `root-config --cflags --glibs` -lMinuit
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
using namespace std;

#include "TROOT.h"
#include "TMath.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TProfile2D.h"
#include "TString.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TF1.h"
#include "TMarker.h"

#include "constants.hxx"
#include "FitFeElecResponse_multiPulse.hxx"
#include "FitFeElecResponse_multiPulse.cxx"

using namespace std;

//global TApplication object declared here for simplicity
TApplication *theApp;

class Analyze {
	public:
	Analyze(std::string inputFileName);
	int processFileName(std::string inputFileName, std::string &baseFileName);

	void doAnalysis();
	void measurePulserTimes();

	void analyzeChannel(unsigned int chan);
	void initializeTree( unsigned int chan );
	void dumpTreeInfo();
	void drawPulse();
	void drawFit();

	void getChannelBaseline();
	void measurePulseDist();
	void measureChannelDist();
	bool getPulseStatus();
	bool getChannelStatus(unsigned int chan);

	void getAvgPulseShape();
	void doSimpleMeasurement();

	void getPulses();
	void doChannelFit(unsigned int chan);

	void getAvgFitResidual();

	void getQRes();
	void doSinglePulseFit(unsigned int eventNum, unsigned int firstNum);
	void drawFit_singlePulse();
	void fillOutputTree();

	//Files
	TFile* inputFile;
	TFile *gOut;

	//data objects
	TCanvas* c0;
	TMarker *mMark0;
	TGraph *gCh;
	TGraph *gFit;

	//input tree
  	unsigned short fRun, fSubrun, fChan, fNum, fSample, fMaxValue, fMinValue, fMaxSample, fMinSample;
	unsigned int fEvent;
  	float fStartSample, fMean, fRms, fSumVal;
	std::vector<unsigned short> *fWf = 0;
  	TTree *tr_fit;

	//pulser tree
	unsigned short fRun_pulser, fSubrun_pulser;
	unsigned int fEvent_pulser;
	std::vector<double> *fPulseStartSamples_pulser = 0;
  	TTree *tr_pulser;

	//Fit objects
	FitFeElecResponse_multiPulse *fitResponse;
	FitFeElecResponse_multiPulse *fitResponse_singlePulse;
	FeElecResponse *fitSig;//drawing purposes only
	
	//constants
	int constant_numChan = 8256;
	double constant_baseFitRange = 2.;
	double constant_pulseFitRange = 4.0;
	unsigned int constant_maxNumberEventsProcessed = 50;
	bool constant_showFitInfo = 0;
	int constant_maxNumberChannels = 10000;
	bool constant_printResults = 1;
	bool constant_doDrawAvgWaveform = 0;
	double constant_ampToHeightFactor = 0.0988165;
	bool constant_fixStartTimes = 1;

	//variables
	double avgPulserPeriod;
	std::vector<double> evPulserStartSamplesErrors;
	std::vector<double> evPulserStartSamples;

	double chAvgMean;
	double chAvgMaxValue;
	double chAvgRms;
	double chAvgSum;
	double chAvgRiseTime;

	double chBase;
	double chRms;

	//histograms
	TH1F *hMean;
  	TH1F *hRms;
	TH1F *hBaseVsChan;
	TH1F *hRmsVsChan;

	TH1F *hPulseRms;
	TH1F *hPulseRmsU;
	TH1F *hPulseRmsV;
	TH1F *hPulseRmsY;
  	TH1F *hPulseHeight;
  	TH1F *hPulseMaxValue;
	TH1F *hPulseMaxPosition;

  	TH1F *hChAvgPulseHeight;

	TH1F *hSamp;
	TH1F *hMaxValue;

	TProfile *pChStatusVsChan;
	TProfile *pPulseStatusVsChan;

	TProfile2D *pAvgSignalVsChan;
	TH1F *hAvgSignalPulseHeightVsChan;

	TH1F *hWidth;
	TH1F *hSum;
	TProfile *pPulseHeightVsChan;
	TProfile *pWidthVsChan;
	TProfile *pSumVsChan;

	TH1F *hFitPulseHeightVsChan;
	TH1F *hFitAmpVsChan;
	TH1F *hFitShapeVsChan;
	TH1F *hFitStartVsChan;
	TProfile2D *pAvgFitSignalVsChan;
	TProfile2D *pAvgFitResidualVsChan;

	TProfile *pSinglePulseAmpVsChan;

	//output event time tree
  	TTree *tr_eventTime;

	//output pulse tree
  	unsigned short oRun, oSubrun, oChan;
	int oStatus;
  	float oBase, oShape, oAmp, oBaseErr, oShapeErr, oAmpErr;
	std::vector<double> oStartTimes;
  	TTree *tr_out;
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

	//make output file
  	std::string outputFileName = "output_processGetPulserSignals.root";
	//if( processFileName( inputFileName, outputFileName ) )
	//	outputFileName = "output_processGetPulserSignals_" + outputFileName;
  	gOut = new TFile(outputFileName.c_str() , "RECREATE");

  	//initialize canvas
  	c0 = new TCanvas("c0", "c0",1400,800);

	mMark0 = new TMarker();
	mMark0->SetMarkerStyle(20);
	mMark0->SetMarkerSize(2);
	mMark0->SetMarkerColor(kRed);

	//initialize graphs
	gCh = new TGraph();
	gFit = new TGraph();

	//initialize histograms

	//fit objects	
	fitSig = new FeElecResponse();
	fitResponse = new FitFeElecResponse_multiPulse();
	fitResponse_singlePulse = new FitFeElecResponse_multiPulse();

	//initialize histograms
	hMean = new TH1F("hMean","",4100,0-0.5,4100-0.5);
	hRms = new TH1F("hRms","",2000,0,200);
	hBaseVsChan = new TH1F("hBaseVsChan","",8256,0-0.5,8256-0.5);
	hRmsVsChan = new TH1F("hRmsVsChan","",8256,0-0.5,8256-0.5);

	hPulseRms = new TH1F("hPulseRms","",1000,0,100);
	hPulseRmsU = new TH1F("hPulseRmsU","",1000,0,100);
	hPulseRmsV = new TH1F("hPulseRmsV","",1000,0,100);
	hPulseRmsY = new TH1F("hPulseRmsY","",1000,0,100);
  	hPulseHeight = new TH1F("hPulseHeight","",4100,0,4100);
	hPulseMaxValue = new TH1F("hPulseMaxValue","",4100,0-0.5,4100-0.5);
	hPulseMaxPosition = new TH1F("hPulseMaxPosition","",300,-10,20);

  	hChAvgPulseHeight = new TH1F("hChAvgPulseHeight","",4100,0,4100);

	hSamp = new TH1F("hSamp","",4100,0-0.5,4100-0.5);
	hMaxValue = new TH1F("hMaxValue","",4100,0-0.5,4100-0.5);

	pChStatusVsChan = new TProfile("pChStatusVsChan","",8256,0-0.5,8256-0.5);
	pPulseStatusVsChan = new TProfile("pPulseStatusVsChan","",8256,0-0.5,8256-0.5);

	pAvgSignalVsChan = new TProfile2D("pAvgSignalVsChan","",8256,0-0.5,8256-0.5,500,-10,40);

	hWidth = new TH1F("hWidth","",1000,0,10);
	hSum = new TH1F("hSum","",2000,0,20000);
	pPulseHeightVsChan = new TProfile("pPulseHeightVsChan","",8256,0-0.5,8256-0.5);
	pWidthVsChan = new TProfile("pWidthVsChan","",8256,0-0.5,8256-0.5);
	pSumVsChan = new TProfile("pSumVsChan","",8256,0-0.5,8256-0.5);

	hFitPulseHeightVsChan = new TH1F("hFitPulseHeightVsChan","",8256,0-0.5,8256-0.5);
	hFitAmpVsChan = new TH1F("hFitAmpVsChan","",8256,0-0.5,8256-0.5);
	hFitShapeVsChan = new TH1F("hFitShapeVsChan","",8256,0-0.5,8256-0.5);
	hFitStartVsChan = new TH1F("hFitStartVsChan","",8256,0-0.5,8256-0.5);
	//pAvgFitSignalVsChan = new TProfile2D("pAvgFitSignalVsChan","",8256,0-0.5,8256-0.5,500,-10,40);
	pAvgFitResidualVsChan = new TProfile2D("pAvgFitResidualVsChan","",8256,0-0.5,8256-0.5,500,-10,40);

	pSinglePulseAmpVsChan = new TProfile("pSinglePulseAmpVsChan","",8256,0-0.5,8256-0.5);
	hAvgSignalPulseHeightVsChan = new TH1F("hAvgSignalPulseHeightVsChan","",8256,0-0.5,8256-0.5);

	//make output tree(s)
	tr_out = new TTree( "tr_out" , "tr_out" );
  	tr_out->Branch("run", &oRun, "run/s");
  	tr_out->Branch("subrun", &oSubrun, "subrun/s");
  	tr_out->Branch("chan", &oChan, "chan/s");
	tr_out->Branch("status", &oStatus, "status/I");
  	tr_out->Branch("base", &oBase, "base/F");
  	tr_out->Branch("shape", &oShape, "shape/F");
  	tr_out->Branch("amp", &oAmp, "amp/F");
  	tr_out->Branch("baseErr", &oBaseErr, "baseErr/F");
  	tr_out->Branch("shapeErr", &oShapeErr, "shapeErr/F");
  	tr_out->Branch("ampErr", &oAmpErr, "ampErr/F");
	tr_out->Branch("startTimes", "vector<double>", &oStartTimes);
}

int Analyze::processFileName(std::string inputFileName, std::string &baseFileName){
        //check if filename is empty
        if( inputFileName.size() == 0 ){
                std::cout << "processFileName : Invalid filename " << std::endl;
                return 0;
        }

        //remove path from name
        size_t pos = 0;
        std::string delimiter = "/";
        while ((pos = inputFileName.find(delimiter)) != std::string::npos)
                inputFileName.erase(0, pos + delimiter.length());

	if( inputFileName.size() == 0 ){
                std::cout << "processFileName : Invalid filename " << std::endl;
                return 0;
        }

        //replace / with _
        std::replace( inputFileName.begin(), inputFileName.end(), '/', '_'); // replace all 'x' to 'y'
        std::replace( inputFileName.begin(), inputFileName.end(), '-', '_'); // replace all 'x' to 'y'

	baseFileName = inputFileName;
	
	return 1;
}

void Analyze::doAnalysis(){

	measurePulserTimes();
	if( avgPulserPeriod <= 0 || evPulserStartSamples.size() == 0 )
		return;

	/*
	analyzeChannel(1200);
	analyzeChannel(3600);
	analyzeChannel(6200);
	*/
	
	for(unsigned int ch = 0 ; ch < constant_numChan ; ch++ ){
	//for(unsigned int ch = 2400 ; ch < 2450 ; ch++ ){
		if( ch % 100 == 0 ) std::cout << "Channel " << ch << std::endl;

		std::clock_t start;
	        start = std::clock();
	
		analyzeChannel(ch);

		//std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC ) << " s" << std::endl;

		if( ch > constant_maxNumberChannels )
			break; 
	}
	

    	gOut->Cd("");
	hMean->Write();
	hRms->Write();
	hBaseVsChan->Write();
	hRmsVsChan->Write();

	hPulseRms->Write();
	hPulseRmsU->Write();
	hPulseRmsV->Write();
	hPulseRmsY->Write();
	hPulseHeight->Write();
	hPulseMaxValue->Write();
	hPulseMaxPosition->Write();

	hChAvgPulseHeight->Write();

	pChStatusVsChan->Write();
	pPulseStatusVsChan->Write();

	pAvgSignalVsChan->Write();
	hAvgSignalPulseHeightVsChan->Write();
	
	pPulseHeightVsChan->Write();
	pWidthVsChan->Write();
	pSumVsChan->Write();

	hFitPulseHeightVsChan->Write();
	hFitAmpVsChan->Write();
	hFitShapeVsChan->Write();
	hFitStartVsChan->Write();
	//pAvgFitSignalVsChan->Write();
	pAvgFitResidualVsChan->Write();

	pSinglePulseAmpVsChan->Write();

	tr_out->Write();
  	gOut->Close();
}

void Analyze::measurePulserTimes(){
	avgPulserPeriod = -1;
	evPulserStartSamples.clear();
	evPulserStartSamplesErrors.clear();

	//initialize tree branches
	std::string title = "GetPulserSignals/tr_pulserStartSamples";
	tr_pulser = (TTree*) inputFile->Get( title.c_str() );
	if( !tr_pulser ){
		std::cout << "Error opening input file tree, exiting" << std::endl;
		gSystem->Exit(0);
  	}

  	tr_pulser->SetBranchAddress("run", &fRun_pulser);
  	tr_pulser->SetBranchAddress("subrun", &fSubrun_pulser);
  	tr_pulser->SetBranchAddress("event", &fEvent_pulser);
	tr_pulser->SetBranchAddress("pulserStartSamples", &fPulseStartSamples_pulser);

	unsigned int numEntries = tr_pulser->GetEntries();
	if( numEntries == 0 ) return;
	tr_pulser->GetEntry(0);

	TF1 *func = new TF1("fit","pol1",-0.5,100);
	avgPulserPeriod = 0;
	int count = 0;

	for(unsigned int entry = 0 ; entry<numEntries; entry++) {
		tr_pulser->GetEntry(entry);
		gCh->Set(0);
		for(unsigned p = 0 ; p < fPulseStartSamples_pulser->size() ; p++ ){
			//std::cout << "\t" << fEvent_pulser << "\t" << p << "\t" << fPulseStartSamples_pulser->at(p) << std::endl;
			gCh->SetPoint( gCh->GetN(), p , fPulseStartSamples_pulser->at(p) );
		}

		if( gCh->GetN() == 0 )
			continue;
		if( gCh->GetN() == 1 ){
			avgPulserPeriod = 0;
			continue;
		}

		func->SetRange(-0.5, fPulseStartSamples_pulser->size() + 0.5 );

		gCh->Fit("fit","QR");
		if( func->GetParError(0) < 0.01 ){
			avgPulserPeriod += func->GetParameter(1);
			count++;
		}
		//fitPulserStartSamples.push_back( func->GetParameter(0) );

		if(0){
			c0->Clear();
			gCh->SetMarkerStyle(21);
			gCh->Draw("AP");
			c0->Update();
			char ct;
			std::cin >> ct;
		}
		
	}
	if( count <= 0 )
		return;

	avgPulserPeriod = avgPulserPeriod / double( count );

	TF1 *newfunc = new TF1("newfit","pol1",-0.5,100);
	newfunc->FixParameter(1, avgPulserPeriod );

	for(unsigned int entry = 0 ; entry<numEntries; entry++) {
		tr_pulser->GetEntry(entry);
		gCh->Set(0);
		for(unsigned p = 0 ; p < fPulseStartSamples_pulser->size() ; p++ ){
			//std::cout << "\t" << fEvent_pulser << "\t" << p << "\t" << fPulseStartSamples_pulser->at(p) << std::endl;
			gCh->SetPoint( gCh->GetN(), p , fPulseStartSamples_pulser->at(p) );
		}

		if( gCh->GetN() == 0 )
			continue;
		if( gCh->GetN() == 1 ){
			avgPulserPeriod = 0;
			continue;
		}

		newfunc->SetRange(-0.5, fPulseStartSamples_pulser->size() + 0.5 );
		gCh->Fit("newfit","QBR");
		evPulserStartSamples.push_back( newfunc->GetParameter(0) );
		evPulserStartSamplesErrors.push_back( newfunc->GetParError(0) );

		if(0){
			std::string title = "Event " + to_string(fEvent_pulser);
			gCh->SetTitle( title.c_str() );
			c0->Clear();
			gCh->SetMarkerStyle(21);
			gCh->Draw("AP");
			c0->Update();
			char ct;
			std::cin >> ct;
		}
	}

	//gCh->GetFunction("fit")->Delete();
	gCh->GetFunction("newfit")->Delete();

	delete func;
	delete newfunc;
}

void Analyze::analyzeChannel(unsigned int chan){

	if( chan >= constant_numChan ) return;
	initializeTree( chan );

	unsigned int numEntries = tr_fit->GetEntries();
	if( numEntries == 0 ) return;
	//dumpTreeInfo();

	getChannelBaseline();
	if( chBase <= 0 || chRms < 0 )
		return;

	measurePulseDist();
	measureChannelDist();
	getAvgPulseShape();

	bool channelStatus = getChannelStatus(chan);
	pChStatusVsChan->Fill(chan, channelStatus);
	if( channelStatus == 0 )
		return;

	doSimpleMeasurement();

	getPulses();
	if( fitResponse->numPulses == 0 )
		return;

	doChannelFit(chan);

	getAvgFitResidual();

	//test outcome of fit here

	fillOutputTree();

	//getQRes();
	
	//char ct;
	//std::cin >> ct;

	return;
}

void Analyze::initializeTree( unsigned int chan ){
	//initialize tree branches
	std::string title = "GetPulserSignals/tr_" + to_string( chan );
	tr_fit = (TTree*) inputFile->Get( title.c_str() );
	if( !tr_fit ){
		std::cout << "Error opening input file tree, exiting" << std::endl;
		gSystem->Exit(0);
  	}

  	tr_fit->SetBranchAddress("run", &fRun);
  	tr_fit->SetBranchAddress("subrun", &fSubrun);
  	tr_fit->SetBranchAddress("event", &fEvent);
  	tr_fit->SetBranchAddress("chan", &fChan);
	tr_fit->SetBranchAddress("num", &fNum);
  	tr_fit->SetBranchAddress("pulsestartsample", &fStartSample);
  	tr_fit->SetBranchAddress("firstsample", &fSample);
  	tr_fit->SetBranchAddress("maxvalue", &fMaxValue);
  	tr_fit->SetBranchAddress("minvalue", &fMinValue);
  	tr_fit->SetBranchAddress("maxsample", &fMaxSample);
  	tr_fit->SetBranchAddress("minsample", &fMinSample);
  	tr_fit->SetBranchAddress("mean", &fMean);
  	tr_fit->SetBranchAddress("rms", &fRms);
  	tr_fit->SetBranchAddress("sum", &fSumVal);
 	tr_fit->SetBranchAddress("wf", &fWf);
}

void Analyze::dumpTreeInfo(){
	unsigned int numEntries = tr_fit->GetEntries();
	if( numEntries == 0 ) return;
	tr_fit->GetEntry(0);

	//check if a good channel
	for(unsigned int entry = 0 ; entry<numEntries; entry++) { 
		tr_fit->GetEntry(entry);
		//tr_fit->Show();
		std::cout << fEvent << "\t" << fNum << std::endl;
		//drawPulse();
	}
	std::cout << "Done" << std::endl;
	char ct;
	std::cin >> ct;
	return;
}

void Analyze::getChannelBaseline(){
	chBase = -1;
	chRms = -1;

	unsigned int numEntries = tr_fit->GetEntries();
	if( numEntries == 0 ) return;
	tr_fit->GetEntry(0);

	hSamp->Reset();
	for(unsigned int entry = 0 ; entry<numEntries; entry++) { 
		tr_fit->GetEntry(entry);

		gCh->Set(0);

		for( int s = 0 ; s < fWf->size() ; s++ ){
			double sample = fSample + s;
			if( sample < fStartSample - 4 ){
				gCh->SetPoint(gCh->GetN() , sample , fWf->at(s) );
				hSamp->Fill( fWf->at(s) );
			}
		}

		if(0){
			c0->Clear();
			std::string title = "Get Baseline: Channel " + to_string( fChan ) + ", Event " + to_string(fEvent);
			gCh->SetTitle( title.c_str() );
			gCh->GetXaxis()->SetTitle("Sample Number");
			gCh->GetYaxis()->SetTitle("Sample Value (ADC counts)");
			gCh->SetMarkerStyle(21);
			gCh->Draw("AP");
			c0->Update();	
			//char cr;
			//std::cin >> cr;
		}
	}

	//extract baseline mean from pedestal sample distribution
	hSamp->GetXaxis()->SetRangeUser(0.5,4094.5);
	double histMean = hSamp->GetMean();
	double histPeak = hSamp->GetBinCenter( hSamp->GetMaximumBin() );
	hSamp->GetXaxis()->SetRangeUser(histPeak-50,histPeak+50.);
	double histRms = hSamp->GetRMS();
	if( histRms < 1 )
		histRms = 1;
	double lowerBound = histPeak-3.5*histRms;
	if( lowerBound < 0.5 )
		lowerBound = 0.5;
	double upperBound = histPeak+3.5*histRms;
	if( upperBound > 4094.5 )
		upperBound = 4094.5;
	hSamp->GetXaxis()->SetRangeUser(lowerBound,upperBound);
	
	//set class chBase variable to the pedestal mean
	chBase = hSamp->GetMean();
	chRms = hSamp->GetRMS();
	hMean->Fill(chBase);
	hRms->Fill(chRms);
	hBaseVsChan->SetBinContent(fChan+1, chBase);
	hRmsVsChan->SetBinContent(fChan+1, chRms);

	if( 0 ){
		std::cout << chAvgMean << "\t" << chBase << std::endl;
		c0->Clear();
		std::string title = "Get Baseline: Channel " + to_string( fChan );
		hSamp->SetTitle( title.c_str() );
		hSamp->Draw();
		c0->Update();	
		char cr;
		std::cin >> cr;
	}

	return;
}

void Analyze::measurePulseDist(){
	unsigned int numEntries = tr_fit->GetEntries();
	if( numEntries == 0 ) return;
	tr_fit->GetEntry(0);

	//loop over channel pulses, measure properties
	int currEvent = -1;
	int numEvent = 0;
	for(unsigned int entry = 0 ; entry<numEntries; entry++) { 
		tr_fit->GetEntry(entry);
		if( fEvent != currEvent ){
			currEvent = fEvent;
			numEvent++;
		}

		hPulseRms->Fill(fRms);
		if( fChan >= 0 && fChan < 2400 )
			hPulseRmsU->Fill(fRms);
		if( fChan >= 2400 && fChan < 4800 )
			hPulseRmsV->Fill(fRms);
		if( fChan >= 4800 && fChan < 8256 )
			hPulseRmsY->Fill(fRms);
  		hPulseHeight->Fill(fMaxValue - chBase);
  		hPulseMaxValue->Fill(fMaxValue);
		hPulseMaxPosition->Fill(fMaxSample - fStartSample);
	}
}

bool Analyze::getPulseStatus(){

	//pulse selection here
	bool isGood = 1;
	if( fChan >= 0 && fChan < 2400 && fRms > 12 )
		isGood = 0;
	if( fChan >= 2400 && fChan < 4800 && fRms > 8 )
		isGood = 0;
	if( fChan >= 4800 && fChan < 8256 && fRms > 5 )
		isGood = 0;
	if( fMaxValue > 3900 )
		isGood = 0;
	if( fMaxValue - chBase < 50 )
		isGood = 0;
	if( fMaxSample - fStartSample < 0 || fMaxSample - fStartSample > 8 )
		isGood = 0;

	return isGood;
}

void Analyze::measureChannelDist(){
	chAvgMaxValue = 0;

	unsigned int numEntries = tr_fit->GetEntries();
	if( numEntries == 0 ) return;
	tr_fit->GetEntry(0);

	int count = 0;
	hMaxValue->Reset();
	for(unsigned int entry = 0 ; entry<numEntries; entry++) { 
		tr_fit->GetEntry(entry);

		bool isGood = getPulseStatus();
		pPulseStatusVsChan->Fill(fChan,isGood);
		if( isGood == 0 )
			continue;

		hMaxValue->Fill(fMaxValue);
		chAvgMaxValue += fMaxValue;
		count++;
	}

	double numPulse = double( count );
	if( numPulse <= 0 )
		return;

	hMaxValue->GetXaxis()->SetRangeUser(0.5,4094.5);
	double histMean = hMaxValue->GetMean();
	double histPeak = hMaxValue->GetBinCenter( hMaxValue->GetMaximumBin() );
	hMaxValue->GetXaxis()->SetRangeUser(histMean-50,histMean+50.);

	//chAvgMaxValue = chAvgMaxValue / numPulse;
	chAvgMaxValue = hMaxValue->GetMean();
  	hChAvgPulseHeight->Fill(chAvgMaxValue - chBase);

	return;
}

bool Analyze::getChannelStatus(unsigned int chan){
	//channel selection here
	bool isGood = 1;
	if( chRms > 12 )
		isGood = 0;
	if( chAvgMaxValue - chBase < 50 )
		isGood = 0;
	double goodPulseFrac = pPulseStatusVsChan->GetBinContent(chan+1);
	if( goodPulseFrac < 0.95 )
		isGood = 0;

	return isGood;
}

void Analyze::getAvgPulseShape(){

	unsigned int numEntries = tr_fit->GetEntries();
	if( numEntries == 0 ) return;
	tr_fit->GetEntry(0);

	gCh->Set(0);
	gCh->Clear();

	int currEvent = -1;
	int eventNum = -1;
	int firstNum = 0;
	for(unsigned int entry = 0 ; entry<numEntries; entry++) { 
		tr_fit->GetEntry(entry);
		if( fEvent != currEvent ){
			currEvent = fEvent;
			eventNum++;
			firstNum = fNum;
		}
		
		//if(0 && !getPulseStatus() )
		//	continue;
		if( evPulserStartSamplesErrors.at(eventNum) > 0.1 ) continue;

		double startSample = evPulserStartSamples.at(eventNum) + fNum*avgPulserPeriod;
		for( int s = 0 ; s < fWf->size() ; s++ ){
			if( constant_doDrawAvgWaveform )
				gCh->SetPoint( gCh->GetN(), fSample + s - startSample , fWf->at(s) );
			pAvgSignalVsChan->Fill(fChan, fSample + s - startSample , fWf->at(s) - chBase );
		}
	}

	//get pulse height
	double maxValue = pAvgSignalVsChan->GetBinContent(fChan+1, 1);
	for( int s = 0 ; s < pAvgSignalVsChan->GetNbinsY() ; s++ ){
		if( pAvgSignalVsChan->GetBinContent(fChan+1, s+1) > maxValue )
			maxValue = pAvgSignalVsChan->GetBinContent(fChan+1, s+1);
	}
	hAvgSignalPulseHeightVsChan->SetBinContent(fChan+1, maxValue);
	hAvgSignalPulseHeightVsChan->SetBinError(fChan+1,0);

	if( constant_doDrawAvgWaveform ){

		c0->Clear();

		std::string title = "Avg Pulse: Channel " + to_string( fChan );
		gCh->SetTitle( title.c_str() );
		gCh->GetXaxis()->SetTitle("Sample Number");
		gCh->GetYaxis()->SetTitle("Sample Value (ADC counts)");
		gCh->SetMarkerStyle(1);
		gCh->Draw("AP");
		c0->Update();

		//sleep(0.1);
		char ct;
		std::cin >> ct;
	}

	return;
}

void Analyze::doSimpleMeasurement(){
	unsigned int numEntries = tr_fit->GetEntries();
	if( numEntries == 0 ) return;
	tr_fit->GetEntry(0);

	hWidth->Reset();
	hSum->Reset();

	gCh->Set(0);

	for(unsigned int entry = 0 ; entry<numEntries; entry++) { 
		tr_fit->GetEntry(entry);
		if( getPulseStatus() == 0 ) continue;

		//integral
		double sum = 0;
		for( int s = 0 ; s < fWf->size() ; s++ ){
			//if( fSample + s > fStartSample - 2 )
			if( fSample + s > fStartSample - 2 && fSample + s < fStartSample + 20 )  
				sum += fWf->at(s) - chBase;
		}
		hSum->Fill( sum );
		pSumVsChan->Fill(fChan, sum);
		
		double thresholdSampleRise = -1;
		double thresholdSampleFall = -1;

		double threshold = 0.5*(chAvgMaxValue-chBase) + chBase;
		if( fMaxValue < threshold )
			continue;

		//rising edge
		for( int samp = fMaxSample - 1 ; samp > fMaxSample - 10 ; samp-- ){
			int s = samp - fSample;
			if( s < 0 || s >= fWf->size()) break;
			if( fWf->at(s) <= threshold && fWf->at(s+1) > threshold ){
				double slope = fWf->at(s+1) - fWf->at(s);
				if( slope <= 0 ) break;
				thresholdSampleRise = fSample + s + (threshold - fWf->at(s)) / slope;
				break;
			}
		}
		if( thresholdSampleRise < 0 )
			continue;

		//falling edge
		for( int samp = fMaxSample ; samp < fMaxSample + 10 ; samp++ ){
			int s = samp - fSample;
			if( s < 0 || s >= fWf->size()-1) break;
			if( fWf->at(s) >= threshold && fWf->at(s+1) < threshold ){
				double slope = fWf->at(s+1) - fWf->at(s);
				if( slope >= 0 ) break;
				thresholdSampleFall = fSample + s + (threshold - fWf->at(s)) / slope;
				break;
			}
		}
		if( thresholdSampleFall < 0 )
			continue;

		double width = thresholdSampleFall - thresholdSampleRise;
		hWidth->Fill( width );
		
	}

	double meanWidth = hWidth->GetBinCenter( hWidth->GetMaximumBin() );
	pWidthVsChan->Fill(fChan, meanWidth );

	if(0){
		mMark0->SetX( meanWidth );
		mMark0->SetY( hWidth->GetMaximum()  );

		c0->Clear();
		hWidth->Draw();
		mMark0->Draw("same");
		c0->Update();
		//char ct;
		//std::cin >> ct;
	}

	if(0){
		//extract integral mean from integral distribution
		hSum->GetXaxis()->SetRangeUser(0.5,20000-0.5);
		double histMean = hSum->GetMean();
		double histPeak = hSum->GetBinCenter( hSum->GetMaximumBin() );
		double histRms = hSamp->GetRMS();
		if( histRms < 100 )
			histRms = 100;
		hSum->GetXaxis()->SetRangeUser(histPeak-3*histRms,histPeak+3*histRms);
		std::cout << pSumVsChan->GetBinContent(fChan + 1 ) << std::endl;
		mMark0->SetX( pSumVsChan->GetBinContent(fChan + 1 ) );
		mMark0->SetY( hSum->GetMaximum()  );
		c0->Clear();
		hSum->Draw();
		mMark0->Draw("same");
		c0->Update();
		//char ct;
		//std::cin >> ct;
	}

	pPulseHeightVsChan->Fill(fChan, chAvgMaxValue-chBase );

	return;
}

void Analyze::getPulses(){
	unsigned int numEntries = tr_fit->GetEntries();
	if( numEntries == 0 ) return;
	tr_fit->GetEntry(0);

	fitResponse->clearData();

	int currEvent = -1;
	int eventNum = -1;
	int firstNum = 0;
	for(unsigned int entry = 0 ; entry<numEntries; entry++) { 
		tr_fit->GetEntry(entry);
		if( fEvent != currEvent ){
			currEvent = fEvent;
			eventNum++;
			firstNum = fNum;
		}
		if( eventNum >= constant_maxNumberEventsProcessed )
			break;

		//do pulse selection here
		bool isGood = getPulseStatus();
		//require event has accurate pulser time estimate
		bool isGoodEvent = 1;
		if( evPulserStartSamplesErrors.at(eventNum) > 0.1 )
			isGoodEvent = 0;
		//double startSample = evPulserStartSamples.at(eventNum);
		double startSample = fStartSample;
		fitResponse->addData(fEvent, isGoodEvent, evPulserStartSamples.at(eventNum)*SAMP_PERIOD, fNum, fStartSample*SAMP_PERIOD, fSample, *fWf, isGood);

		if( 0 ){
			std::cout << startSample << "\t" << avgPulserPeriod << std::endl;
			tr_fit->Show();
			drawPulse();
		}
	}
}

void Analyze::doChannelFit(unsigned int chan){
	fitResponse->setSampleError( chRms );
	fitResponse->setBaseFitRange( constant_baseFitRange );
	fitResponse->setPulseFitRange( constant_pulseFitRange );
	fitResponse->showOutput = constant_showFitInfo;
	fitResponse->fixStartTimes = constant_fixStartTimes;

	//fix specific variables
	fitResponse->fixFitVars.push_back(2); //fix baseline
	fitResponse->fixFitVars.push_back(3); //fix period

	//fitResponse->doFit( (chAvgMaxValue - chAvgMean )/constant_ampToHeightFactor, chAvgRiseTime+0.4, chAvgMean, 106.871);
	//fitResponse->doFit( (chAvgMaxValue - chBase )/constant_ampToHeightFactor, 2.2, chBase, chPulserPeriod);
	fitResponse->doFit( (chAvgMaxValue - chBase )/constant_ampToHeightFactor, 2.2, chBase, avgPulserPeriod*SAMP_PERIOD);

	if( fitResponse->fitVals.size() < 6 || fitResponse->fitValErrs.size() < 6 )
		return;

	hFitPulseHeightVsChan->SetBinContent(chan + 1,fitResponse->fitVals.at(0)*constant_ampToHeightFactor );
	hFitPulseHeightVsChan->SetBinError(chan + 1,fitResponse->fitValErrs.at(0)*constant_ampToHeightFactor );

	hFitAmpVsChan->SetBinContent(chan + 1,fitResponse->fitVals.at(0) );
	hFitAmpVsChan->SetBinError(chan + 1,fitResponse->fitValErrs.at(0) );

	hFitShapeVsChan->SetBinContent(chan + 1,fitResponse->fitVals.at(1) );
	hFitShapeVsChan->SetBinError(chan + 1,fitResponse->fitValErrs.at(1) );

	hFitStartVsChan->SetBinContent(chan + 1,fitResponse->fitVals.at(5) );
	hFitStartVsChan->SetBinError(chan + 1,fitResponse->fitValErrs.at(5) );

	if(constant_printResults){
		std::cout << std::endl;
		std::cout << "Channel " << fChan << std::endl;
		std::cout << std::setprecision(6) << std::endl;
		std::cout << "\tAmp\t" << fitResponse->fitVals.at(0) << "\t+/-\t" << fitResponse->fitValErrs.at(0) << std::endl;
		std::cout << "\tShape\t" << fitResponse->fitVals.at(1) << "\t+/-\t" << fitResponse->fitValErrs.at(1) << std::endl;
		std::cout << "\tBase\t" << fitResponse->fitVals.at(2) << "\t+/-\t" << fitResponse->fitValErrs.at(2) << std::endl;
		std::cout << "\tPeriod\t" << fitResponse->fitVals.at(3) << "\t+/-\t" << fitResponse->fitValErrs.at(3) << std::endl;
		std::cout << "\tOffset\t" << fitResponse->fitVals.at(4) << "\t+/-\t" << fitResponse->fitValErrs.at(4) << std::endl;
		std::cout << "\tStart 0\t" << fitResponse->fitVals.at(5) << "\t+/-\t" << fitResponse->fitValErrs.at(5) << std::endl;

		//drawFit();
		//char ct;
		//std::cin >> ct;
	}
}

void Analyze::drawPulse(){
	//load pulse into graph object, do NOT include stuck codes, convert sample number to time (us)
	//tr_fit->Show();
	gCh->Set(0);
	for( int s = 0 ; s < fWf->size() ; s++ )
		gCh->SetPoint(gCh->GetN() , fSample + s , fWf->at(s) );
		//gCh->SetPoint(gCh->GetN() , fSample + s , fWf->at(s) - chBase );
	c0->Clear();
	std::string title = "Channel " + to_string( fChan ) + ", Event " + to_string(fEvent);
	gCh->SetTitle( title.c_str() );
	gCh->GetXaxis()->SetTitle("Sample Number");
	gCh->GetYaxis()->SetTitle("Sample Value (ADC counts)");
	//gCh->GetXaxis()->SetRangeUser(0,128);
	//gCh->GetYaxis()->SetRangeUser(500,1000);
	gCh->SetMarkerStyle(21);
	gCh->Draw("AP");
	c0->Update();
	char ct;
	std::cin >> ct;
}

void Analyze::drawFit(){

	unsigned int numEntries = tr_fit->GetEntries();
	if( numEntries == 0 ) return;
	tr_fit->GetEntry(0);

	int numEvents = fitResponse->fitVals.size()-5;
	if( numEvents <= 0 )
		return;

	double amp = fitResponse->fitVals.at(0);
	double shape = fitResponse->fitVals.at(1);
	double base = fitResponse->fitVals.at(2);
	double period = fitResponse->fitVals.at(3);
	double offset = fitResponse->fitVals.at(4);
	
	for(unsigned int ev = 0 ; ev < fitResponse->fitEventData->size() ; ev++ ){
		if( fitResponse->fitEventData->at(ev).pulseData.size() == 0 || fitResponse->fitEventData->at(ev).isGoodEvent == 0 )
			continue;

  		for(unsigned int p = 0 ; p < fitResponse->fitEventData->at(ev).pulseData.size() ; p ++ ){
			gCh->Set(0);
			gFit->Set(0);

			int pulseEvent = fitResponse->fitEventData->at(ev).eventNumber;
			int pulseFirstSample = fitResponse->fitEventData->at(ev).pulseData.at(p).firstSample;
			int pulseNum = fitResponse->fitEventData->at(ev).pulseData.at(p).num;
			std::vector<unsigned short> *pulseWf = &fitResponse->fitEventData->at(ev).pulseData.at(p).wf;

			for( int s = 0 ; s < pulseWf->size() ; s++ )
				gCh->SetPoint( gCh->GetN(), pulseFirstSample + s , pulseWf->at(s) );

			double startTime = fitResponse->fitVals.at(5+ev) + pulseNum*period + offset;
			for( int s = 0 ; s < pulseWf->size() ; s++ ){
				for(int step = 0 ; step < 10 ; step++){
					double sample = pulseFirstSample + s + step/10.;
					double time = sample*SAMP_PERIOD;
					double fitVal = -10000;
					fitSig->getSignalValue(time, base, startTime,shape, amp, fitVal);
					gFit->SetPoint( gFit->GetN(), sample, fitVal );
				}
			}

			if( 1 ){
				c0->Clear();

				std::string title = "Event Fit Result: Channel " + to_string( fChan ) + ", Event " + to_string(pulseEvent);
				title += ", Pulse # " + to_string(pulseNum);
				gCh->SetTitle( title.c_str() );
				gCh->GetXaxis()->SetTitle("Sample Number");
				gCh->GetYaxis()->SetTitle("Sample Value (ADC counts)");
	
				gCh->SetMarkerStyle(21);
				gCh->Draw("AP");
		
				gFit->SetLineColor(kRed);
				gFit->Draw("L");
		
				c0->Update();

				//sleep(0.1);
				char ct;
				std::cin >> ct;
			}
		}
	}

	return;
}

void Analyze::getAvgFitResidual(){

	unsigned int numEntries = tr_fit->GetEntries();
	if( numEntries == 0 ) return;
	tr_fit->GetEntry(0);

	if( fitResponse->fitVals.size() < 6 || fitResponse->fitValErrs.size() < 6 )
		return;

	int numEvents = fitResponse->fitVals.size()-5;
	if( numEvents <= 0 )
		return;

	double amp = fitResponse->fitVals.at(0);
	double shape = fitResponse->fitVals.at(1);
	double base = fitResponse->fitVals.at(2);
	double period = fitResponse->fitVals.at(3);
	double offset = fitResponse->fitVals.at(4);

	gCh->Set(0);
	gFit->Set(0);

	for(unsigned int ev = 0 ; ev < fitResponse->fitEventData->size() ; ev++ ){
		if( fitResponse->fitEventData->at(ev).pulseData.size() == 0 || fitResponse->fitEventData->at(ev).isGoodEvent == 0 )
			continue;

  		for(unsigned int p = 0 ; p < fitResponse->fitEventData->at(ev).pulseData.size() ; p ++ ){
			if( fitResponse->fitEventData->at(ev).pulseData.at(p).isGood == 0 )
				continue;

			int pulseEvent = fitResponse->fitEventData->at(ev).eventNumber;
			int pulseFirstSample = fitResponse->fitEventData->at(ev).pulseData.at(p).firstSample;
			int pulseNum = fitResponse->fitEventData->at(ev).pulseData.at(p).num;
			std::vector<unsigned short> *pulseWf = &fitResponse->fitEventData->at(ev).pulseData.at(p).wf;

			double startTime = fitResponse->fitVals.at(5+ev) + pulseNum*period + offset;
			double startSample = startTime/SAMP_PERIOD;

			for( int s = 0 ; s < pulseWf->size() ; s++ ){
				double sample = pulseFirstSample + s;
				double time = sample*SAMP_PERIOD;
				double fitVal = -10000;
				fitSig->getSignalValue(time, base, startTime, shape, amp, fitVal);
				//pAvgFitSignalVsChan->Fill(fChan, sample - startSample , pulseWf->at(s) - base );
				pAvgFitResidualVsChan->Fill(fChan, sample - startSample , pulseWf->at(s) - fitVal );
				//gCh->SetPoint( gCh->GetN(), sample - startSample , pulseWf->at(s) - fitVal );
				//gCh->SetPoint( gCh->GetN(), sample - startSample , pulseWf->at(s) );
			}
		}
	}

	if( 0 ){
		c0->Clear();

		std::string title = "Event Fit Residual: Channel " + to_string( fChan );
		gCh->SetTitle( title.c_str() );
		gCh->GetXaxis()->SetTitle("Sample Number");
		gCh->GetYaxis()->SetTitle("Sample Value (ADC counts)");

		gCh->SetMarkerStyle(21);
		gCh->Draw("AP");

		//gFit->SetLineColor(kRed);
		//gFit->Draw("L");
	
		c0->Update();
		//sleep(0.1);
		char ct;
		std::cin >> ct;
	}

	return;
}

void Analyze::fillOutputTree(){
	unsigned int numEntries = tr_fit->GetEntries();
	if( numEntries == 0 ) return;
	tr_fit->GetEntry(0);

	if( fitResponse->fitVals.size() < 6 || fitResponse->fitValErrs.size() < 6 )
		return;

	int numEvents = fitResponse->fitVals.size()-5;
	if( numEvents <= 0 )
		return;

	oRun = fRun;
	oSubrun = fSubrun;
	oChan = fChan;
	oStatus = fitResponse->status;
	oBase = fitResponse->fitVals.at(2);
	oShape = fitResponse->fitVals.at(1);
	oAmp = fitResponse->fitVals.at(0);
	oBaseErr = fitResponse->fitValErrs.at(2);
	oShapeErr = fitResponse->fitValErrs.at(1);
	oAmpErr = fitResponse->fitValErrs.at(0);

	oStartTimes.clear();
	for( unsigned int ev = 0 ; ev < numEvents ; ev++ ){
		oStartTimes.push_back( fitResponse->fitVals.at( 5+ ev ) );
	}

	tr_out->Fill();

	return;
}

void Analyze::getQRes(){

	unsigned int numEntries = tr_fit->GetEntries();
	if( numEntries == 0 ) return;
	tr_fit->GetEntry(0);

	if( fitResponse->fitVals.size() < 6 || fitResponse->fitValErrs.size() < 6 )
		return;

	int numEvents = fitResponse->fitVals.size()-5;
	if( numEvents <= 0 )
		return;

	int currEvent = -1;
	int eventNum = -1;
	int firstNum = 0;
	for(unsigned int entry = 0 ; entry<numEntries; entry++) { 
		tr_fit->GetEntry(entry);
		if( fEvent != currEvent ){
			currEvent = fEvent;
			eventNum++;
			firstNum = fNum;
		}
		if(eventNum >= numEvents)	
			break;

		if( !getPulseStatus() )
			continue;

		if( evPulserStartSamplesErrors.at(eventNum) > 0.1 ) 
			continue;

		doSinglePulseFit(eventNum, firstNum);
	}

	return;
}

void Analyze::doSinglePulseFit(unsigned int eventNum, unsigned int firstNum){

	if( fitResponse->fitVals.size() < 6 || fitResponse->fitValErrs.size() < 6 )
		return;

	int numEvents = fitResponse->fitVals.size()-5;
	if( eventNum >= numEvents )
		return;

	double amp = fitResponse->fitVals.at(0);
	double shape = fitResponse->fitVals.at(1);
	double base = fitResponse->fitVals.at(2);
	double period = fitResponse->fitVals.at(3);
	double offset = fitResponse->fitVals.at(4);

	fitResponse_singlePulse->clearData();	

	double startTime = fitResponse->fitVals.at(5+eventNum) + fNum*period + offset;
	//fitResponse_singlePulse->addData(fEvent, 1, fNum, startTime, fSample, *fWf, 1);

	//fix everything except pulse amplitude
	fitResponse_singlePulse->fixFitVars.push_back(1);
	fitResponse_singlePulse->fixFitVars.push_back(2);
	fitResponse_singlePulse->fixFitVars.push_back(3);
	fitResponse_singlePulse->fixFitVars.push_back(4);
	fitResponse_singlePulse->fixFitVars.push_back(5);

	fitResponse_singlePulse->setSampleError( chAvgRms );
	fitResponse_singlePulse->showOutput = 1;
	fitResponse_singlePulse->setBaseFitRange( constant_baseFitRange );
	fitResponse_singlePulse->setPulseFitRange( constant_pulseFitRange );
	fitResponse_singlePulse->doFit( amp, shape, base, period);

	pSinglePulseAmpVsChan->Fill(fChan, fitResponse_singlePulse->fitVals.at(0) );

	drawFit_singlePulse();

	return;
}

void Analyze::drawFit_singlePulse(){
	double amp_singlePulse = fitResponse_singlePulse->fitVals.at(0);
	double shape_singlePulse = fitResponse_singlePulse->fitVals.at(1);
	double base_singlePulse = fitResponse_singlePulse->fitVals.at(2);
	double period_singlePulse = fitResponse_singlePulse->fitVals.at(3);
	double offset_singlePulse = fitResponse_singlePulse->fitVals.at(4);
	double startTime_singlePulse = fitResponse_singlePulse->fitVals.at(5);

	gCh->Set(0);
	gFit->Set(0);
	for(int s = 0 ; s < fWf->size() ; s++ ){
		for(int step = 0 ; step < 10 ; step++){
			double sample = fSample + s + step/10.;
			double time = sample*SAMP_PERIOD;
			double fitVal = -10000;
			fitSig->getSignalValue(time, base_singlePulse, startTime_singlePulse,shape_singlePulse, amp_singlePulse, fitVal);
			gFit->SetPoint( gFit->GetN(), sample, fitVal );
		}
	}

	for( int s = 0 ; s < fWf->size() ; s++ )
		gCh->SetPoint( gCh->GetN(), fSample + s , fWf->at(s) );

	if( 1 ){
		c0->Clear();
		//tr_fit->Show();

		std::string title = "Single Fit Result: Channel " + to_string( fChan ) + ", Event " + to_string(fEvent);
		gCh->SetTitle( title.c_str() );
		gCh->GetXaxis()->SetTitle("Sample Number");
		gCh->GetYaxis()->SetTitle("Sample Value (ADC counts)");
		gCh->SetMarkerStyle(21);
		gCh->Draw("AP");
	
		gFit->SetLineColor(kRed);
		gFit->Draw("L");
		
		c0->Update();
		//sleep(0.1);
		char ct;
		std::cin >> ct;
	}
	return;
}

void processGetPulserSignals(std::string inputFileName) {

  Analyze ana(inputFileName);
  ana.doAnalysis();

  return;
}

int main(int argc, char *argv[]){
  if(argc!=2){
    cout<<"Usage: processGetPulserSignals [inputFilename]"<<endl;
    return 0;
  }

  std::string inputFileName = argv[1];
  std::cout << "inputFileName " << inputFileName << std::endl;

  //define ROOT application object
  theApp = new TApplication("App", &argc, argv);
  processGetPulserSignals(inputFileName); 

  //return 1;
  gSystem->Exit(0);
}
