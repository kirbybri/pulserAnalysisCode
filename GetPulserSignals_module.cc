#ifndef GETPULSERSIGNALS_H
#define GETPULSERSIGNALS_H

#include <string>
#include <vector>
#include <iostream>
#include <TH1I.h>
#include <TProfile.h>
#include <TTree.h>
#include <TMath.h>

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "lardataobj/RawData/RawDigit.h"

namespace getpulsersignals {

  class GetPulserSignals : public art::EDAnalyzer {

  public:
    explicit GetPulserSignals(fhicl::ParameterSet const& pset);
    virtual ~GetPulserSignals();

    void analyze(art::Event const& evt);
    void reconfigure(fhicl::ParameterSet const& pset);
    
    void beginJob();

    void calculateAverageWaveform( std::vector<raw::RawDigit> const& rawDigitVector);
    void findPulses();
    void processEventPulses( art::Event const& evt, std::vector<raw::RawDigit> const& rawDigitVector);
    void processPulse(art::Event const& evt, unsigned int pulseNum, double pulseStartSample, raw::RawDigit const& rawDigit );
    double calculatePedestalMean( unsigned int sampleNum, raw::RawDigit const& rawDigit );
    double calculatePedestalRms(unsigned int sampleNum, double mean, raw::RawDigit const& rawDigit);
    void calculatePulsePeak( unsigned int sampleNum, double mean, raw::RawDigit const& rawDigit, short &max, short &min, short &maxSamp, short &minSamp, double &sum);
    void recordPulserSamples(art::Event const& evt);
    
  private:

    //******************************
    //Variables Taken from FHICL File
    std::string fRawDigitModuleLabel;
    unsigned int fMinSubRun, fMaxSubRun;
    unsigned int fPedRange, fPreRange, fPostRange, fMaxNumPulse, fSaveWf;

    //Tree variables
    const unsigned int numChan = 8256;
    TTree* fPulse[8256];
    unsigned short fRun, fSubrun, fChan, fNum, fFirstSample, fMaxValue, fMinValue, fMaxSample, fMinSample;
    unsigned int fEvent;
    float fPulseStartSample, fMean, fRms, fSum;
    std::vector<short> fWf;

    TTree *tPulserStartSamples;
    std::vector<double> fPulserStartSamples;

    //Other variables
    TProfile *pWave;
    TH1I *hSamp;
    std::vector<double> pulseStartSamples;
  }; //end class Noise


  //-------------------------------------------------------------------
  GetPulserSignals::GetPulserSignals(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset){ 
    this->reconfigure(pset); 
  }

  //-------------------------------------------------------------------
  GetPulserSignals::~GetPulserSignals(){}

  //-------------------------------------------------------------------
  void GetPulserSignals::reconfigure(fhicl::ParameterSet const& pset){
    fMinSubRun 			= pset.get<unsigned int>("minSubRun");
    fMaxSubRun 			= pset.get<unsigned int>("maxSubRun");
    fRawDigitModuleLabel 	= pset.get<std::string>("RawDigitModuleLabel");
    fPedRange             	= pset.get<unsigned int>("pedRange");
    fPreRange             	= pset.get<unsigned int>("preRange");
    fPostRange             	= pset.get<unsigned int>("postRange");
    fMaxNumPulse            	= pset.get<unsigned int>("maxNumPulse");
    fSaveWf	 		= pset.get<unsigned int>("savewf");
  }

  //-------------------------------------------------------------------
  void GetPulserSignals::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;

    for(unsigned int ch = 0 ; ch < numChan ; ch++ ){
	std::string title = "tr_" + std::to_string( ch );
    	fPulse[ch] = tfs->make<TTree>(title.c_str(),title.c_str());
    	fPulse[ch]->Branch("run", &fRun, "run/s");
    	fPulse[ch]->Branch("subrun", &fSubrun, "subrun/s");
    	fPulse[ch]->Branch("event", &fEvent, "event/i");
    	fPulse[ch]->Branch("chan", &fChan, "chan/s");
	fPulse[ch]->Branch("num", &fNum, "num/s");
    	fPulse[ch]->Branch("pulsestartsample", &fPulseStartSample, "pulsestartsample/F");
    	fPulse[ch]->Branch("firstsample", &fFirstSample, "firstsample/s");
    	fPulse[ch]->Branch("maxvalue", &fMaxValue, "maxvalue/s");
    	fPulse[ch]->Branch("minvalue", &fMinValue, "minvalue/s");
    	fPulse[ch]->Branch("maxsample", &fMaxSample, "maxsample/s");
    	fPulse[ch]->Branch("minsample", &fMinSample, "minsample/s");
    	fPulse[ch]->Branch("mean", &fMean, "mean/F");
    	fPulse[ch]->Branch("rms", &fRms, "rms/F");
    	fPulse[ch]->Branch("sum", &fSum, "sum/F");
    	fPulse[ch]->Branch("wf", "vector<short>", &fWf);
    }

    tPulserStartSamples = tfs->make<TTree>("tr_pulserStartSamples","tr_pulserStartSamples");
    tPulserStartSamples->Branch("run", &fRun, "run/s");
    tPulserStartSamples->Branch("subrun", &fSubrun, "subrun/s");
    tPulserStartSamples->Branch("event", &fEvent, "event/i");
    tPulserStartSamples->Branch("pulserStartSamples", "vector<double>", &fPulserStartSamples);

    pWave = tfs->make<TProfile>("pWave","",9594,-0.5,9594-0.5);
    hSamp = tfs->make<TH1I>("hSamp","",10000,-5000-0.5,5000-0.5);
  }

  //-------------------------------------------------------------------
  void GetPulserSignals::analyze(art::Event const& evt){

    art::ServiceHandle<art::TFileService> tfs;    
    art::Handle< std::vector<raw::RawDigit> > rawDigitHandle;
    evt.getByLabel(fRawDigitModuleLabel,rawDigitHandle);
    std::vector<raw::RawDigit> const& rawDigitVector(*rawDigitHandle);

    //check that event variables are reasonable
    if( evt.subRun() < fMinSubRun || evt.subRun() > fMaxSubRun)
	return;
    if( pWave->GetNbinsX() <= 1 )
	return;
    if( fPedRange == 0 || fPreRange == 0 || fPostRange == 0 )
	return;

    calculateAverageWaveform( rawDigitVector );

    findPulses();

    //ignore noisy events
    if( pulseStartSamples.size() > fMaxNumPulse )
	return;

    //loop over channels, for each global pulse time save waveform section info
    processEventPulses( evt, rawDigitVector );

    //record pulse times
    recordPulserSamples(evt);

    return;
  }//end analyze function

  //-------------------------------------------------------------------
  void GetPulserSignals::calculateAverageWaveform( std::vector<raw::RawDigit> const& rawDigitVector){
    //loop over channels, get overall average waveform for event
    const unsigned int n_channels = rawDigitVector.size();
    if( n_channels == 0 )
	return;
    pWave->Reset();
    for( unsigned int ich = 0 ; ich < n_channels ; ich++ ){
	const size_t n_samp = rawDigitVector.at(ich).NADC();
	if( n_samp == 0 ) 
		continue;
	for( unsigned int s = 0 ; s < n_samp ; s++ )
		pWave->Fill( s , double( rawDigitVector.at(ich).ADC(s) ) );
    }
  }

  //-------------------------------------------------------------------
  void GetPulserSignals::findPulses(){
    //get overall average waveform mean and max sample
    double maxVal = -5000;
    hSamp->Reset();
    for(int s = 0 ; s < pWave->GetNbinsX() ; s++){
	hSamp->Fill( pWave->GetBinContent(s+1) );
	if( pWave->GetBinContent(s+1) > maxVal )
		maxVal = pWave->GetBinContent(s+1);
    }
    double mean = double( hSamp->GetBinCenter( hSamp->GetMaximumBin() ) );
    if( mean <= 0 || mean >= 4095 )
	return;

    //calculate pulser signal threshold
    double fThreshold = (maxVal - mean)/2.;
    if( fThreshold < 50 )
	fThreshold = 50;
    fThreshold = mean + fThreshold;

    //identify unbiased pulse times
    pulseStartSamples.clear();
    for( int s = 0 ; s < pWave->GetNbinsX() - 1 ; s++){
	double samp = pWave->GetBinContent(s+1);
	double sampNext = pWave->GetBinContent(s+2);
	//rising edge
	if( samp <= fThreshold && sampNext > fThreshold ){
		double slope = ( sampNext - samp );
		if( slope > 0 ){
			double pulseStart = s - (samp - mean) / slope;
			if( pulseStart > 0 && pulseStart < pWave->GetNbinsX()-1 ){
				pulseStartSamples.push_back( pulseStart );
			}
		}
	}
    }
  }//end find pulses

  //-------------------------------------------------------------------
  void GetPulserSignals::processEventPulses( art::Event const& evt, std::vector<raw::RawDigit> const& rawDigitVector){
    const unsigned int n_channels = rawDigitVector.size();
    if( n_channels == 0 )
	return;

    for(unsigned int ich=0; ich<n_channels; ich++){
	const size_t n_samp = rawDigitVector.at(ich).NADC();
	if( n_samp == 0 ) 
		continue;

	//loop over waveforms sections
	unsigned int numPulse = 0;
    	for(unsigned int p = 0 ; p < pulseStartSamples.size() ; p++){
		processPulse( evt, p, pulseStartSamples.at(p), rawDigitVector.at(ich) );
		//count number of pulses
		numPulse++;
		if( numPulse > fMaxNumPulse )
			break;
	}//end loop over pulse times
    }//end loop over channels
  }

  //-------------------------------------------------------------------
  void GetPulserSignals::processPulse(art::Event const& evt, unsigned int pulseNum, double pulseStartSample, raw::RawDigit const& rawDigit ){
	const size_t n_samp = rawDigit.NADC();
	unsigned int sampleNum = (unsigned int) TMath::FloorNint( pulseStartSample );
	if ( sampleNum < fPedRange + fPreRange + 10 || sampleNum >= n_samp - fPostRange - 10 )
		return;

	//get mean estimate from waveform sections before pulses
	double mean = calculatePedestalMean(sampleNum, rawDigit);
	if( mean < 0 || mean > 4095)
		return;

	//get rms estimate
	double rms = calculatePedestalRms(sampleNum, mean, rawDigit);
	if( rms <= 0 )
		return;

	//find pulse peak
	short max = 0;
	short min = 5000;
	short maxSamp = -1;
	short minSamp = -1;
	double sum = 0;
        calculatePulsePeak(sampleNum, mean, rawDigit, max, min, maxSamp, minSamp, sum);
	if( max <= 0 || max >= 4095 )
		return;
	if( min < 0 || min > 4095 )
		return;
	if( maxSamp < 0 || minSamp < 0 )
		return;

	fRun = evt.run();
    	fSubrun = evt.subRun();
    	fEvent = evt.event();
	fChan = rawDigit.Channel();
	fNum = pulseNum;
	fPulseStartSample =  pulseStartSample;
	fFirstSample = sampleNum - fPreRange; //index of first sample in save waveform
	fMaxValue = max;
	fMinValue = min;
 	fMaxSample = maxSamp;
 	fMinSample = minSamp;
	fMean = mean;
 	fRms = rms;
	fSum = sum;

	fWf.clear();
	if( fSaveWf == 1 ){
		for(unsigned int s = fFirstSample ; s < fFirstSample + fPreRange + fPostRange ; s++){
			if( s >= rawDigit.NADC() )
				break;
			fWf.push_back( rawDigit.ADC(s) );
		}
	}

	//save the pulse info for this channel to tree
	if( fChan < numChan )
		fPulse[fChan]->Fill();
  }//end processPulse

  //-------------------------------------------------------------------
  double GetPulserSignals::calculatePedestalMean( unsigned int sampleNum, raw::RawDigit const& rawDigit ){
	if( sampleNum - fPedRange - fPreRange < 0 || sampleNum - fPreRange >= rawDigit.NADC()  )
		return -1;
	double mean = 0;
	int count = 0;
	for(unsigned int s = sampleNum - fPedRange - fPreRange ; s < sampleNum - fPreRange ; s++){
		if( s >= rawDigit.NADC() )
			break;
		mean += rawDigit.ADC(s);
		count++;
	}
	if( count <= 0 )
		return -1;
	mean = mean / (double) count;
	return mean;
  }

  //-------------------------------------------------------------------
  double GetPulserSignals::calculatePedestalRms(unsigned int sampleNum, double mean, raw::RawDigit const& rawDigit){
	if( sampleNum - fPedRange - fPreRange < 0 || sampleNum - fPreRange >= rawDigit.NADC()  )
		return -1;
	double rms = 0;
	int count = 0;
	for(unsigned int s = sampleNum - fPedRange - fPreRange ; s < sampleNum - fPreRange ; s++){
		if( s >= rawDigit.NADC() )
			break;
		rms += ( rawDigit.ADC(s) - mean )*( rawDigit.ADC(s) - mean );
		count++;
	}
	if( count <= 1 )
		return -1;
	rms = TMath::Sqrt( rms / (double)(count - 1 ) );
	return rms;
  }

  //-------------------------------------------------------------------
  void GetPulserSignals::calculatePulsePeak( unsigned int sampleNum, double mean, raw::RawDigit const& rawDigit, short &max, short &min, short &maxSamp, short &minSamp, double &sum){
	max = 0;
	min = 5000;
	maxSamp = -1;
	minSamp = -1;
	sum = 0;
	if( sampleNum - fPreRange < 0 || sampleNum + fPostRange >= rawDigit.NADC()  )
		return;
	for(unsigned int s = sampleNum - fPreRange ; s < sampleNum + fPostRange ; s++){
		if( s >= rawDigit.NADC() )
			break;
		if( rawDigit.ADC(s) > max ){
			max = rawDigit.ADC(s);
			maxSamp = s;
		}
		if( rawDigit.ADC(s) < min){
			min = rawDigit.ADC(s);
			minSamp = s;
		}
		sum += rawDigit.ADC(s) - mean;
	}
	return;
  }

  //-------------------------------------------------------------------
  void GetPulserSignals::recordPulserSamples(art::Event const& evt){

	fRun = evt.run();
    	fSubrun = evt.subRun();
    	fEvent = evt.event();

	fPulserStartSamples.clear();
	//loop over waveforms sections
	unsigned int numPulse = 0;
    	for(unsigned int p = 0 ; p < pulseStartSamples.size() ; p++){
		fPulserStartSamples.push_back( pulseStartSamples.at(p) );
		//count number of pulses
		numPulse++;
		if( numPulse > fMaxNumPulse )
			break;
	}//end loop over pulse times

	tPulserStartSamples->Fill();
  }

  DEFINE_ART_MODULE(GetPulserSignals)

} //end namespace GetPulserSignals

#endif //GETPULSERSIGNALS_H
