#ifndef GETPULSERSHAPE_H
#define GETPULSERSHAPE_H

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

namespace getpulsershape {

  class GetPulserShape : public art::EDAnalyzer {

  public:
    explicit GetPulserShape(fhicl::ParameterSet const& pset);
    virtual ~GetPulserShape();

    void analyze(art::Event const& evt);
    void reconfigure(fhicl::ParameterSet const& pset);
    
    void beginJob();

    void calculateAverageWaveform( std::vector<raw::RawDigit> const& rawDigitVector);
    void findPulses();
    void checkRisingEdge(double threshold, int s, double samp, double sampNext);
    void recordPulserSamples(art::Event const& evt);
    
  private:

    //******************************
    //Variables Taken from FHICL File
    std::string fRawDigitModuleLabel;
    unsigned int fMinSubRun, fMaxSubRun;
    unsigned int fPedRange, fPreRange, fPostRange, fMaxNumPulse;

    //Tree variables
    const unsigned int numChan = 8256;
    const unsigned int numTicks = 9594;
    const int minCode = -5000;
    const int maxCode = 5000;
    const int minThresholdVal = 50;
    TTree *tPulserStartSamples;
    unsigned int fRun, fSubrun, fEvent;
    std::vector<double> fPulserStartSamples;

    //Other variables
    TProfile *pWave;
    TH1I *hSamp;
    std::vector<double> pulseStartSamples;
  }; //end class Noise


  //-------------------------------------------------------------------
  GetPulserShape::GetPulserShape(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset){ 
    this->reconfigure(pset); 
  }

  //-------------------------------------------------------------------
  GetPulserShape::~GetPulserShape(){}

  //-------------------------------------------------------------------
  void GetPulserShape::reconfigure(fhicl::ParameterSet const& pset){
    fMinSubRun 			= pset.get<unsigned int>("minSubRun");
    fMaxSubRun 			= pset.get<unsigned int>("maxSubRun");
    fRawDigitModuleLabel 	= pset.get<std::string>("RawDigitModuleLabel");
    fPedRange             	= pset.get<unsigned int>("pedRange");
    fPreRange             	= pset.get<unsigned int>("preRange");
    fPostRange             	= pset.get<unsigned int>("postRange");
    fMaxNumPulse            	= pset.get<unsigned int>("maxNumPulse");
  }

  //-------------------------------------------------------------------
  void GetPulserShape::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;

    //make output tree for simple pulser time measurements
    tPulserStartSamples = tfs->make<TTree>("tr_pulserStartSamples","tr_pulserStartSamples");
    tPulserStartSamples->Branch("run", &fRun, "run/s");
    tPulserStartSamples->Branch("subrun", &fSubrun, "subrun/s");
    tPulserStartSamples->Branch("event", &fEvent, "event/i");
    tPulserStartSamples->Branch("pulserStartSamples", "vector<double>", &fPulserStartSamples);

    //make generic histograms
    pWave = tfs->make<TProfile>("pWave","",numTicks,-0.5,numTicks-0.5);
    hSamp = tfs->make<TH1I>("hSamp","",maxCode-minCode,-minCode-0.5,maxCode-0.5);
  }

  //-------------------------------------------------------------------
  void GetPulserShape::analyze(art::Event const& evt){
    art::ServiceHandle<art::TFileService> tfs;    
    art::Handle< std::vector<raw::RawDigit> > rawDigitHandle;
    evt.getByLabel(fRawDigitModuleLabel,rawDigitHandle);
    std::vector<raw::RawDigit> const& rawDigitVector(*rawDigitHandle);

    //check that event variables are reasonable
    if( evt.subRun() < fMinSubRun || evt.subRun() > fMaxSubRun)
      return;
    //check that input fhicl parameters are reasonable
    if( fPedRange == 0 || fPreRange == 0 || fPostRange == 0 )
      return;

    //average all channel waveforms for easy pulser signal ID
    calculateAverageWaveform( rawDigitVector );

    //find pulser signal times
    findPulses();

    //record pulse times
    recordPulserSamples(evt);

    return;
  }//end analyze function

  //-------------------------------------------------------------------
  void GetPulserShape::calculateAverageWaveform( std::vector<raw::RawDigit> const& rawDigitVector){
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
  void GetPulserShape::findPulses(){
    //get overall average waveform mean and max/mins samples
    double maxVal = minCode;
    double minVal = maxCode;
    hSamp->Reset();
    for(int s = 0 ; s < pWave->GetNbinsX() ; s++){
      hSamp->Fill( pWave->GetBinContent(s+1) );
      if( pWave->GetBinContent(s+1) > maxVal )
        maxVal = pWave->GetBinContent(s+1);
      if( pWave->GetBinContent(s+1) < minVal )
        minVal = pWave->GetBinContent(s+1);
    }
    double mean = double( hSamp->GetBinCenter( hSamp->GetMaximumBin() ) );
    if( mean <= 0 || mean >= maxCode )
      return;

    //calculate pulser signal threshold
    double threshold = (maxVal - mean)/2.;
    if( threshold < minThresholdVal )
      threshold = minThresholdVal;
    threshold = mean + threshold;

    //identify unbiased pulser times
    pulseStartSamples.clear();
    for( int s = 0 ; s < pWave->GetNbinsX() - 1 ; s++){
      double samp = pWave->GetBinContent(s+1);
      double sampNext = pWave->GetBinContent(s+2);
      //rising edge
      checkRisingEdge(threshold,s,samp,sampNext);
    }
  }//end find pulses

  //-------------------------------------------------------------------
  void GetPulserShape::checkRisingEdge(double threshold, int s, double samp, double sampNext){
    if( !(samp <= threshold && sampNext > threshold) )
      return;

    //linear interpolation between neighboring samples to get the pulser time
    double slope = ( sampNext - samp );
    if( slope <= 0 )
      return;

    //double pulseStart = s - (samp - mean) / slope; //project to start time
    double pulseStart = s + (threshold - samp) / slope; //threshold time
    if( pulseStart > 0 && pulseStart < pWave->GetNbinsX()-1 )
      pulseStartSamples.push_back( pulseStart );
  }//end check rising edge

  //-------------------------------------------------------------------
  void GetPulserShape::recordPulserSamples(art::Event const& evt){
    fRun = evt.run();
    fSubrun = evt.subRun();
    fEvent = evt.event();
    fPulserStartSamples.clear();

    //record pulser times, stop if max number of pulses per event exceeded
    unsigned int numPulse = 0;
    for(unsigned int p = 0 ; p < pulseStartSamples.size() ; p++){
      fPulserStartSamples.push_back( pulseStartSamples.at(p) );
      if( p > fMaxNumPulse )
        break;
    }//end loop over pulse times

    tPulserStartSamples->Fill();
    return;
  }

  DEFINE_ART_MODULE(GetPulserShape)

} //end namespace GetPulserShape

#endif //GETPULSERSHAPE_H