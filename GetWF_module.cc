#ifndef GETWF_H
#define GETWF_H

#include <string>
#include <vector>
#include <iostream>
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

namespace calibration {

  class GetWF : public art::EDAnalyzer {

  public:
    explicit GetWF(fhicl::ParameterSet const& pset);
    virtual ~GetWF();

    void analyze(art::Event const& evt);
    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void endJob();

    //likely we will need begin/end run and subrun functions
    void beginRun(art::Run const& run);
    void endRun(art::Run const& run);
    void beginSubRun(art::SubRun const& subrun);
    void endSubRun(art::SubRun const& subrun);
    
  private:

    //******************************
    //Variables Taken from FHICL File
    std::string       fRawDigitModuleLabel;   //label for rawdigit module

    //Other variables
    TTree* fTree;
    short fRun, fSubrun, fEvent, fChan;
    std::vector<short> fWf;

  }; //end class Noise


  //-------------------------------------------------------------------
  GetWF::GetWF(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset){ 
    this->reconfigure(pset); 
  }

  //-------------------------------------------------------------------
  GetWF::~GetWF(){}

  //-------------------------------------------------------------------
  void GetWF::reconfigure(fhicl::ParameterSet const& pset){
    fRawDigitModuleLabel = pset.get<std::string>("RawDigitModuleLabel");
  }

  //-------------------------------------------------------------------
  void GetWF::beginJob(){    
  }

  //-------------------------------------------------------------------
  void GetWF::endJob(){
  }

  //-------------------------------------------------------------------
  void GetWF::beginRun(art::Run const& run){

    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("fTree","_wf_tree");
    fTree->Branch("_run", &fRun, "_run/s");
    fTree->Branch("_subrun", &fSubrun, "_subrun/s");
    fTree->Branch("_event", &fEvent, "_event/s");
    fTree->Branch("_chan", &fChan, "_chan/s");
    fTree->Branch("adc_v", "vector<short>", &fWf);
  }

  //-------------------------------------------------------------------
  void GetWF::endRun(art::Run const& run){

    art::ServiceHandle<art::TFileService> tfs;
    return;
  }


  //-------------------------------------------------------------------
  void GetWF::beginSubRun(art::SubRun const& subrun){

    art::ServiceHandle<art::TFileService> tfs;
  }

  //-------------------------------------------------------------------
  void GetWF::endSubRun(art::SubRun const& subrun){

    art::ServiceHandle<art::TFileService> tfs;
  }
  
  //-------------------------------------------------------------------
  void GetWF::analyze(art::Event const& evt){

    art::ServiceHandle<art::TFileService> tfs;    
    art::Handle< std::vector<raw::RawDigit> > rawDigitHandle;
    evt.getByLabel(fRawDigitModuleLabel,rawDigitHandle);
    std::vector<raw::RawDigit> const& rawDigitVector(*rawDigitHandle);

    //loop over channels
    const unsigned int n_channels = rawDigitVector.size();
    for(unsigned int ich=0; ich<n_channels; ich++){
	const size_t n_samp = rawDigitVector.at(ich).NADC();
	if( n_samp == 0 ) 
		continue;

	fRun = evt.run();
    	fSubrun = evt.subRun();
    	fEvent = evt.event();
	fChan = rawDigitVector.at(ich).Channel();
	
	fWf.clear();
	for(unsigned int s = 0 ; s < n_samp ; s++)
		fWf.push_back( rawDigitVector.at(ich).ADC(s) );

	//save the pulse info for this channel to tree
	fTree->Fill();
    }
  }

  DEFINE_ART_MODULE(GetWF)

} //end namespace GetPulserSignals

#endif //GETWF_H
