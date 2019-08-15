#include "FitFeElecResponse.hxx"

void fitExample(){
	//declare fit object, ROOT objects
	FeElecResponse *elecRespSig = new FeElecResponse();
	FitFeElecResponse_multiPulse *elecRespFit = new FitFeElecResponse_multiPulse();
	TGraph* gData = new TGraph();
	TGraph* gSim = new TGraph();
	TGraph* gFit = new TGraph();
	TCanvas* c0 = new TCanvas("c0", "c0",1400,800);
	TRandom3* rand = new TRandom3(0);

	//define simulated wf + pulse properties
	int simNum = 50; //number of samples
	double simPeriod = 0.5; //sample time base ie sampling period in AU
	double simBase = 100.; //waveform baseline
	double simStart = 12.; //pulse start TIME, not start sample number
	double simShape = 2.2; //shaping TIME, not expressed in # of samples
	double simAmp = 25000./191.6/0.0988165; //amp factor is not pulse height
	double simNoise = 5.; //Guassian distribution width for random noise added to waveform

	//vectors to hold samples, sample valid flag
	std::vector<unsigned short> wf;
	std::vector<bool> wfQuality;

	//simulate some waveform data, add to sample vectors
	for(int i = 0 ; i < simNum ; i++ ){
		double sample = i;
		double time = i*simPeriod;

		//simulate noise to add to sample
		double noise = rand->Gaus(0,simNoise);

		//Get the simulated signal value
		double sig = 0;
		elecRespSig->getSignalValue(time, simBase, simStart, simShape, simAmp, sig);
		double sim = noise + sig;

		//load simulated value into wf vector
		if(sim < 0 ) sim = 0; //require positive ADC values
		unsigned short simVal = noise + sig; //convert doubles to unsigned short ADC code
		gData->SetPoint(gData->GetN(),time,simVal);
		wf.push_back(simVal);
		wfQuality.push_back(1);
	}

	//FIT simulated data section

	//get initial values for fitter
	double initAmp = simAmp;
	double initShape = simShape/simPeriod;
	double initBase = simBase;
	double initTime = simStart/simPeriod - 2;

	//add waveform to fit object, fit	
	elecRespFit->clearData();
	//addPulseData(evt #, evt sample #, pulse #, pulse start sample #, wf first sample #, wf vector, quality vector);
	elecRespFit->addPulseData(0, 0, 0, 0, 0, wf, wfQuality); 
	elecRespFit->setSampleError(2.5);
	elecRespFit->showOutput = 1;
	//elecRespFit->fixFitVars.push_back(2); //fix parameter 2 ie baseline
	elecRespFit->doFit(initAmp , initShape, initBase, 0, initTime);

	//get fit values from fit object        
  	double fitAmp = elecRespFit->fitVals.at(0);
  	double fitShape = elecRespFit->fitVals.at(1)*simPeriod;
  	double fitBase = elecRespFit->fitVals.at(2);
  	double fitStart = elecRespFit->fitVals.at(4)*simPeriod;
	//std::cout << fitAmp << "\t" << fitShape*simPeriod << "\t" << fitStart*simPeriod << std::endl;

	//load fit shape into tgraph
	gFit->Set(0);
	for(unsigned int s = 0 ; s < wf.size() ; s++ ){
		for(int step = 0 ; step < 10 ; step++){
			double sample = s + step/10.;
			double time = sample*simPeriod;
			double val = 0;
			elecRespSig->getSignalValue(time, fitBase, fitStart,fitShape, fitAmp, val);
			gFit->SetPoint( gFit->GetN(), time, val );
		}
	}
        

	//load truth shape into tgraph
	gSim->Set(0);
	for(unsigned int s = 0 ; s < wf.size() ; s++ ){
		for(int step = 0 ; step < 10 ; step++){
			double sample = s + step/10.;
			double time = sample*simPeriod;
			double val = 0;
			elecRespSig->getSignalValue(time, simBase, simStart,simShape, simAmp, val);
			gSim->SetPoint( gSim->GetN(), time, val );
		}
	}

	//draw graphs
	gData->SetMarkerStyle(21);
	gData->SetMarkerColor(kBlack);
	gFit->SetLineColor(kBlue);
	gSim->SetLineColor(kRed);
	gData->GetXaxis()->SetTitle("Sample Time");
	gData->GetYaxis()->SetTitle("Sample Value");

        auto legend = new TLegend(0.6,0.7,0.9,0.9);
        legend->AddEntry(gData,"Data","p");
        legend->AddEntry(gFit,"Fit","l");
        legend->AddEntry(gSim,"Truth","l");

	c0->Clear();
	gData->Draw("AP");
	gFit->Draw("LP");
	gSim->Draw("LP");
        legend->Draw("same");
	c0->Update();

}
