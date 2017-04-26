#include "FeElecResponse.hxx"

using namespace std;

FitFeElecResponse_multiPulse::FitFeElecResponse_multiPulse(){
  	sig = new FeElecResponse();
	showOutput = 0;
	sampleErr = 0.1;
   	minLnL = 0;
	fixStartTimes = 1;
	status = -1;
	numPulses = 0;
	baseFitRange = 100;//default, us
	pulseFitRange = 6; //default, us
	samplePeriod = 0.5; //default, us
	fitEventData = &eventData;
}

FitFeElecResponse_multiPulse::~FitFeElecResponse_multiPulse(){
	delete sig;
}

void FitFeElecResponse_multiPulse::clearData(){
	eventData.clear();
	numPulses = 0;
	fixFitVars.clear();
}

void FitFeElecResponse_multiPulse::addData(unsigned int event, bool isGoodEvent, double eventTime, unsigned short num, float startTime, unsigned short firstSample, 
const std::vector<unsigned short>& wf, const std::vector<bool>& wfQuality, bool isGood){
	if( eventData.size() == 0 ){
		Events evTemp;
		evTemp.eventNumber = event;
		evTemp.isGoodEvent= isGoodEvent;
		evTemp.eventTime = eventTime;
		eventData.push_back( evTemp );
	}
	else if( event != eventData.back().eventNumber ){
		Events evTemp;
		evTemp.eventNumber = event;
		evTemp.isGoodEvent = isGoodEvent;
		evTemp.eventTime = eventTime;
		eventData.push_back( evTemp );
	}
	
 	Pulses pTemp;
	pTemp.num = num;
	pTemp.startTime = startTime;
	pTemp.firstSample = firstSample;
	pTemp.wf = wf;
	pTemp.wfQuality = wfQuality;
	pTemp.isGood = isGood;
	eventData.back().pulseData.push_back(pTemp);
	numPulses++;
	return;
}

//wrapper function for TMinuit
void FitFeElecResponse_multiPulse::doFit(double initAmp, double initShape, double initBase, double initPeriod, double initOffset){

  if( numPulses == 0 )
	return;

  //test to see if there are data point
  status = -1;

  //unsigned int numParameters = 4 + eventData.size();
  unsigned int numParameters = 5 + eventData.size();
  //std::unique_ptr<TMinuit> minimizer (new TMinuit(numParameters) );
  TMinuit *minimizer = new TMinuit(numParameters);

  //Set print level , -1 = suppress, 0 = info
  minimizer->SetPrintLevel(-1);
  if( showOutput == 1 )
  	minimizer->SetPrintLevel(3);

  //define parameters
  minimizer->SetFCN(fitFuncML);
  minimizer->DefineParameter(0, "Amp", initAmp, initAmp/1000.,0,0);
  minimizer->DefineParameter(1, "Shape", initShape, initShape/1000.,0,0);
  minimizer->DefineParameter(2, "Base", initBase, initBase/1000.,0,0);
  minimizer->DefineParameter(3, "Period", initPeriod, 0.01,0,0);
  minimizer->DefineParameter(4, "Offset", initOffset, 0.02,0,0);

  //event start fit
  for(unsigned int ev = 0 ; ev < eventData.size() ; ev++ ){
	char name[200];
	memset(name,0,sizeof(char)*100 );
        sprintf(name,"Start_%i",ev);
	//minimizer->DefineParameter(4+ev, name, eventData.at(ev).pulseData.front().startTime, 0.01,0,0); //10us error, note hardcoded
	//minimizer->DefineParameter(5+ev, name, eventData.at(ev).pulseData.front().startTime, 0.01,0,0); //10us error, note hardcoded
	minimizer->DefineParameter(5+ev, name, eventData.at(ev).eventTime, 0.01,0,0); //10us error, note hardcoded
  }

  //optionally fix parameter 
  for( unsigned int i = 0 ; i <  fixFitVars.size() ; i++ ){
	if( fixFitVars.at(i) < numParameters )
  		minimizer->FixParameter( fixFitVars.at(i) );
  }
  //fix period variable if only one pulse
  if( numPulses == 1 )
	minimizer->FixParameter(3);
  //fix event start variables if event contains no valid pulses
  for(unsigned int ev = 0 ; ev < eventData.size() ; ev++ ){
	int numGoodPulses = 0;
	for(unsigned int p = 0 ; p < eventData.at(ev).pulseData.size() ; p ++ ){
		if( eventData.at(ev).pulseData.at(p).isGood == 1)
			numGoodPulses++;
	}
	if( numGoodPulses == 0 )
		minimizer->FixParameter(5+ev);
  }
  //fix all pulse start times
  if( fixStartTimes == 1){
  	for(unsigned int ev = 0 ; ev < eventData.size() ; ev++ ){
		minimizer->FixParameter(5+ev);
  	}
  }
  else{
	minimizer->FixParameter(4);
  }

  //Set Minuit flags
  Double_t arglist[10];
  arglist[0] = 0.5;
  Int_t ierflg = 0;
  minimizer->mnexcm("SET ERR", arglist ,1,ierflg);  //command, arguments, # arguments, error flag
        
  //MIGRAD minimization
  Double_t tmp[1];
  tmp[0] = 100000;
  Int_t err;
  minimizer->mnexcm("MIG", tmp ,1,err);
  status = err;

  //HESSE
  //Double_t hesstmp[1];
  //hesstmp[0] = 100000;
  //Int_t hesserr;
  //minimizer->mnexcm("HES", hesstmp ,1,hesserr);
  //status = hesserr;
  
  fitVals.clear();
  fitValErrs.clear();
  for(unsigned int i = 0 ; i < numParameters ; i++ ){
	double fitVal, fitValErr;
	minimizer->GetParameter(i, fitVal, fitValErr);
	fitVals.push_back( fitVal );
	fitValErrs.push_back( fitValErr );
  }

  //int npars = minimizer->GetNumPars(); 
  //minimizer->mnemat(fitParEmat,npars);

  //calculate ln likelihood of fit
  //calcLnL(fitPar,minLnL);

  delete minimizer; //memory leak?

  return;
}

void FitFeElecResponse_multiPulse::setSampleError(double err){
  sampleErr = 0.1;
  if(err > 0.1 )
	sampleErr = err;
  return;
}

void FitFeElecResponse_multiPulse::setPulseFitRange(double range){
  pulseFitRange = 6;
  if( range > 0.5 )
	pulseFitRange = range;
}

void FitFeElecResponse_multiPulse::setBaseFitRange(double range){
  baseFitRange = 100;
  if( range > 0.5 )
	baseFitRange = range;
}

//likelihood calc - note not included in class
void calcLnL(double par[], double& result){
  double diffSq = 0;
  double dataX,fitY;
  unsigned short dataY;

  double amp = par[0];
  double shape = par[1];
  double base = par[2];

  double norm = 1./(sampleErr)/(sampleErr);

  //TCanvas *cFit = new TCanvas("cFit", "cFit",1400,800);
  //TGraph *gData = new TGraph();
  //TGraph *gFit = new TGraph();

  //loop over events, pulses, samples
  for(unsigned int ev = 0 ; ev < eventData.size() ; ev++ ){
	if( eventData[ev].isGoodEvent == 0 )
		continue;
  	for(unsigned int p = 0 ; p < eventData[ev].pulseData.size() ; p ++ ){
		//skip bad pulses
		if( eventData[ev].pulseData[p].isGood == 0)
			continue;
		//double startSample = evPulserStartSamples.at(eventNum) + fNum*avgPulserPeriod;
		//double start = par[3]*(eventData[ev].pulseData[p].num - eventData[ev].pulseData[0].num) + par[4+ev];
		//double start = par[3]*(eventData[ev].pulseData[p].num - eventData[ev].pulseData[0].num) + par[5+ev] + par[4];
		double start = par[3]*eventData[ev].pulseData[p].num + par[5+ev] + par[4];
		//double minFitTime =  par[3]*eventData[ev].pulseData[p].num + eventData[ev].pulseData[p].startTime - baseFitRange;
		//double maxFitTime =  par[3]*eventData[ev].pulseData[p].num + eventData[ev].pulseData[p].startTime + pulseFitRange;
		//double minFitTime =  eventData[ev].pulseData[p].startTime - baseFitRange;
		//double maxFitTime =  eventData[ev].pulseData[p].startTime + pulseFitRange;
	
		//std::cout << amp << "\t" << shape << "\t" << base << std::endl;
		//gData->Set(0);
		//gFit->Set(0);
		for(unsigned int s = 0 ; s < eventData[ev].pulseData[p].wf.size() ; s++ ){ //using data vector directly
			dataX = eventData[ev].pulseData[p].firstSample + s;
	
			//only allow fit to include certain part of pulse waveform
			//if( dataX*samplePeriod < minFitTime )
			//	continue;
			//if( dataX*samplePeriod > maxFitTime )
			//	break;
			if( eventData[ev].pulseData[p].wfQuality[s] == 0 )
				continue;

			dataY = eventData[ev].pulseData[p].wf[s];

			//void getSignalValueFromVector(double time, double base, double pulseStart, double shapeTime, double amp,double& simVal);
			sig->getSignalValueFromVector(dataX*samplePeriod, base, start, shape, amp,fitY);

			//double y = (dataY - fitY)*(dataY - fitY)*norm - cumul;
			//double t = diffSq + y;
			//cumul = (t - diffSq) - y;
			//diffSq = t;

			diffSq = diffSq + (dataY - fitY)*(dataY - fitY)*norm; //gauss err assumed, noise is correlated but assume small

			//std::cout << dataX << "\t" << dataY << "\t" << fitY << std::endl;
			//gData->SetPoint(gData->GetN() , dataX*samplePeriod, dataY);
			//gFit->SetPoint(gFit->GetN() , dataX*samplePeriod, fitY);
  		}
		/*
		cFit->Clear();
		gData->SetMarkerStyle(21);
		gData->Draw("AP");
		gFit->SetLineColor(kRed);
		gFit->Draw("LP");
		cFit->Update();
		*/
		//char ct;
		//std::cin >> ct;
		
  	}//end of pulse loop
  }//end of event loop

  //delete cFit;
  //delete gData;
  //delete gFit;

  //calculate value to minimize
  result = -0.5*diffSq;
  return;
}

//Fit wrapper function - used by Minuit - has to be static void, annoying
static void fitFuncML(int& npar, double* gout, double& result, double par[], int flg){
  calcLnL(par, result);
  result = -1.*result;//Minuit is minimizing result ie maximizing LnL
  return;
}
