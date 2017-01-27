#include "FeElecResponse.hxx"

using namespace std;

FeElecResponse::FeElecResponse(){
	//load response function into vector
	sigNum = 4000;
	sigPeriod = 0.5/200.;
	for(int i = 0 ; i < sigNum ; i++){
		double time = i*sigPeriod; //us
		double val = 0;
		getSignalValue(time, 0., 0., 1., 1., val); //"unit" response
		sig.push_back(val);
  	}
}

FeElecResponse::~FeElecResponse(){
}

//function to get response function value from response function
void FeElecResponse::getSignalValue(double time, double base, double pulseStart, double shapeTime, double amp,double& simVal){
  double pulseTime = 0;
  if( time > pulseStart )
  	pulseTime = time - pulseStart;
  simVal  = 0;
  responseFunc( pulseTime, shapeTime, amp, simVal);
  simVal += base;
  return;
}

//get response function value quickly from vector interpolation
void FeElecResponse::getSignalValueFromVector(double time, double base, double pulseStart, double shapeTime, double amp,double& simVal){
  double pulseTime = 0;
  if( time > pulseStart )
  	pulseTime = time - pulseStart;
  simVal  = base;
  if( shapeTime <= 0 )
	return;

  //determine position of time value in signal vector
  double pulseTimeSamp = pulseTime/shapeTime/sigPeriod;
  unsigned int pulseSamp = floor(pulseTimeSamp);
  if( pulseSamp >= sigNum-1 )
	return;
  
  //do linear interpolation
  double sigVal = sig[pulseSamp] + ( sig[pulseSamp+1] - sig[pulseSamp] )*(pulseTimeSamp - pulseSamp);
  //scale by amplitude factor
  sigVal = amp*sigVal;

  simVal += sigVal;
  return;
}

//FE ASIC response
void FeElecResponse::responseFunc( double time, double To, double Ao, double& val){
	val = 4.31054*exp(-2.94809*time/To)*Ao - 2.6202*exp(-2.82833*time/To)*cos(1.19361*time/To)*Ao
	-2.6202*exp(-2.82833*time/To)*cos(1.19361*time/To)*cos(2.38722*time/To)*Ao
	+0.464924*exp(-2.40318*time/To)*cos(2.5928*time/To)*Ao
	+0.464924*exp(-2.40318*time/To)*cos(2.5928*time/To)*cos(5.18561*time/To)*Ao
	+0.762456*exp(-2.82833*time/To)*sin(1.19361*time/To)*Ao
	-0.762456*exp(-2.82833*time/To)*cos(2.38722*time/To)*sin(1.19361*time/To)*Ao
	+0.762456*exp(-2.82833*time/To)*cos(1.19361*time/To)*sin(2.38722*time/To)*Ao
	-2.6202*exp(-2.82833*time/To)*sin(1.19361*time/To)*sin(2.38722*time/To)*Ao
	-0.327684*exp(-2.40318*time/To)*sin(2.5928*time/To)*Ao
	+0.327684*exp(-2.40318*time/To)*cos(5.18561*time/To)*sin(2.5928*time/To)*Ao
	-0.327684*exp(-2.40318*time/To)*cos(2.5928*time/To)*sin(5.18561*time/To)*Ao
	+0.464924*exp(-2.40318*time/To)*sin(2.5928*time/To)*sin(5.18561*time/To)*Ao;
}
