#ifndef __LIB_FITFEELECRESPONSE_MULTIPULSE__
#define __LIB_FITFEELECRESPONSE_MULTIPULSE__

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "FeElecResponse.hxx"
#include "FeElecResponse.cxx"
#include "TMinuit.h"

//Stupid TMinuit global variables/functions
static void fitFuncML(int& npar, double* gout, double& result, double par[], int flg);
void calcLnL(double par[], double& result);
FeElecResponse *sig;
double sampleErr;
double baseFitRange;
double pulseFitRange;

struct Pulses{
    unsigned short num;
    float startTime;
    unsigned short firstSample;
    std::vector<unsigned short> wf;
    std::vector<bool> wfQuality;
    bool isGood;
};

struct Events{
    unsigned int eventNumber;
    bool isGoodEvent;
    double eventTime;
    std::vector<Pulses> pulseData;
};

std::vector<Events> eventData;

class FitFeElecResponse_multiPulse {
private:

public:

  FitFeElecResponse_multiPulse();
  ~FitFeElecResponse_multiPulse();

  void clearData();
  void addData(unsigned int event, bool isGoodEvent, double eventTime, unsigned short num, float startTime, unsigned short firstSample, const std::vector<unsigned short>& wf, const std::vector<bool>& wfQuality, bool isGood);
  void doFit(double initAmp, double initShape, double initBase, double initPeriod, double initOffset);
  void setSampleError(double err);
  void setBaseFitRange(double range);
  void setPulseFitRange(double range);

  unsigned int numPulses;
  bool showOutput;
  bool fixStartTimes;
  double minLnL;
  int status;
  std::vector<Events> *fitEventData;
  std::vector<double> fitVals;
  std::vector<double> fitValErrs;
  std::vector<unsigned int> fixFitVars;
};

#endif

