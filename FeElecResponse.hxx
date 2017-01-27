#ifndef __LIB_FEELECRESPONSE__
#define __LIB_FEELECRESPONSE__

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>

class FeElecResponse {
private:

public:

  FeElecResponse();
  ~FeElecResponse();
  
  void getSignalValue(double time, double base, double pulseStart, double shapeTime, double amp,double& simVal);
  void getSignalValueFromVector(double time, double base, double pulseStart, double shapeTime, double amp,double& simVal);
  void responseFunc( double time, double To, double Ao, double& val);
  std::vector<double> sig;
  unsigned int sigNum;
  double sigPeriod;
};

#endif

