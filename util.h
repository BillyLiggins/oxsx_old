#ifndef __UTIL__
#define __UTIL__
#include <TH1D.h>
#include <iostream>

// TH1D* diffHist(TH1D * h1,TH1D * h2);
//TCanvas* doublePlot( const vector<TH1D>& histList);
double find2n2bRate(double percentage=0.5,double fidRad=6000);

double scaleForTime(double yearRate,double runtime);

std::vector<double> normRates(std::vector<double>& rates , double normConstant);
TH1D* diffHist(TH1D * h1,TH1D * h2);

void normChecker(std::vector<double>& expectedrate, std::vector<double>& fit_result);

#endif
