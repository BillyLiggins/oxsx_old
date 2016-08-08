#ifndef __UTIL__
#define __UTIL__

#include <TH1D.h>
#include <BinnedPdf.h>        
#include <CutCollection.h>        

#include <string>
#include <iostream>

namespace	UTIL{

				std::vector<std::string> glob( const std::string& path, const std::string& start );

				void makeCutCollection(CutCollection& cutCol);
				TH1D* diffHist(TH1D * h1,TH1D * h2);

				double find2n2bRate(double percentage=0.5,double fidRad=6000);
				//After chatting to Jeanne she said you may need to apply the decaying forumla.



				double scaleForTime(double yearRate,double runtime);


				std::vector<double> normRates(std::vector<double>& rates , double normConstant);

				void normChecker(std::vector<double>& expectedrate, std::vector<double>& fit_result);

}




#endif
