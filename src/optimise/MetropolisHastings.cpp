#include <MetropolisHastings.h>
#include <Rand.h>
#include <math.h> //sqrt
#include <TestStatistic.h>
#include <iostream>

const std::vector<double>&
MetropolisHastings::GetSigmas() const{
    return fSigmas;
}

void
MetropolisHastings::SetSigmas(const std::vector<double>& sigmas_){
    fSigmas = sigmas_;
}

unsigned
MetropolisHastings::GetBurnIn() const{
    return fBurnIn;
}

void
MetropolisHastings::SetBurnIn(unsigned burnIn_){
    fBurnIn = burnIn_;
}

unsigned
MetropolisHastings::GetThinFactor() const{
    return fThinFactor;
}

void
MetropolisHastings::SetThinFactor(unsigned thinFac_){
    fThinFactor = thinFac_;
}

unsigned
MetropolisHastings::GetMaxIter() const{
    return fMaxIter;
}

void
MetropolisHastings::SetMaxIter(unsigned maxIter_){
    fMaxIter = maxIter_;
}

double
MetropolisHastings::GetRejectionRate() const{
    return fRejectionRate;
}

std::vector<double>
MetropolisHastings::GetMinima() const{
    return fMinima;
}

void
MetropolisHastings::SetMinima(const std::vector<double>& minima_){
    fMinima = minima_;
}


std::vector<double>
MetropolisHastings::GetMaxima() const{
    return fMaxima;
}

void
MetropolisHastings::SetMaxima(const std::vector<double>& maxima_){
    fMaxima = maxima_;
}


void
MetropolisHastings::SetFlipSign(bool b_){
    fFlipSign = b_;
}

bool
MetropolisHastings::GetFlipSign() const{
    return fFlipSign;
}

bool
MetropolisHastings::GetTestStatLogged() const{
    return fTestStatLogged;
}

void
MetropolisHastings::SetTestStatLogged(bool b_){
    fTestStatLogged = b_;
}

void
MetropolisHastings::SetInitialTrial(const std::vector<double>& trial_){
    fInitialTrial = trial_;
}

std::vector<double>
MetropolisHastings::GetInitialTrial() const{
    return fInitialTrial;
}

const FitResult&
MetropolisHastings::Optimise(TestStatistic* testStat_){
    testStat_->RegisterFitComponents();
    pTestStatistic = testStat_;
    fNDims = fMinima.size();
    fRejectionRate = 0;
    fBestFit.resize(fNDims, 0);
    fMaxVal = 0;
    
    // Set up the histogram 
    AxisCollection axes;
    std::vector<std::string> paramNames = pTestStatistic->GetParameterNames();
    for(size_t i = 0; i < fMinima.size(); i++){
        axes.AddAxis(PdfAxis(paramNames.at(i), fMinima.at(i), fMaxima.at(i), 
                             int(pow(fMaxIter, 1./fMinima.size()))));
    }
    fHist = Histogram(axes);

    // 1. Choose a random starting point or the user defined one
    std::vector<double> currentStep;
    if(fInitialTrial.size())
      currentStep = fInitialTrial;
    else
      for(size_t i = 0; i < fMinima.size(); i++){
          currentStep.push_back(fMinima.at(i) + Rand::Uniform() * (fMaxima.at(i) - fMinima.at(i)));
      }

    std::cout << "Metropolis Hastings::Initial Position @ \t\t";
    for(size_t i = 0; i < currentStep.size(); i++)
        std::cout << currentStep.at(i) << ",";
    std::cout << std::endl;

    fSample.empty();
    fSample.reserve(fMaxIter - fBurnIn);
   
    // 2. Loop step through the space a fixed number of times and
    for(unsigned i = 0; i < fMaxIter; i++){
        // a. Save the point in question if you are past burn-in phase and according to thinning
        if (i > fBurnIn && !(i%fThinFactor)){
            fSample.push_back(currentStep);
            fHist.Fill(currentStep);
        }
      
        if(!(i%1000))
            std::cout << i << "  /  " << fMaxIter 
		      << "\t" << (i - fRejectionRate) / static_cast<double>(i) << std::endl;

        // b. Propose a new step according to a random jump distribution
        std::vector<double> proposedStep = JumpDraw(currentStep);

        // c. Test for acceptance
        if (StepAccepted(currentStep, proposedStep))
            currentStep = proposedStep;
        else
            fRejectionRate++;
    }
    fRejectionRate/= static_cast<double>(fMaxIter);
    std::cout << "Metropolis Hastings:: acceptance rate = " << 1 - fRejectionRate << std::endl;

    fFitResult.SetBestFit(fBestFit);
    fFitResult.SetParameterNames(pTestStatistic->GetParameterNames());
    fFitResult.SetStatSample(fSample);
    fFitResult.SetStatSpace(fHist);
    fFitResult.SetValid(true);
    return fFitResult;
}

bool
MetropolisHastings::StepAccepted(const std::vector<double>& thisStep_, 
                                 const std::vector<double>& proposedStep_){

    // dont step outside of the fit region 
    for (size_t i = 0; i < fMinima.size(); i++){
        if (proposedStep_.at(i) < fMinima.at(i) || proposedStep_.at(i) > fMaxima.at(i))
            return false;
    }

    pTestStatistic -> SetParameters(thisStep_);
    double thisVal = pTestStatistic -> Evaluate();
    
    pTestStatistic -> SetParameters(proposedStep_);
    double proposedVal = pTestStatistic -> Evaluate();

    if(fFlipSign){
        thisVal = -thisVal;
        proposedVal = -proposedVal;
    }

    if(thisVal > fMaxVal || fMaxVal == 0.0){
        fMaxVal = thisVal;
        fBestFit = thisStep_;
    }
    
    double acceptanceParam = 0;
    if(fTestStatLogged)
        acceptanceParam = JumpProbRatio(thisStep_, proposedStep_) * exp(proposedVal - thisVal);

    else
        acceptanceParam = JumpProbRatio(thisStep_, proposedStep_) * proposedVal/thisVal;
    
    if (acceptanceParam > 1)
        return true;

    else if (acceptanceParam >= Rand::Uniform())
        return true;

    else{
        return false;
    }
}

std::vector<double> 
MetropolisHastings::JumpDraw(const std::vector<double>& thisStep_) const{
    std::vector<double> newStep(thisStep_);
    for(size_t i = 0; i < newStep.size(); i++){
        newStep[i] = thisStep_.at(i) + Rand::Gaus(0, fSigmas.at(i));
    }
    return newStep;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
// std::vector<double>                                                                             //
// MetropolisHastings::JumpDraw(const std::vector<double>& thisStep_) const{                       //
//     double stepSize = Rand::Gaus(0, fSigma);                                                    //
//                                                                                                 //
//     // choose random direction on NParam -1 Sphere according to                                 //
//     // http://mathworld.wolfram.com/HyperspherePointPicking.html (Muller 1959, Marsaglia 1972)  //
//                                                                                                 //
//     std::vector<double> newPoint;                                                               //
//     double norm = 0;                                                                            //
//     for(size_t i = 0; i < fNDims; i++){                                                         //
//         double component = Rand::Gaus();                                                        //
//         newPoint.push_back(component);                                                          //
//         norm += component * component;                                                          //
//     }                                                                                           //
//                                                                                                 //
//     norm = sqrt(norm);                                                                          //
//     for(size_t i = 0; i < fNDims; i++){                                                         //
//         newPoint[i] *= stepSize/norm;                                                           //
//     }                                                                                           //
//                                                                                                 //
//     // Add to orginal vector and return                                                         //
//     for(size_t i = 0; i < fNDims; i++)                                                          //
//         newPoint[i] += thisStep_.at(i);                                                         //
//                                                                                                 //
//     return newPoint;                                                                            //
// }                                                                                               //
/////////////////////////////////////////////////////////////////////////////////////////////////////
