#include <AnalyticPdf.h>
#include <Exceptions.h>
#include <IntegrableFunction.h>


AnalyticPdf::AnalyticPdf(IntegrableFunction* f_){
    fFunction = dynamic_cast<IntegrableFunction*>(f_->Clone());
    fNorm     = 1;
}

AnalyticPdf::~AnalyticPdf(){
    delete fFunction;
}

AnalyticPdf::AnalyticPdf(const AnalyticPdf& other_){
    fNorm  = other_.fNorm;
    fDataRep = other_.fDataRep;
    fFunction = dynamic_cast<IntegrableFunction*>(other_.fFunction->Clone());
}

Pdf*
AnalyticPdf::Clone() const{
    return static_cast<Pdf*>(new AnalyticPdf(*this));
}

double
AnalyticPdf::operator()(const std::vector<double>& vals_) const{    
    try{
        return fFunction->operator()(vals_)/fNorm;
    }
    catch(const DimensionError& e_){
        throw DimensionError(std::string("AnalyticPdf internal function ::") + e_.what());
                             
                             
    }
}

double
AnalyticPdf::Probability(const EventData& event_) const{
    try{
        return operator()(event_.ToRepresentation(fDataRep));
    }

    catch(const RepresentationError& e_){
        throw RepresentationError("AnalyticPdf::Probability() failed with  "
                                  + std::string(e_.what()) + 
                                  " is the rep set correctly?");
    }
}

double 
AnalyticPdf::Integral() const{
    return fFunction->Integral();
}

void
AnalyticPdf::Normalise(){
    fNorm = Integral();
}

void
AnalyticPdf::SetDataRep(const DataRepresentation& rep_){
    fDataRep = rep_;
}

DataRepresentation
AnalyticPdf::GetDataRep() const {
    return fDataRep;
}

unsigned
AnalyticPdf::GetNDims() const{
    return fFunction->GetNDims();
}

// Fitting this pdf to data means adjusting the underlying function
void
AnalyticPdf::MakeFittable(){
    fFunction -> MakeFittable();
}

std::vector<std::string> 
AnalyticPdf::GetParameterNames() const{
    std::vector<std::string> funcNames = fFunction->GetParameterNames();
    for(size_t i = 0; i < funcNames.size(); i++)
        funcNames[i] = "Analytic PDF: " + funcNames[i];
    return funcNames;
}

std::vector<double>
AnalyticPdf::GetParameters() const{
    return fFunction->GetParameters();
}

size_t 
AnalyticPdf::GetParameterCount() const{
    return fFunction->GetParameterCount();
}

void
AnalyticPdf::SetParameters(const std::vector<double>& params_){
    try{
        fFunction->SetParameters(params_);
    }
    catch(const ParameterCountError& e_){
        throw ParameterCountError(std::string("AnalyticPdf internal function : ") + e_.what());

    }
}
