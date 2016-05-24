#include <PdfManager.h>
#include <Pdf.h>
#include <EventData.h>
#include <Exceptions.h>

PdfManager::~PdfManager(){
    for(size_t i = 0; i < fPdfs.size(); i++)
        delete fPdfs[i];
}

void 
PdfManager::AddPdf(Pdf * pdf_){
    if (!fPdfs.size())
        fNDims = pdf_->GetNDims();

    else if(pdf_->GetNDims() != fNDims)
        throw DimensionError("PdfManager::AddPdf", fNDims, pdf_->GetNDims(),
                             " dimensions in added pdf");

    fPdfs.push_back(pdf_->Clone());    
    fNPdfs++;
    fNormalisations.resize(fNPdfs, 0);
}

void
PdfManager::AddPdfs(const std::vector<Pdf*>& pdfs_){
    for(size_t i = 0; i < pdfs_.size(); i++)
        AddPdf(pdfs_.at(i));
}

double
PdfManager::Probability(const EventData& event_) const{
    double prob = 0;
    for(size_t i = 0; i < fPdfs.size(); i++)
        prob += fNormalisations.at(i) * fPdfs.at(i)->Probability(event_);
    
    return prob;
}

const std::vector<double>&
PdfManager::GetNormalisations() const{
    return fNormalisations;
}

void
PdfManager::SetNormalisations(const std::vector<double>& norms_){
    if (norms_.size() != fNPdfs)
        throw DimensionError("PdfManager::SetNormalisations", fNPdfs, 
                             norms_.size());
    fNormalisations = norms_;
}

size_t 
PdfManager::GetNDims() const{
    return fNDims;
}

size_t 
PdfManager::GetNPdfs() const{
    return fNPdfs;
}

// Make a fittable component - i.e. rescale pdfs inside to fit
void
PdfManager::MakeFittable(){
    fParameterManager.Clear();
    fNormalisations.resize(fNPdfs);
    fParameterManager.AddContainer<std::vector<double> >(fNormalisations, 
                                                         "Pdf Normalisation");
}

std::vector<std::string>
PdfManager::GetParameterNames() const{
    return fParameterManager.GetParameterNames();
}
std::vector<double>
PdfManager::GetParameters() const{
    return fParameterManager.GetParameters();
}

size_t 
PdfManager::GetParameterCount() const{
    return fParameterManager.GetParameterCount();
}

void
PdfManager::SetParameters(const std::vector<double>& params_){
    try{
        fParameterManager.SetParameters(params_);
    }
    catch(const ParameterCountError& e_){
        throw ParameterCountError(std::string("PdfManager:: ") + e_.what());
    }
}
