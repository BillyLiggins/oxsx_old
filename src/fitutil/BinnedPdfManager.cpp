#include <BinnedPdfManager.h>
#include <SystematicManager.h>
#include <BinnedPdf.h>
#include <PdfConverter.h>
#include <Exceptions.h>
#include <sstream>


unsigned 
BinnedPdfManager::GetNPdfs() const{
    return fOriginalPdfs.size();
}

size_t
BinnedPdfManager::GetNDims() const{
    return fNDims;
}

double 
BinnedPdfManager::Probability(const EventData& data_) const{
    double sum = 0;

    for(size_t i = 0; i < fWorkingPdfs.size(); i++){
        sum += fNormalisations.at(i) * fWorkingPdfs[i].Probability(data_);
    }

    return sum;
}

double
BinnedPdfManager::BinProbability(size_t bin_) const{
    double sum = 0;
    try{
        for(size_t i = 0; i < fWorkingPdfs.size(); i++){
            sum += fNormalisations.at(i) * fWorkingPdfs.at(i).GetBinContent(bin_);
    
        }
    }

    catch(const std::out_of_range&){
        throw LogicError("BinnedPdfManager:: Normalisation vector doesn't match pdf vector - are the normalisations set?");
    }
    return sum;
}


void
BinnedPdfManager::SetNormalisations(const std::vector<double>& normalisations_){
    if (normalisations_.size() != fOriginalPdfs.size())
        throw LogicError("BinnedPdfManager: number of norms doesn't match #pdfs");
    fNormalisations = normalisations_;
}

void
BinnedPdfManager::ApplySystematics(const SystematicManager& sysMan_){
    // If there are no systematics or they haven't changed, dont transform the working pdfs 
    //  ( = original pdfs from initialisagtion)
    
    // FIXME: Do you need some equivilent to the code below for saving time on a grid search??
//     if((!sysMan_.GetSystematics().size()) || sysMan_.GetParameters() == fCachedParams)
//         return;
//    fCachedParams = sysMan_.GetParameters();
    if(!sysMan_.GetSystematics().size())
        return;
    
    for(size_t j = 0; j < fOriginalPdfs.size(); j++)
      fWorkingPdfs[j] = sysMan_.GetTotalResponse().operator()(fOriginalPdfs[j]);
}

const BinnedPdf&
BinnedPdfManager::GetOriginalPdf(size_t index_) const{
    return fOriginalPdfs.at(index_);
}

void
BinnedPdfManager::AddPdf(const BinnedPdf& pdf_){
    fOriginalPdfs.push_back(pdf_);
    fWorkingPdfs.push_back(pdf_);
    fNPdfs++;
}

void 
BinnedPdfManager::AddPdfs(const std::vector<BinnedPdf>& pdfs_){
    for(size_t i = 0; i < pdfs_.size(); i++){
        AddPdf(pdfs_.at(i));
    }
    fNPdfs += pdfs_.size();
}

const std::vector<double>&
BinnedPdfManager::GetNormalisations() const{
    return fNormalisations;
}

void
BinnedPdfManager::ApplyShrink(const BinnedPdfShrinker& shrinker_){
    if (!shrinker_.GetBuffers().size())
        return;
        
    // only shrink if not already shrunk! FIXME: more obvious behaviour
    if (!fWorkingPdfs.size() || fWorkingPdfs.at(0).GetNBins() != fOriginalPdfs.at(0).GetNBins())
        return;

    for (size_t i = 0; i < fWorkingPdfs.size(); i++){
        fWorkingPdfs[i] = shrinker_.ShrinkPdf(fWorkingPdfs.at(i));
        fWorkingPdfs[i].Normalise();
    }
    
}


////////////////////////////////
// Make this a fit component! //
////////////////////////////////

void
BinnedPdfManager::MakeFittable(){
    fParameterManager.Clear();
    if(fNormalisations.size() < fNPdfs)
        fNormalisations.resize(fNPdfs, 0);
    fParameterManager.AddContainer(fNormalisations, "Pdf Normalisation");
}

std::vector<std::string>
BinnedPdfManager::GetParameterNames() const {
    return fParameterManager.GetParameterNames();
}

std::vector<double>
BinnedPdfManager::GetParameters() const{
    return fParameterManager.GetParameters();
}

size_t
BinnedPdfManager::GetParameterCount() const{
    return fParameterManager.GetParameterCount();
}

void
BinnedPdfManager::SetParameters(const std::vector<double>& params_){
    try{
        fParameterManager.SetParameters(params_);
    }
    catch(const ParameterCountError& e_){
        throw ParameterCountError("BinnedPdfManager:: " + 
                                  std::string(e_.what())
                                  );
    }
}
