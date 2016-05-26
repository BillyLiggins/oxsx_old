#include <ROOTNtuple.h>
#include <EventData.h>
#include <iostream>
#include <Exceptions.h>
#include <TNtuple.h>

ROOTNtuple::ROOTNtuple(const std::string& fileName_, const std::string& treeName_){
    fROOTFile = new TFile(fileName_.c_str());

    if (fROOTFile->IsZombie()){
        delete fROOTFile;
        throw IOError("ROOTNtuple::File Does not Exist! or is Zombie " + fileName_);
    }

    fNtuple = dynamic_cast<TNtuple*>(fROOTFile -> Get(treeName_.c_str()));

    if(!fNtuple){
        delete fROOTFile;
        throw IOError("ROOTNtuple::Tree does not exist, or isn't an ntuple! " + treeName_);
    }        

}

ROOTNtuple::~ROOTNtuple(){
    if (fROOTFile)
        fROOTFile -> Close();
    delete fROOTFile;
}

EventData 
ROOTNtuple::Assemble(size_t iEvent_) const{
    if (iEvent_ >= GetNEntries())
        throw NotFoundError("Exceeded end of ROOT NTuple");

    fNtuple -> GetEntry(iEvent_);
    float* vals = fNtuple -> GetArgs();
    return EventData(std::vector<double> (vals, vals + GetNObservables()));
    
}

EventData
ROOTNtuple::GetEntry(size_t iEvent_) const{
    if(iEvent_ >= GetNEntries())
        throw NotFoundError("Exceeded end of ROOT NTuple");
    return Assemble(iEvent_);
}

std::vector<std::string>
ROOTNtuple::GetObservableNames() const{
    unsigned nObs = GetNObservables();
    std::vector<std::string> names(nObs);
    for(unsigned i = 0; i < nObs; i++)
        names[i] = (fNtuple -> GetListOfBranches()->At(i) -> GetName());
    return names;
    
}

void
ROOTNtuple::LoadBaskets(){
  if(!fNtuple)
    return;
  fNtuple->LoadBaskets();
}


void
ROOTNtuple::DropBaskets(){
  if(!fNtuple)
    return;
  fNtuple->DropBaskets();
}

unsigned
ROOTNtuple::GetNEntries() const{
    return fNtuple->GetEntries();
}

unsigned
ROOTNtuple::GetNObservables() const{
    return fNtuple->GetNvar();
}
