#ifndef __RooUtils_C__
#define __RooUtils_C__

#include "TFile.h"
#include "RooWorkspace.h"

using namespace RooFit;

RooRealVar*  getVarFromWorkspace(RooWorkspace* ws, std::string varName, double& value, double &error) {
    assert ( 0 != ws );
    RooRealVar* var = ws->var(varName.c_str());
    assert ( 0 != var );
    value = var->getVal();
    error = var->getError();
    return var;
}


/// Get an object from a TFile
/// @param fileName Name of TFile to open
/// @param keName Name of object in TFile to get  
/// @return fetched object
/// @throws An error string if file not there or key not available
template <class T>
T* getFromTFile(const std::string &fileName, const std::string &keyName) {
    TFile *inFile = new TFile(fileName.c_str(), "R");
    if (inFile->IsZombie()) {
	throw std::string("Error opening file " + fileName);
    }
    T *t = dynamic_cast<T*>(inFile->Get(keyName.c_str()));
    if (t == 0) {
	throw std::string("Cannot get object " + keyName + " from file " + fileName);
    }
    delete inFile;
    return t;
}

#endif
