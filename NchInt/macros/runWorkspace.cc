#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "TROOT.h"

#include "createWorkspace.C"

using namespace onia;

//========================================================
// code to read input arguments
template<typename T>
void fromSplit(const std::string& key, const std::string &arg, T& out)
{
  const char delim = '=';
  // Skip if key or delimiter not there
  if ((arg.find(key) == std::string::npos) ||
      (arg.find(delim) == std::string::npos))
    return;

  std::string skey, sval;
  std::stringstream sstr(arg);
  std::getline(sstr, skey, delim); // Dummy read to skip key
  std::getline(sstr, sval, delim); // Get value
  T tout;
  if (!(std::istringstream(sval) >> std::boolalpha >> tout))
    return;
  out = tout;
  std::cout << std::boolalpha << skey << ": "  << out << std::endl;
}

// Special version for string without the conversion 
template<>
void fromSplit(const std::string& key, const std::string &arg, std::string &out)
{
  const char delim = '=';
  // Skip if key or delimiter not there
  if ((arg.find(key) == std::string::npos) ||
      (arg.find(delim) == std::string::npos))
    return;
  std::string skey, sval;
  std::stringstream sstr(arg);
  std::getline(sstr, skey, delim); // Dummy read to skip key
  std::getline(sstr, sval, delim); // Get value
  out = sval;
  std::cout << skey << ": "  << out << std::endl;
}

//=====================================================================
int main(int argc, char* argv[]){

  // set default values
  int nState = 999;
	bool correctCtau  = false;
	bool drawRapPt2D  = false;

  // Loop over argument list                                                                                                                                                       
  for (int i=1; i < argc; i++)
    {
      std::string arg = argv[i];
      fromSplit("nState", arg, nState);
      fromSplit("correctCtau", arg, correctCtau);
      fromSplit("drawRapPt2D", arg, drawRapPt2D);
    }
  
  const std::string infilename = "tmpFiles/selEvents_data.root";
  createWorkspace(infilename.c_str(), nState, correctCtau, drawRapPt2D);
  
  return 0;
}
