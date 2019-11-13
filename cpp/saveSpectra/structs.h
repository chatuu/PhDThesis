#include "headers.h"

struct Spectra
{
    std::string intName;
    std::string varName;
    std::string cutName;
    Spectrum *spectrum;
};


struct Cuts
{
    std::string name;
    Cut cut;
};

struct TruthCuts
{
    std::string name;
    NuTruthCut cut;
};

struct Vars
{   
    std::string name;
    HistAxis axis;
};

struct TruthVars
{   
    std::string name;
    NuTruthHistAxis axis;
};