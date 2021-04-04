#ifndef SWITCHES_H
#define SWITCHES_H

//#include "structs.h"
//#include "switches.cxx"
class switches
{
    private:
        std::vector<SwitchInfo> switchList;
    public:
        switches();
        ~switches();
        bool checkToSaveSpectrum(std::string key);
        void InsertSwitch(SwitchInfo Switch);
    
    //ClassDef(switches,1)
};
//#if !defined(__CINT__)
//ClassImp(switches);
//#endif
#endif //SWITCHES_H