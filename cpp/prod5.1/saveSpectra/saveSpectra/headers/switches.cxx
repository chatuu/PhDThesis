//#include "headers.h"
#include "switches.h"

switches::switches()
{

}
switches::~switches()
{

}

void switches::InsertSwitch(SwitchInfo Switch)
{
    switchList.push_back(Switch);

}

bool switches::checkToSaveSpectrum(std::string key)
{
    for (auto it = switchList.begin(); it != switchList.end(); ++it)
    {
        std::string testString = it->interaction+it->variable+it->cut;
        if (testString == key)
            return true;
        else
            continue;
    }
    return false;
}
    
