#ifndef SWITCHES_H
#define SWITCHES_H

class switches
{
    private:
        std::vector<SwitchInfo> switchList;
    public:
        switches();
        ~switches();
        bool checkToSaveSpectrum(std::string key);
        void InsertSwitch(SwitchInfo Switch);
};
#endif //SWITCHES_H