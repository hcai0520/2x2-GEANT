#ifndef MyPhysicsList_h
#define MyPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "FTFP_BERT.hh"

class MyPhysicsList : public FTFP_BERT {
public:
    MyPhysicsList();
    ~MyPhysicsList() override = default;

    void ConstructParticle() override;
    void ConstructProcess() override;

private:
    void DefineMillichargedParticle();
};

#endif
