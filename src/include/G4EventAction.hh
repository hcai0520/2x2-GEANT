/// \file G4EventAction.hh
/// \brief Definition of the G4EventAction class

#ifndef G4EventAction_h
#define G4EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4Event.hh"
/// Event action class
///
/// It defines data members to hold the energy deposit and track lengths
/// of charged particles in Absober and Gap layers:
/// - fEnergyAbs, fEnergyGap, fTrackLAbs, fTrackLGap
/// which are collected step by step via the functions
/// - AddAbs(), AddGap()

class G4EventAction : public G4UserEventAction
{
public:
    G4EventAction();
    virtual ~G4EventAction();
    
    virtual void    BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);
    void ResetMuonTrack();
    void UpdateMuonTrack(const G4ThreeVector& prePos,
                         const G4ThreeVector& postPos);
    G4ThreeVector GetMuonStart() const { return muonStart; }
    G4ThreeVector GetMuonEnd()   const { return muonEnd; }
    
public:
    G4int nOfReflections;
    G4int nOfDetections;

private:
    G4ThreeVector muonStart;
    G4ThreeVector muonEnd;
    G4bool hasMuonTrack;
    
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

