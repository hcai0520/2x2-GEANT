#ifndef G4RunAction_h
#define G4RunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4RootAnalysisManager.hh"
#include "DetectorConfig.hh"

#include "G4NeutronHPBuilder.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Timer;
class G4Run;

class G4RunAction : public G4UserRunAction
{
  public:
    G4RunAction(const DetectorConfig& config);
   ~G4RunAction();
  //static G4NeutronHPManager *   GetInstance ();


  public:
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);
    G4int nOfReflections_Total;
    G4int nOfDetections_Total;
// Soma do tempo de voo de todos os fotons detectados
    G4double TOF_Detections_Total;

  private:
    G4Timer* timer;
    const DetectorConfig fConfig;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*G4RunAction_h*/
