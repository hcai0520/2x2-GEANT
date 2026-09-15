#include "G4ActionInitialization.hh"
#include "G4DetectorConstruction.hh"
#include "G4EventAction.hh"

#include "G4PrimaryGeneratorAction.hh"
#include "G4RunAction.hh"
#include "G4SteppingAction.hh"

// ============================================================================

G4ActionInitialization::G4ActionInitialization
                                (G4DetectorConstruction* detConstruction, DetectorConfig& config)
: G4VUserActionInitialization(), fDetConstruction(detConstruction),fConfig(config)
{;}

// ============================================================================

G4ActionInitialization::~G4ActionInitialization() {;}

// ============================================================================

void G4ActionInitialization::BuildForMaster() const {
    
    SetUserAction(new G4RunAction(fConfig));
    
}

// ============================================================================

void G4ActionInitialization::Build() const {
    
    auto primaryGenAction = new G4PrimaryGeneratorAction();
    auto runAction = new G4RunAction(fConfig);
    auto eventAction = new G4EventAction();
    auto steppingAction = new G4SteppingAction(eventAction, runAction, primaryGenAction);

    SetUserAction(primaryGenAction);
    SetUserAction(runAction);
    SetUserAction(eventAction);
    SetUserAction(steppingAction);
} 

// ============================================================================
