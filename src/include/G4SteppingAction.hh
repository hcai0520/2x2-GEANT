#ifndef G4SteppingAction_h
#define G4SteppingAction_h 1

#include "G4RunAction.hh"
#include "G4EventAction.hh"
#include "G4PrimaryGeneratorAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4SteppingManager.hh"
#include "G4UserSteppingAction.hh"
#include "G4SystemOfUnits.hh"

class G4EventAction;
class G4RunAction;
class G4PrimaryGeneratorAction;

//Stepping action class
class G4SteppingAction : public G4UserSteppingAction {

   public:
       G4SteppingAction(G4EventAction *eventAction, G4RunAction* runAction, G4PrimaryGeneratorAction* PrimaryGeneratorAction);
	virtual ~G4SteppingAction ();

	virtual void UserSteppingAction (const G4Step*);

   private:
	G4int fEventNumber;
    G4EventAction* fEventAction;
    G4RunAction* fRunAction;
    G4PrimaryGeneratorAction* fPrimaryGeneratorAction;

};
#endif
