#include "G4MTRunManager.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4DetectorConstruction.hh"
#include "G4ActionInitialization.hh"
#include "G4PhysicsList.hh"

#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "DetectorConfig.hh"
#include "MyPhysicsList.hh"
#include "G4SystemOfUnits.hh"
// ============================================================================

int main(int argc,char** argv) {
    
    G4double RIndex = 1;
    if (argc==2) {
      RIndex = atof(argv [1]);
    }
    if (argc==3) {
      RIndex = atof(argv [2]);
    }
    // sizeX, sizeY, sizeZ;
    // pixelSizeY, pixelSizeZ;
    DetectorConfig GeoConf = {60,140,60,.4,.4};
    // ============================================================================
    // Run manager
    // this first step sets the number of cores available.
    G4RunManager* runManager = new G4RunManager;
    //runManager->SetNumberOfThreads(1);
    runManager->SetVerboseLevel(0);

    // declares the geometry of the simulation
    G4DetectorConstruction* detConstruction = new G4DetectorConstruction(RIndex, GeoConf );
    runManager->SetUserInitialization(detConstruction);

    // declares the physics of the simulation
    //G4VModularPhysicsList* physicsList = new FTFP_BERT();
    //physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
    //G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();

    //physicsList->RegisterPhysics(opticalPhysics);
    //runManager->SetUserInitialization(physicsList);
    runManager->SetUserInitialization(new MyPhysicsList());


    // declares the first action of the simulation
    G4ActionInitialization* actionInitialization = new G4ActionInitialization(detConstruction, GeoConf);
    runManager->SetUserInitialization(actionInitialization);
    
    // visualization manager
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
    
    // Initialize G4 kernel
    runManager->Initialize();
    
    // ============================================================================
    // User Interface
    
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    
    // Define UI session for interactive mode
    if (argc < 3) {

      G4UIExecutive * ui = new G4UIExecutive(argc,argv);

      UImanager->ApplyCommand("/control/execute vis.mac");

      ui->SessionStart();
      delete ui;

    } else if (argc == 3) {
      
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
      
    }
    
    // ============================================================================
    // Delete

    delete visManager;

    delete runManager;

    return 0;
}

// ============================================================================
