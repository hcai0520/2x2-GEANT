#include "G4Timer.hh"

#include "G4RunAction.hh"
#include "G4Constantes.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4RootAnalysisManager.hh"
#include "G4Analyser.hh"
#include "G4AnalysisManager.hh"
#include <fstream>
#include "FastBlipModel.hh"
#include <string>

// #include "G4NeutronHPManager.hh"

using namespace std;

//================================================================================

G4RunAction::G4RunAction(const DetectorConfig& config) : G4UserRunAction(), fConfig(config) {

    
    timer = new G4Timer;
    nOfReflections_Total = 0;
    nOfDetections_Total = 0;
    TOF_Detections_Total = 0;
}

//================================================================================

G4RunAction::~G4RunAction() {
    
    delete timer;
  
}

//================================================================================

void G4RunAction::BeginOfRunAction(const G4Run* aRun) {
    
    // G4NeutronHPManager::GetInstance()->SetVerboseLevel(0);


    G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

    auto analysisManager = G4AnalysisManager::Instance();
    
    G4cout << "Using " << analysisManager->GetType() << " analysis manager." << G4endl;
    timer->Start();
    analysisManager->SetDefaultFileType("root");
    if (FastBlipModel::Instance().Enabled()) {
        // Keep separate files for multiple beamOn commands in the same process.
        const G4String filename = aRun->GetRunID() == 0 ? "mcp_fast"
            : "mcp_fast_run" + std::to_string(aRun->GetRunID());
        analysisManager->SetVerboseLevel(0);
        analysisManager->OpenFile(filename);
        if (!fastNtupleBooked) {
            analysisManager->CreateNtuple("FastBlips", "Weighted two-point MCP response");
            for (const auto name : {"RunID", "EventID", "Status", "TrackID", "ActiveSegments"})
                analysisManager->CreateNtupleIColumn(name);
            for (const auto name : {"ActiveLength_cm", "MeanBlips", "PairWeight",
                                   "s1_cm", "s2_cm", "x1_cm", "y1_cm", "z1_cm",
                                   "x2_cm", "y2_cm", "z2_cm", "Distance_cm"})
                analysisManager->CreateNtupleDColumn(name);
            analysisManager->FinishNtuple();
            analysisManager->CreateNtuple("RunSummary", "Geant4 event accounting");
            analysisManager->CreateNtupleIColumn("RunID");
            analysisManager->CreateNtupleIColumn("GeantEventCount");
            analysisManager->FinishNtuple();
            fastNtupleBooked = true;
        }
        return;
    }
    analysisManager->OpenFile("mcp_02");
    analysisManager->SetVerboseLevel(0);


    // Indica o id do primeiro histograma criado -- default = 0
    analysisManager->SetFirstHistoId(1);
    // Cria histogramas    -- (nome, titulo, nbins, xmin, xmax)
    analysisManager->CreateH1("0","Time of arrival",5000,0.,10);


    // Indica o id da primeira ntuple criada -- default = 0
    //analysisManager->SetFirstNtupleId(1);
    //Declara ntuples
    analysisManager->CreateNtuple("PixelEDep", "Energy per pixel");

    analysisManager->CreateNtupleIColumn("EventID");
    analysisManager->CreateNtupleIColumn("PixelID");
    analysisManager->CreateNtupleIColumn("VolName");

    analysisManager->CreateNtupleDColumn("EnergyDep");
    analysisManager->CreateNtupleDColumn("EnergyDep_Primary");
    analysisManager->CreateNtupleDColumn("EnergyDep_Primary_only");
    analysisManager->CreateNtupleDColumn("EnergyDep_Secondary_only");

    analysisManager->CreateNtupleDColumn("x");
    analysisManager->CreateNtupleDColumn("y");
    analysisManager->CreateNtupleDColumn("z");

    analysisManager->CreateNtupleDColumn("x_start");
    analysisManager->CreateNtupleDColumn("y_start");
    analysisManager->CreateNtupleDColumn("z_start");

    analysisManager->CreateNtupleDColumn("x_end");
    analysisManager->CreateNtupleDColumn("y_end");
    analysisManager->CreateNtupleDColumn("z_end");

    analysisManager->CreateNtupleIColumn("IOGroup");
    analysisManager->CreateNtupleIColumn("Plane");
    analysisManager->CreateNtupleIColumn("PixelY");
    analysisManager->CreateNtupleIColumn("PixelZ");
    analysisManager->CreateNtupleIColumn("ThresholdLocalKey");

    analysisManager->FinishNtuple();


    //analysisManager->CreateNtuple("ntuple", "data");
    //analysisManager->CreateNtupleIColumn("evt");
    //analysisManager->CreateNtupleDColumn("x");
    //analysisManager->CreateNtupleDColumn("y");
    //analysisManager->CreateNtupleDColumn("z");
    //analysisManager->CreateNtupleDColumn("t");
    //analysisManager->CreateNtupleDColumn("dedx");
    //analysisManager->CreateNtupleDColumn("theta");
    //analysisManager->CreateNtupleIColumn("VolName");
    //analysisManager->CreateNtupleDColumn("dy");
    //analysisManager->CreateNtupleDColumn("dz");
    //analysisManager->CreateNtupleDColumn("momentum");
   // analysisManager->CreateNtupleIColumn("copyno");
   // analysisManager->CreateNtupleDColumn("accumlated_en");


    //analysisManager->FinishNtuple();
}

//================================================================================

void G4RunAction::EndOfRunAction(const G4Run* aRun) {
    G4cout << "End of run." << G4endl;
    timer->Stop();

    auto analysisManager = G4AnalysisManager::Instance();
    if (FastBlipModel::Instance().Enabled()) {
        analysisManager->FillNtupleIColumn(1, 0, aRun->GetRunID());
        analysisManager->FillNtupleIColumn(1, 1, aRun->GetNumberOfEvent());
        analysisManager->AddNtupleRow(1);
    }
    analysisManager->Write();
    analysisManager->CloseFile();

    //delete G4AnalysisManager::Instance();
    //analysisManager->Clear();    
}

//================================================================================
