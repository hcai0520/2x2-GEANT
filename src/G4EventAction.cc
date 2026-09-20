/// \file G4EventAction.cc
/// \brief Implementation of the G4EventAction class

#include "G4EventAction.hh"
#include "G4RunAction.hh"
#include "G4Constantes.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "MyGlobals.hh"
#include "G4Exception.hh"
#include <limits>
#include <stdexcept>

//================================================================================

G4EventAction::G4EventAction() : G4UserEventAction(), fastModel(FastBlipModel::Instance()) {}


//================================================================================

G4EventAction::~G4EventAction() {}

//================================================================================

void G4EventAction::BeginOfEventAction(const G4Event* /*event*/) {
    
    nOfReflections = 0;
    nOfDetections = 0;
    ResetMuonTrack();
    activePath.clear();
    primaryTrackID = -1;
}

//================================================================================

void G4EventAction::ResetMuonTrack() {
    muonStart = G4ThreeVector(0,0,0);
    muonEnd   = G4ThreeVector(0,0,0);
    hasMuonTrack = false;
}

void G4EventAction::UpdateMuonTrack(const G4ThreeVector& prePos,
                                    const G4ThreeVector& postPos) {
    if (!hasMuonTrack) {
        muonStart = prePos;
        hasMuonTrack = true;
    }

    muonEnd = postPos;
}
//================================================================================
void G4EventAction::RecordFastStep(G4int trackID, const G4ThreeVector& pre,
                                  const G4ThreeVector& post, G4bool active) {
    if (primaryTrackID != -1 && primaryTrackID != trackID) {
        G4Exception("G4EventAction::RecordFastStep", "FastBlip002", FatalException,
                    "Fast blips require exactly one MCP primary per event.");
        throw std::runtime_error("Multiple MCP primaries in fast mode");
    }
    primaryTrackID = trackID;
    if (active) fastModel.AddStep(pre, post, activePath);
}

void G4EventAction::WriteFastBlips(const G4Event* event) {
    FastBlipModel::Pair pair;
    if (!event->IsAborted()) pair = fastModel.Sample(activePath);
    else for (const auto& segment : activePath) pair.length += segment.length;
    const G4int status = event->IsAborted() ? 2 : (pair.valid ? 1 : 0);
    const G4double missing = std::numeric_limits<G4double>::quiet_NaN();
    auto analysis = G4AnalysisManager::Instance();
    const auto run = G4RunManager::GetRunManager()->GetCurrentRun();
    analysis->FillNtupleIColumn(0, 0, run->GetRunID());
    analysis->FillNtupleIColumn(0, 1, event->GetEventID());
    analysis->FillNtupleIColumn(0, 2, status);
    analysis->FillNtupleIColumn(0, 3, primaryTrackID);
    analysis->FillNtupleIColumn(0, 4, static_cast<G4int>(activePath.size()));
    const G4double values[] = {
        pair.length / CLHEP::cm, pair.mean, pair.weight,
        pair.valid ? pair.s1 / CLHEP::cm : missing,
        pair.valid ? pair.s2 / CLHEP::cm : missing,
        pair.valid ? pair.first.x() / CLHEP::cm : missing,
        pair.valid ? pair.first.y() / CLHEP::cm : missing,
        pair.valid ? pair.first.z() / CLHEP::cm : missing,
        pair.valid ? pair.second.x() / CLHEP::cm : missing,
        pair.valid ? pair.second.y() / CLHEP::cm : missing,
        pair.valid ? pair.second.z() / CLHEP::cm : missing,
        pair.valid ? (pair.second - pair.first).mag() / CLHEP::cm : missing
    };
    for (G4int i = 0; i < 12; ++i) analysis->FillNtupleDColumn(0, i + 5, values[i]);
    // Include no-path and aborted events; never normalize by candidate rows alone.
    analysis->AddNtupleRow(0);
}

void G4EventAction::EndOfEventAction(const G4Event* event) {
    if (FastBlipsEnabled()) {
        WriteFastBlips(event);
        return;
    }
    
    // Print number of reflections
    if (IMPRIMIR_EVENTOS == true) {
        G4cout << "Number of reflections [Event] = " << nOfReflections << G4endl;
        G4cout << "Number of detections [Event] = " << nOfDetections << G4endl;
    }
    auto analysisManager = G4AnalysisManager::Instance();

    G4int eventID = event->GetEventID();
    G4ThreeVector startPos = GetMuonStart();
    G4ThreeVector endPos   = GetMuonEnd();
    G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
    if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
        G4cout << "---> End of event: " << eventID << G4endl;
    }
    for (const auto& volumeEntry : accumulatedEnergy) {
        const G4String& volumeName = volumeEntry.first;
        const std::map<G4int, G4double>& pixelMap = volumeEntry.second;

        for (const auto& entry : pixelMap) {
            G4int pixelID = entry.first;
            G4double energyDeposit = entry.second; // In MeV
            if (energyDeposit <= 0.0) continue;

            G4int volname = -1;
            if      (volumeName == "Prisms_M0") volname = 2;
            else if (volumeName == "Prisms_M1") volname = 0;
            else if (volumeName == "Prisms_M2") volname = 3;
            else if (volumeName == "Prisms_M3") volname = 1;

            const G4int nPixelY = 280;
            const G4int nPixelZ = 140;
            const G4int pixelsPerIOGroup = nPixelY * nPixelZ;
            const G4int plane = pixelID / pixelsPerIOGroup;
            const G4int thresholdLocalKey = pixelID;
            const G4int rem = pixelID % pixelsPerIOGroup;
            const G4int pixelZ = rem % nPixelZ;
            const G4int pixelY = rem / nPixelZ;
            const G4int ioGroup = (volname >= 0 && plane >= 0 && plane < 2)
                                      ? (2 * volname + plane + 1)
                                      : -1;
            
            //primary pixel energy
            //secondary energy
            G4double ePrim = 0.0, eSec = 0.0;
            if (auto vp = accumulatedEnergy_Primary.find(volumeName); vp != accumulatedEnergy_Primary.end()) {
                if (auto ip = vp->second.find(pixelID); ip != vp->second.end()) ePrim = ip->second;
            }
            if (auto vs = accumulatedEnergy_Secondary.find(volumeName); vs != accumulatedEnergy_Secondary.end()) {
                if (auto is = vs->second.find(pixelID); is != vs->second.end()) eSec = is->second;
            }



            const bool hasPrim = (ePrim > 0.0);
            const bool hasSec  = (eSec  > 0.0);

            const G4double primaryPixelEnergy   = hasPrim ? ePrim : 0.0;
            const G4double secondaryOnlyEnergy  = (!hasPrim && hasSec) ? energyDeposit : 0.0;
            const G4double primaryOnlyEnergy  = (hasPrim && !hasSec) ? energyDeposit : 0.0;

            //position y and z
            G4double x=0.0, y = 0.0, z = 0.0;
            if (auto vx = positionX.find(volumeName); vx != positionX.end()) {
                if (auto ix = vx->second.find(pixelID); ix != vx->second.end()) x = ix->second;
            }
            if (auto vy = positionY.find(volumeName); vy != positionY.end()) {
                if (auto iy = vy->second.find(pixelID); iy != vy->second.end()) y = iy->second;
            }
            if (auto vz = positionZ.find(volumeName); vz != positionZ.end()) {
                if (auto iz = vz->second.find(pixelID); iz != vz->second.end()) z = iz->second;
            }
            analysisManager->FillNtupleIColumn(0, 0, eventID);
            analysisManager->FillNtupleIColumn(0, 1, pixelID);
            analysisManager->FillNtupleIColumn(0, 2, volname);

            analysisManager->FillNtupleDColumn(0, 3, energyDeposit / CLHEP::MeV);
            analysisManager->FillNtupleDColumn(0, 4, primaryPixelEnergy / CLHEP::MeV);
            analysisManager->FillNtupleDColumn(0, 5, primaryOnlyEnergy / CLHEP::MeV);
            analysisManager->FillNtupleDColumn(0, 6, secondaryOnlyEnergy / CLHEP::MeV);

            analysisManager->FillNtupleDColumn(0, 7, x / CLHEP::cm);
            analysisManager->FillNtupleDColumn(0, 8, y / CLHEP::cm);
            analysisManager->FillNtupleDColumn(0, 9, z / CLHEP::cm);

            analysisManager->FillNtupleDColumn(0,10, startPos.x() / CLHEP::cm);
            analysisManager->FillNtupleDColumn(0,11, startPos.y() / CLHEP::cm);
            analysisManager->FillNtupleDColumn(0,12, startPos.z() / CLHEP::cm);

            analysisManager->FillNtupleDColumn(0,13, endPos.x() / CLHEP::cm);
            analysisManager->FillNtupleDColumn(0,14, endPos.y() / CLHEP::cm);
            analysisManager->FillNtupleDColumn(0,15, endPos.z() / CLHEP::cm);

            analysisManager->FillNtupleIColumn(0,16, ioGroup);
            analysisManager->FillNtupleIColumn(0,17, plane);
            analysisManager->FillNtupleIColumn(0,18, pixelY);
            analysisManager->FillNtupleIColumn(0,19, pixelZ);
            analysisManager->FillNtupleIColumn(0,20, thresholdLocalKey);
            analysisManager->AddNtupleRow(0);            
        }
    }

    accumulatedEnergy_Primary.clear();
    accumulatedEnergy_Secondary.clear();
    accumulatedEnergy.clear();
    positionY.clear();
    positionZ.clear();
    positionX.clear();
}

//================================================================================
