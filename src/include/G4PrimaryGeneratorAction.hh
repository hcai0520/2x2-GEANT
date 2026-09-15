#ifndef G4PrimaryGeneratorAction_h
#define G4PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <string>
#include <vector>

class G4ParticleGun;
class G4Event;

class G4PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
  public:
    G4PrimaryGeneratorAction();
    ~G4PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event*) override;

    void SetOptPhotonPolar();
    void SetOptPhotonPolar(G4double);

    G4long GetSelectedTrackId() const { return selectedTrackId; }
    G4long GetSelectedDataEventIndex() const { return selectedDataEventIndex; }
    G4long GetSelectedDataGroupId() const { return selectedDataGroupId; }

    bool flag_alpha;

  private:
    struct McpTrack {
        G4long trackId;
        G4long eventIndex;
        G4long groupId;
        G4ThreeVector start;
        G4ThreeVector momentum;
        G4double totalEnergy;
        G4int pythiaPdg;
    };

    void LoadMcpCsv(const std::string& path);

    G4ParticleGun* particleGun;
    std::vector<McpTrack> mcpTracks;

    G4long selectedTrackId;
    G4long selectedDataEventIndex;
    G4long selectedDataGroupId;
};

#endif
