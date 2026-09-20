#include "G4PrimaryGeneratorAction.hh"
#include "G4Constantes.hh"

#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <cmath>
#include <cstdlib>

namespace {
constexpr G4double kUpstreamZ = -70.0 * cm;

std::vector<std::string> SplitCsv(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) {
        fields.push_back(field);
    }
    return fields;
}

const char* InputPath() {
    const char* path = std::getenv("MCP_INPUT_CSV");
    if (!path || !*path) {
        G4Exception("G4PrimaryGeneratorAction::InputPath",
                    "MCPInput001", FatalException,
                    "MCP_INPUT_CSV is not set.");
        return "";
    }
    return path;
}
}

G4PrimaryGeneratorAction::G4PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(),
      particleGun(new G4ParticleGun(1)),
      selectedTrackId(-1),
      selectedDataEventIndex(-1),
      selectedDataGroupId(-1) {
    flag_alpha = true;
    const char* path = InputPath();
    LoadMcpCsv(path);

    G4cout << "Loaded accepted MCP spectra: " << mcpTracks.size()
           << " tracks from " << path << G4endl;
}

G4PrimaryGeneratorAction::~G4PrimaryGeneratorAction() {
    delete particleGun;
}

void G4PrimaryGeneratorAction::LoadMcpCsv(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        G4Exception("G4PrimaryGeneratorAction::LoadMcpCsv",
                    "MCPInput002", FatalException,
                    "Could not open MCP input CSV.");
        return;
    }

    std::string line;
    G4long lineNumber = 0;

    while (std::getline(input, line)) {
        ++lineNumber;
        if (line.empty()) continue;

        const auto fields = SplitCsv(line);
        if (fields.empty()) continue;
        if (fields[0] == "event_index") continue;

        if (fields.size() < 8) {
            G4ExceptionDescription message;
            message << "CSV line " << lineNumber
                    << " has fewer than 8 columns.";
            G4Exception("G4PrimaryGeneratorAction::LoadMcpCsv",
                        "MCPInput003", FatalException, message);
            return;
        }

        try {
            const G4long eventIndex = std::stol(fields[0]);
            const G4int mcpPdg = std::stoi(fields[1]);
            const G4double px = std::stod(fields[2]);
            const G4double py = std::stod(fields[3]);
            const G4double pz = std::stod(fields[4]);
            const G4double energy = std::stod(fields[5]);
            const G4double xDetector = std::stod(fields[6]);
            const G4double yDetector = std::stod(fields[7]);

            const G4ThreeVector momentum(px * GeV, py * GeV, pz * GeV);

            if ((mcpPdg != 1000222 && mcpPdg != -1000222) ||
                !std::isfinite(momentum.mag()) || momentum.mag2() <= 0.0 ||
                !std::isfinite(energy) || energy <= 0.0 ||
                !std::isfinite(xDetector) || !std::isfinite(yDetector)) {
                G4ExceptionDescription message;
                message << "Invalid values on CSV line " << lineNumber;
                G4Exception("G4PrimaryGeneratorAction::LoadMcpCsv",
                            "MCPInput004", FatalException, message);
                return;
            }

            McpTrack track;
            track.trackId = static_cast<G4long>(mcpTracks.size());
            track.eventIndex = eventIndex;
            track.groupId = mcpPdg;
            track.start =
                G4ThreeVector(xDetector * m, yDetector * m, kUpstreamZ);
            track.momentum = momentum;
            track.totalEnergy = energy * GeV;
            track.pythiaPdg = mcpPdg;
            mcpTracks.push_back(track);
        } catch (const std::exception&) {
            G4ExceptionDescription message;
            message << "Could not parse CSV line " << lineNumber;
            G4Exception("G4PrimaryGeneratorAction::LoadMcpCsv",
                        "MCPInput005", FatalException, message);
            return;
        }
    }

    if (mcpTracks.empty()) {
        G4Exception("G4PrimaryGeneratorAction::LoadMcpCsv",
                    "MCPInput006", FatalException,
                    "No MCP tracks were found in the CSV.");
    }
}


void G4PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    const auto eventID = static_cast<std::size_t>(event->GetEventID());
    const auto nTracks = mcpTracks.size();
    const auto replica = eventID / nTracks;
    const auto index = eventID % nTracks;
    if (eventID == 0) {
        const G4long seed = 3911247;
        CLHEP::HepRandom::setTheSeed(seed);
    }

    const auto& track = mcpTracks[index];

    selectedTrackId = track.trackId;
    selectedDataEventIndex = track.eventIndex;
    selectedDataGroupId = track.groupId;

    const char* particleName = track.pythiaPdg < 0 ? "antimillicharged" : "millicharged";
    auto* particle = G4ParticleTable::GetParticleTable()->FindParticle(particleName);
    if (!particle) {
        G4Exception("G4PrimaryGeneratorAction::GeneratePrimaries",
                    "MCPInput007", FatalException,
                    "Could not find the MCP particle definition.");
        return;
    }

    const G4double kineticEnergy = track.totalEnergy - particle->GetPDGMass();
    if (!std::isfinite(kineticEnergy) || kineticEnergy < 0.0) {
        G4Exception("G4PrimaryGeneratorAction::GeneratePrimaries",
                    "MCPInput008", FatalException,
                    "Input total energy is below the configured MCP mass.");
        return;
    }

    particleGun->SetParticleDefinition(particle);
    particleGun->SetParticlePosition(track.start);
    particleGun->SetParticleMomentumDirection(track.momentum.unit());
    particleGun->SetParticleEnergy(kineticEnergy);
    particleGun->GeneratePrimaryVertex(event);

    flag_alpha = true;
}

void G4PrimaryGeneratorAction::SetOptPhotonPolar() {
    SetOptPhotonPolar(G4UniformRand() * 360.0 * deg);
}

void G4PrimaryGeneratorAction::SetOptPhotonPolar(G4double angle) {
    if (particleGun->GetParticleDefinition()->GetParticleName()
        != "opticalphoton") {
        G4cout << "--> warning from PrimaryGeneratorAction::SetOptPhotonPolar(): "
               << "the particleGun is not an opticalphoton" << G4endl;
        return;
    }

    const G4ThreeVector normal(1., 0., 0.);
    const G4ThreeVector direction =
        particleGun->GetParticleMomentumDirection();
    const G4ThreeVector product = normal.cross(direction);
    const G4double modulus2 = product * product;

    G4ThreeVector perpendicular(0., 0., 1.);
    if (modulus2 > 0.) perpendicular = product / std::sqrt(modulus2);
    const G4ThreeVector parallel = perpendicular.cross(direction);

    particleGun->SetParticlePolarization(
        std::cos(angle) * parallel + std::sin(angle) * perpendicular);
}
