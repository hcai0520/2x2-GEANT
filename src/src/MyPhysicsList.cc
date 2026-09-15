#include "MyPhysicsList.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4eMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4CoulombScattering.hh"
#include "G4Exception.hh"

#include <cstdlib>
#include <string>


MyPhysicsList::MyPhysicsList() : FTFP_BERT() {
    ReplacePhysics(new G4EmStandardPhysics_option4());
    RegisterPhysics(new G4OpticalPhysics());
}

void MyPhysicsList::DefineMillichargedParticle() {
    auto* table = G4ParticleTable::GetParticleTable();
    const char* massText = std::getenv("MCP_MASS_GEV");
    if (!massText || !*massText) {
        G4Exception("MyPhysicsList::DefineMillichargedParticle",
                    "MCPMass001", FatalException,
                    "MCP_MASS_GEV is not set.");
        return;
    }
    const G4double mass = std::stod(massText) * CLHEP::GeV;

    if (!table->FindParticle("millicharged")) {
        new G4ParticleDefinition(
            "millicharged",        // name
            mass,                   // mass
            0.0 * CLHEP::MeV,             // width
            0.01*CLHEP::eplus,           // charge
            1,                     // 2*spin
            0, 0,                  // parity, C-conjugation
            0, 0, 0,               // isospin
            "lepton",              // type
            1,                     // lepton number
            0,                     // baryon number
            1000222,               // Pythia MCP PDG code
            true,                  // stable
            -1.0,                  // lifetime
            nullptr,               // decay table
            false,                 // shortlived
            "millicharged",        // subtype
            -1000222                // anti-PDG code
        );
    }

    if (!table->FindParticle("antimillicharged")) {
        new G4ParticleDefinition(
            "antimillicharged",
            mass,
            0.0 * CLHEP::MeV,
            -0.01*CLHEP::eplus,
            1,
            0, 0,
            0, 0, 0,
            "lepton",
            -1,
            0,
            -1000222,
            true,
            -1.0,
            nullptr,
            false,
            "millicharged",
            1000222
        );
    }
}

void MyPhysicsList::ConstructParticle() {
    // Let base FTFP_BERT register all standard particles
    FTFP_BERT::ConstructParticle();
    // Then add mCP
    DefineMillichargedParticle();
}

void MyPhysicsList::ConstructProcess() {
    // Let FTFP_BERT register all standard processes
    FTFP_BERT::ConstructProcess();

    // Add processes for millicharged
    auto* particleIterator = GetParticleIterator();
    particleIterator->reset();

    while ((*particleIterator)()) {
        G4ParticleDefinition* particle = particleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();

        if (particle->GetParticleName() == "millicharged" ||
            particle->GetParticleName() == "antimillicharged") {
            pmanager->AddProcess(new G4eMultipleScattering(), -1, 1, 1);
            pmanager->AddProcess(new G4hIonisation(),         -1, 2, 2);
            pmanager->AddProcess(new G4eBremsstrahlung(),     -1, -1, 3);
            pmanager->AddDiscreteProcess(new G4CoulombScattering());
        }
    }
}
