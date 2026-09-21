#ifndef FastBlipModel_h
#define FastBlipModel_h 1

#include "G4ThreeVector.hh"
#include "globals.hh"
#include <vector>

// An effective visible-blip Poisson model, not microscopic interaction forcing.
class FastBlipModel {
public:
    struct Segment {
        G4ThreeVector start, end;
        G4double length; // Geometric chord length in Geant4 units.
    };
    struct Pair {
        G4ThreeVector first, second;
        G4double length = 0., mean = 0., weight = 0.;
        G4double s1 = 0., s2 = 0.;
        G4bool valid = false;
    };

    static const FastBlipModel& Instance();
    G4bool Enabled() const { return enabled; }
    void AddStep(const G4ThreeVector& pre, const G4ThreeVector& post,
                 std::vector<Segment>& path) const;
    Pair Sample(const std::vector<Segment>& path) const;

private:
    FastBlipModel();
    G4double Mean(G4double length) const;
    G4double InverseMean(G4double mean) const;
    G4bool enabled = false;
    G4double rate = 0.; // Per Geant4 length unit; used without a profile.
    std::vector<G4double> lengths, means;
};
#endif
