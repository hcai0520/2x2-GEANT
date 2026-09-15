#ifndef G4ActionInitialization_h
#define G4ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "DetectorConfig.hh"
class G4DetectorConstruction;

/// Action initialization class.
class G4ActionInitialization : public G4VUserActionInitialization {
    
  public:
    G4ActionInitialization(G4DetectorConstruction*, DetectorConfig& config);
    virtual ~G4ActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

  private:
    G4DetectorConstruction* fDetConstruction;
    DetectorConfig fConfig;
};

#endif

    
