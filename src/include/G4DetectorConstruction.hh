#ifndef G4DetectorConstruction_h
#define G4DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "DetectorConfig.hh"
class G4VPhysicalVolume;

class G4DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
     G4DetectorConstruction(G4double RIndex, DetectorConfig& GeoConf);
    virtual ~G4DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
  
  
  private:
    // methods
    //
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
    DetectorConfig fConfig;
    // World
    G4double world_x, world_y, world_z;

    G4double VLAr_x, VLAr_y, VLAr_z;

    G4double Pixel_x, Pixel_y, Pixel_z;

    G4double clearance;

    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps

    G4double Refr_Index;
};

#endif /*G4DetectorConstruction_h*/
