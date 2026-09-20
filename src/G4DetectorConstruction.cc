    #include "G4DetectorConstruction.hh"
    #include "G4Constantes.hh"

    #include "G4Material.hh"
    #include "G4Element.hh"
    #include "G4NistManager.hh"

    #include "G4Box.hh"
    #include "G4Orb.hh"
    #include "G4Tubs.hh"
    #include "G4Sphere.hh"
    #include "G4Trd.hh"

    #include "G4LogicalVolume.hh"
    #include "G4ThreeVector.hh"
    #include "G4PVPlacement.hh"
    #include "G4AutoDelete.hh"
    #include "G4LogicalBorderSurface.hh"
    #include "G4LogicalSkinSurface.hh"
    #include "G4OpticalSurface.hh"

    #include "G4SubtractionSolid.hh"
    #include "G4VSolid.hh"
    #include "G4UnionSolid.hh"
    #include "G4VPVParameterisation.hh"
    #include "G4PVParameterised.hh"


    #include "G4VisAttributes.hh"
    #include "G4Colour.hh"

    #include "G4PhysicalConstants.hh"
    #include "G4SystemOfUnits.hh"
    #include "globals.hh"

    #include "G4Navigator.hh"
    #include "G4VPhysicalVolume.hh"
    #include "G4TransportationManager.hh"

    #include <cmath>
    #include <fstream>
    using namespace std;

    #include "G4NeutronHPBuilder.hh"
    #include <G4HadronicProcessStore.hh>
    #include "G4UserLimits.hh"
    #include "G4ProductionCuts.hh"
    #include "G4Region.hh"
    #include "G4RegionStore.hh"
    //================================================================================

    G4DetectorConstruction::G4DetectorConstruction (G4double RIndex, DetectorConfig& GeoConf)
    : G4VUserDetectorConstruction(), fCheckOverlaps(true), Refr_Index(RIndex), fConfig(GeoConf) {

        //G4NeutronHPManager::GetInstance()->SetVerboseLevel(0);
        G4HadronicProcessStore::Instance()->SetVerbose(0);
        // World
        world_x = 600 * cm * 0.5;
        world_y = 600 * cm * 0.5;
        world_z = 600 * cm * 0.5;
        //total size of the volume
        //VLAr_x =  GeoConf.sizeX * cm * 0.5;
        //VLAr_y =  GeoConf.sizeY * cm * 0.5;
        //VLAr_z =  GeoConf.sizeZ * cm * 0.5;
        // Pixel geometry is fixed in DefineVolumes() from the detector module
        // dimensions and the real threshold-map indexing.
        Pixel_x = 0.0;
        Pixel_y = 0.0;
        Pixel_z = 0.0;


    }

    //================================================================================

    G4DetectorConstruction::~G4DetectorConstruction (){;}




   
    class Full3DParameterisation : public G4VPVParameterisation {
    public:
        Full3DParameterisation(G4int nPlanes,
                               G4int nY,
                               G4int nZ,
                               G4double pitchX,
                               G4double pitchY,
                               G4double pitchZ)
            : fNPlanes(nPlanes),
              fNY(nY),
              fNZ(nZ),
              fPitchX(pitchX),
              fPitchY(pitchY),
              fPitchZ(pitchZ) {}

        void ComputeTransformation(G4int copyNo, G4VPhysicalVolume* physVol) const override {
            const G4int pixelsPerPlane = fNY * fNZ;
            G4int plane = copyNo / pixelsPerPlane;
            G4int rem = copyNo % pixelsPerPlane;
            G4int iz = rem % fNZ;
            G4int iy = rem / fNZ;

            G4double x = (-fNPlanes/2.0 + plane + 0.5) * fPitchX;
            G4double y = (-fNY/2.0 + iy + 0.5) * fPitchY;
            G4double z = (-fNZ/2.0 + iz + 0.5) * fPitchZ;

            physVol->SetTranslation(G4ThreeVector(x, y, z));
        }

        void ComputeDimensions(G4Box& box, G4int, const G4VPhysicalVolume*) const override {
            // box is already defined — nothing to do here if you're using a fixed size
        }

    private:
        G4int fNPlanes;
        G4int fNY;
        G4int fNZ;
        G4double fPitchX;
        G4double fPitchY;
        G4double fPitchZ;
    };

    //================================================================================

    G4VPhysicalVolume* G4DetectorConstruction::Construct () {

        // Define materials 
        DefineMaterials();
        // Define volumes
        return DefineVolumes();
        
    }

    //================================================================================

    void G4DetectorConstruction::DefineMaterials() {
        
        G4NistManager* nistManager = G4NistManager::Instance();

        nistManager->FindOrBuildMaterial("G4_lAr");
        nistManager->FindOrBuildMaterial("G4_Galactic");
        nistManager->FindOrBuildMaterial("G4_URANIUM_MONOCARBIDE");

        G4cout << *(G4Material::GetMaterialTable()) << G4endl;
        
    }

    //================================================================================

    G4VPhysicalVolume* G4DetectorConstruction::DefineVolumes() {
        // Get materials
       
        G4Material* Vacuo      = G4Material::GetMaterial("G4_Galactic");       
        G4Material* LAr        = G4Material::GetMaterial("G4_lAr");
        G4Material* HighSP     = G4Material::GetMaterial("G4_URANIUM_MONOCARBIDE");

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Construction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //====================== World ===========-===========

        G4Box*              WorldS          = new G4Box (MUNDO_NOME, world_x, world_y, world_z);

        G4LogicalVolume*    WorldLV         = new G4LogicalVolume (WorldS, Vacuo, MUNDO_NOME);

        G4VPhysicalVolume*  WorldPV         = new G4PVPlacement (0, G4ThreeVector (), WorldLV, MUNDO_NOME, 0, true, 0, fCheckOverlaps);

    //======================= Parametric volumes =======================


        const G4double mod_half_x = 30.4310 * cm;
        const G4double mod_half_y = 61.8543 * cm;
        const G4double mod_half_z = 30.8163 * cm;

        const G4double cx = 33.5000 * cm;
        const G4double cz = 33.5000 * cm;

        const G4int nPlanes = 2;
        const G4int nY = 280;
        const G4int nZ = 140;
        const G4int pixelsPerPlane = nY * nZ;
        const G4int nTotal = nPlanes * pixelsPerPlane;

        G4double pitchX = (2.0 * mod_half_x) / nPlanes;
        G4double pitchY = (2.0 * mod_half_y) / nY;
        G4double pitchZ = (2.0 * mod_half_z) / nZ;

        const G4double active_half_y = 0.5 * nY * pitchY;
        const G4double active_half_z = 0.5 * nZ * pitchZ;

        G4cout << "2x2 detector geometry:" << G4endl
               << "  module centers: (+/-" << cx / cm << ", 0, +/-" << cz / cm << ") cm" << G4endl
               << "  module half-size: x=" << mod_half_x / cm
               << " cm, y=" << mod_half_y / cm
               << " cm, z=" << mod_half_z / cm << " cm" << G4endl
               << "  IO-group x pitch: " << pitchX / cm << " cm" << G4endl
               << "  pixel pitch: y=" << pitchY / cm
               << " cm, z=" << pitchZ / cm << " cm" << G4endl
               << "  pixel grid per IO group: nY=" << nY
               << ", nZ=" << nZ
               << ", copies/IO group=" << pixelsPerPlane << G4endl
               << "  IO groups per module=" << nPlanes
               << ", copies/module=" << nTotal << G4endl
               << "  active half-size from pixels: y=" << active_half_y / cm
               << " cm, z=" << active_half_z / cm << " cm" << G4endl;

        auto param = new Full3DParameterisation(nPlanes, nY, nZ, pitchX, pitchY, pitchZ);
    //G4int nY = VLAr_y/Pixel_y;
    //G4int nZ = VLAr_z/Pixel_z;
    //G4int nTotal = nY * nZ;//copy number from 0 to fnY*fnZ - 1    
    
    //G4float gapsize  = 35*cm;
    
    // Slice volume (same dimensions for all copies)
    auto sliceSolid = new G4Box("Slice", 0.5 * pitchX, 0.5 * pitchY, 0.5 * pitchZ);
    auto sliceLogic = new G4LogicalVolume(sliceSolid, LAr, "Slice");
    
    
    // M0_Mother volume for stack
    auto M0_motherSolid = new G4Box("M0_Mother", mod_half_x, mod_half_y, mod_half_z);
    auto M0_motherLogic = new G4LogicalVolume(M0_motherSolid, Vacuo, "M0_Mother");
    new G4PVPlacement(0, G4ThreeVector(+cx, 0, +cz),M0_motherLogic, "M0_Mother", WorldLV, false, 0, fCheckOverlaps);


    new G4PVParameterised("Prisms_M0",
                          sliceLogic,
                          M0_motherLogic,
                          kUndefined,  // Not tied to a single axis
                          nTotal,
                          param);
    

    //M1_Mother volume for stack
    auto M1_motherSolid = new G4Box("M1_Mother", mod_half_x, mod_half_y, mod_half_z);
    auto M1_motherLogic = new G4LogicalVolume(M1_motherSolid, Vacuo, "M1_Mother");
    new G4PVPlacement(0, G4ThreeVector(+cx, 0, -cz), M1_motherLogic, "M1_Mother", WorldLV, false, 0, fCheckOverlaps);


    new G4PVParameterised("Prisms_M1",
                          sliceLogic,
                          M1_motherLogic,
                          kUndefined,  // Not tied to a single axis
                          nTotal,
                          param);

     // M2_Mother volume for stack
    auto M2_motherSolid = new G4Box("M2_Mother", mod_half_x, mod_half_y, mod_half_z);
    auto M2_motherLogic = new G4LogicalVolume(M2_motherSolid, Vacuo, "M2_Mother");
    new G4PVPlacement(0, G4ThreeVector(-cx, 0, +cz), M2_motherLogic, "M2_Mother", WorldLV, false, 0, fCheckOverlaps);


    new G4PVParameterised("Prisms_M2",
                          sliceLogic,
                          M2_motherLogic,
                          kUndefined,  // Not tied to a single axis
                          nTotal,
                          param);


    // M3_Mother volume for stack
    auto M3_motherSolid = new G4Box("M3_Mother", mod_half_x, mod_half_y, mod_half_z);
    auto M3_motherLogic = new G4LogicalVolume(M3_motherSolid, Vacuo, "M3_Mother");
    new G4PVPlacement(0, G4ThreeVector(-cx, 0, -cz), M3_motherLogic, "M3_Mother", WorldLV, false, 0, fCheckOverlaps);


    new G4PVParameterised("Prisms_M3",
                          sliceLogic,
                          M3_motherLogic,
                          kUndefined,  // Not tied to a single axis
                          nTotal,
                          param);




    //sliceLogic->SetUserLimits(new G4UserLimits(5 * mm));

    // Production cuts
    auto cuts = new G4ProductionCuts();
    cuts->SetProductionCut(5 * mm, G4ProductionCuts::GetIndex("gamma"));
    cuts->SetProductionCut(5 * mm, G4ProductionCuts::GetIndex("e-"));
    cuts->SetProductionCut(5 * mm, G4ProductionCuts::GetIndex("e+"));

    for (auto logic : {M0_motherLogic, M1_motherLogic, M2_motherLogic, M3_motherLogic}) {
        auto region = new G4Region("LArRegion_" + logic->GetName());  
        region->AddRootLogicalVolume(logic);
        logic->SetRegion(region);
        region->SetProductionCuts(cuts);  
    }

    G4VisAttributes* grey  = new G4VisAttributes (G4Colour (0.5, 0.5, 0.5, 0.8));
    G4VisAttributes* blue  = new G4VisAttributes (G4Colour (0.5, 0.5, 1.0, 0.8));
    G4VisAttributes* red   = new G4VisAttributes (G4Colour (1.0, 0.5, 0.5, 0.8));
    G4VisAttributes* green = new G4VisAttributes (G4Colour (0.5, 1.0, 0.5, 0.8));
    grey->SetForceSolid (true);
    blue->SetForceSolid (true);
    red->SetForceSolid (true);
    green->SetForceSolid (true);
    M0_motherLogic->SetVisAttributes(grey);
    M1_motherLogic->SetVisAttributes(blue);
    M2_motherLogic->SetVisAttributes(red);
    M3_motherLogic->SetVisAttributes(green);



 

 

    return WorldPV;
        
    }
