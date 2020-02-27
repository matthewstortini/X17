//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Polyhedra.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4VSolid.hh"

#include "G4Torus.hh"
#include "G4UnionSolid.hh"

#include "G4String.hh"
#include "math.h"

#include "G4VisAttributes.hh"
#include "G4Color.hh"


#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct() {

   // Get nist material manager
   G4NistManager* nist = G4NistManager::Instance();
   // Option to switch on/off checking of volumes overlaps
   G4bool checkOverlaps = true;

   G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
   G4Material* germanium = nist->FindOrBuildMaterial("G4_Ge");
   G4Material* lithium = nist->FindOrBuildMaterial("G4_Li");

   //
   // World
   //
   G4Box* solidWorld
      = new G4Box("World",    // its name
            3*m, 3*m, 3*m);   // its size

   G4LogicalVolume* logicWorld
      = new G4LogicalVolume(
            solidWorld,          // its solid
            vacuum,              // its material
            "World");            // its name

   G4VPhysicalVolume* physWorld
      = new G4PVPlacement(
            0,                    // no rotation
            G4ThreeVector(),      // at (0,0,0)
            logicWorld,           // its logical volume                         
            "World",              // its name
            0,                    // its mother volume (0 since it is world)
            false,                // no boolean operation
            0,                    // copy number
            checkOverlaps);       // checking overlaps 

   //
   // Ge bege detector
   //
   G4Tubs* bege_solid
      = new G4Tubs("bege",   // its name
            0*cm,            // innerRadius
            3.5*cm,          // outerRadius
            1.5*cm,          // half-height 
            0,               // start angle
            2*M_PI);         // spanning angle

   G4LogicalVolume* bege_logical
      = new G4LogicalVolume(
            bege_solid,            // its solid
            germanium,             // its material
            "bege");               // its name

   G4VPhysicalVolume* bege_phys
      = new G4PVPlacement(
            0,                        // no rotation
            G4ThreeVector(0,0,1*m),   // at (0,0,1000) mm
            bege_logical,             // its logical volume                         
            "bege",                   // its name
            logicWorld,               // its mother volume
            false,                    // no boolean operation
            0,                        // copy number
            checkOverlaps);           // checking overlaps 

   G4VisAttributes* begeVisAtt = new G4VisAttributes(G4Colour(1,1.0,0));
   bege_logical->SetVisAttributes(begeVisAtt);

   //
   // Lithium target foil
   //
   G4Tubs* foil_solid
      = new G4Tubs("foil",   // its name
            0*cm,            // innerRadius
            3.5*cm,          // outerRadius
            1.0*mm,          // half-height 
            0,               // start angle
            2*M_PI);         // spanning angle

   G4LogicalVolume* foil_logical
      = new G4LogicalVolume(
            foil_solid,            // its solid
            lithium,               // its material
            "foil");               // its name

   G4VPhysicalVolume* foil_phys
      = new G4PVPlacement(
            0,                      // no rotation
            G4ThreeVector(0,0,0),   // at (0,0,0) mm
            foil_logical,           // its logical volume                         
            "foil",                 // its name
            logicWorld,             // its mother volume
            false,                  // no boolean operation
            0,                      // copy number
            checkOverlaps);         // checking overlaps 

   G4VisAttributes* foilVisAtt = new G4VisAttributes(G4Colour(0,1.0,0));
   foil_logical->SetVisAttributes(foilVisAtt);

return physWorld;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4RunManager.hh"

void DetectorConstruction::UpdateGeometry() {
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

