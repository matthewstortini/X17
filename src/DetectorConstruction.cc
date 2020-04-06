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
   // Ge bege detector #1
   //
   G4RotationMatrix* myRotation1 = new G4RotationMatrix();
   myRotation1->rotateX(90*degree);
   myRotation1->rotateY(-(360*0/7)*degree);
   G4double x1 = 150*std::cos(M_PI/2-(2*M_PI*0/7));
   G4double y1 = 150*std::sin(M_PI/2-(2*M_PI*0/7));
   G4ThreeVector position1(x1,y1,0);

   G4Tubs* bege1_solid
      = new G4Tubs("bege1",   // its name
            0*cm,             // innerRadius
            3.5*cm,           // outerRadius
            1.5*cm,           // half-height 
            0,                // start angle
            2*M_PI);          // spanning angle

   G4LogicalVolume* bege1_logical
      = new G4LogicalVolume(
            bege1_solid,            // its solid
            germanium,              // its material
            "bege1");               // its name

   new G4PVPlacement(
       myRotation1,                  // rotation
       position1,                    // position
       bege1_logical,               // its logical volume                         
       "bege1",                     // its name
       logicWorld,                  // its mother volume
       false,                       // no boolean operation
       0,                           // copy number
       checkOverlaps);              // checking overlaps 

   G4VisAttributes* bege1VisAtt = new G4VisAttributes(G4Colour(1,1.0,0));
   bege1_logical->SetVisAttributes(bege1VisAtt);

   //
   // Ge bege detector #2
   //
   G4RotationMatrix* myRotation2 = new G4RotationMatrix();
   myRotation2->rotateX(90*degree);
   myRotation2->rotateY(-(360*1/7)*degree);
   G4double x2 = 150*std::cos(M_PI/2-(2*M_PI*1/7));
   G4double y2 = 150*std::sin(M_PI/2-(2*M_PI*1/7));   
   G4ThreeVector position2(x2,y2,0);

   G4Tubs* bege2_solid
      = new G4Tubs("bege2",   // its name
            0*cm,             // innerRadius
            3.5*cm,           // outerRadius
            1.5*cm,           // half-height 
            0,                // start angle
            2*M_PI);          // spanning angle

   G4LogicalVolume* bege2_logical
      = new G4LogicalVolume(
            bege2_solid,            // its solid
            germanium,              // its material
            "bege2");               // its name

   new G4PVPlacement(
       myRotation2,                  // rotation
       position2,                    // set position
       bege2_logical,               // its logical volume                         
       "bege2",                     // its name
       logicWorld,                  // its mother volume
       false,                       // no boolean operation
       0,                           // copy number
       checkOverlaps);              // checking overlaps 

   G4VisAttributes* bege2VisAtt = new G4VisAttributes(G4Colour(1,1.0,0));
   bege2_logical->SetVisAttributes(bege2VisAtt);


   //
   // Ge bege detector #3
   //
   G4RotationMatrix* myRotation3 = new G4RotationMatrix();
   myRotation3->rotateX(90*degree);
   myRotation3->rotateY(-(360*2/7)*degree);
   G4double x3 = 150*std::cos(M_PI/2-(2*M_PI*2/7));
   G4double y3 = 150*std::sin(M_PI/2-(2*M_PI*2/7));   
   G4ThreeVector position3(x3,y3,0);

   G4Tubs* bege3_solid
      = new G4Tubs("bege3",   // its name
            0*cm,             // innerRadius
            3.5*cm,           // outerRadius
            1.5*cm,           // half-height 
            0,                // start angle
            2*M_PI);          // spanning angle

   G4LogicalVolume* bege3_logical
      = new G4LogicalVolume(
            bege3_solid,            // its solid
            germanium,              // its material
            "bege3");               // its name

   new G4PVPlacement(
       myRotation3,                  // rotation
       position3,                    // set position
       bege3_logical,               // its logical volume                         
       "bege3",                     // its name
       logicWorld,                  // its mother volume
       false,                       // no boolean operation
       0,                           // copy number
       checkOverlaps);              // checking overlaps 

   G4VisAttributes* bege3VisAtt = new G4VisAttributes(G4Colour(1,1.0,0));
   bege3_logical->SetVisAttributes(bege3VisAtt);

   //
   // Ge bege detector #4
   //
   G4RotationMatrix* myRotation4 = new G4RotationMatrix();
   myRotation4->rotateX(90*degree);
   myRotation4->rotateY(-(360*3/7)*degree);
   G4double x4 = 150*std::cos(M_PI/2-(2*M_PI*3/7));
   G4double y4 = 150*std::sin(M_PI/2-(2*M_PI*3/7));   
   G4ThreeVector position4(x4,y4,0);

   G4Tubs* bege4_solid
      = new G4Tubs("bege4",   // its name
            0*cm,             // innerRadius
            3.5*cm,           // outerRadius
            1.5*cm,           // half-height 
            0,                // start angle
            2*M_PI);          // spanning angle

   G4LogicalVolume* bege4_logical
      = new G4LogicalVolume(
            bege4_solid,            // its solid
            germanium,              // its material
            "bege4");               // its name

   new G4PVPlacement(
       myRotation4,                  // rotation
       position4,                    // set position
       bege4_logical,               // its logical volume                         
       "bege4",                     // its name
       logicWorld,                  // its mother volume
       false,                       // no boolean operation
       0,                           // copy number
       checkOverlaps);              // checking overlaps 

   G4VisAttributes* bege4VisAtt = new G4VisAttributes(G4Colour(1,1.0,0));
   bege4_logical->SetVisAttributes(bege4VisAtt);

   //
   // Ge bege detector #5
   //
   G4RotationMatrix* myRotation5 = new G4RotationMatrix();
   myRotation5->rotateX(90*degree);
   myRotation5->rotateY(-(360*4/7)*degree);
   G4double x5 = 150*std::cos(M_PI/2-(2*M_PI*4/7));
   G4double y5 = 150*std::sin(M_PI/2-(2*M_PI*4/7));
   G4ThreeVector position5(x5,y5,0);

   G4Tubs* bege5_solid
      = new G4Tubs("bege5",   // its name
            0*cm,             // innerRadius
            3.5*cm,           // outerRadius
            1.5*cm,           // half-height 
            0,                // start angle
            2*M_PI);          // spanning angle

   G4LogicalVolume* bege5_logical
      = new G4LogicalVolume(
            bege5_solid,            // its solid
            germanium,              // its material
            "bege4");               // its name

   new G4PVPlacement(
       myRotation5,                  // rotation
       position5,                    // set position
       bege5_logical,               // its logical volume                         
       "bege5",                     // its name
       logicWorld,                  // its mother volume
       false,                       // no boolean operation
       0,                           // copy number
       checkOverlaps);              // checking overlaps 

   G4VisAttributes* bege5VisAtt = new G4VisAttributes(G4Colour(1,1.0,0));
   bege5_logical->SetVisAttributes(bege5VisAtt);

   //
   // Ge bege detector #6
   //
   G4RotationMatrix* myRotation6 = new G4RotationMatrix();
   myRotation6->rotateX(90*degree);
   myRotation6->rotateY(-(360*5/7)*degree);
   G4double x6 = 150*std::cos(M_PI/2-(2*M_PI*5/7));
   G4double y6 = 150*std::sin(M_PI/2-(2*M_PI*5/7));
   G4ThreeVector position6(x6,y6,0);

   G4Tubs* bege6_solid
      = new G4Tubs("bege6",   // its name
            0*cm,             // innerRadius
            3.5*cm,           // outerRadius
            1.5*cm,           // half-height 
            0,                // start angle
            2*M_PI);          // spanning angle

   G4LogicalVolume* bege6_logical
      = new G4LogicalVolume(
            bege6_solid,            // its solid
            germanium,              // its material
            "bege6");               // its name

   new G4PVPlacement(
       myRotation6,                  // rotation
       position6,                    // set position
       bege5_logical,               // its logical volume                         
       "bege6",                     // its name
       logicWorld,                  // its mother volume
       false,                       // no boolean operation
       0,                           // copy number
       checkOverlaps);              // checking overlaps 

   G4VisAttributes* bege6VisAtt = new G4VisAttributes(G4Colour(1,1.0,0));
   bege6_logical->SetVisAttributes(bege6VisAtt);

   //
   // Ge bege detector #7
   //
   G4RotationMatrix* myRotation7 = new G4RotationMatrix();
   myRotation7->rotateX(90*degree);
   myRotation7->rotateY(-(360*6/7)*degree);
   G4double x7 = 150*std::cos(M_PI/2-(2*M_PI*6/7));
   G4double y7 = 150*std::sin(M_PI/2-(2*M_PI*6/7));
   G4ThreeVector position7(x7,y7,0);

   G4Tubs* bege7_solid
      = new G4Tubs("bege7",   // its name
            0*cm,             // innerRadius
            3.5*cm,           // outerRadius
            1.5*cm,           // half-height 
            0,                // start angle
            2*M_PI);          // spanning angle

   G4LogicalVolume* bege7_logical
      = new G4LogicalVolume(
            bege7_solid,            // its solid
            germanium,              // its material
            "bege7");               // its name

   new G4PVPlacement(
       myRotation7,                  // rotation
       position7,                    // set position
       bege7_logical,               // its logical volume                         
       "bege7",                     // its name
       logicWorld,                  // its mother volume
       false,                       // no boolean operation
       0,                           // copy number
       checkOverlaps);              // checking overlaps 

   G4VisAttributes* bege7VisAtt = new G4VisAttributes(G4Colour(1,1.0,0));
   bege7_logical->SetVisAttributes(bege7VisAtt);

   //
   // lithium target foil
   //
   G4Tubs* foil_solid
      = new G4Tubs("foil",   // its name
            0*cm,            // innerRadius
            3.5*cm,          // outerRadius
            3.0045*mm,       // half-height 
            0,               // start angle
            2*M_PI);         // spanning angle

   G4LogicalVolume* foil_logical   
      = new G4LogicalVolume(
            foil_solid,            // its solid
            lithium,               // its material
            "foil");               // its name

   new G4PVPlacement(
       0,                      // no rotation
       G4ThreeVector(0,0,0),   // at (0,0,0)
       foil_logical,           // its logical volume                         
       "foil",                 // its name
       logicWorld,             // its mother volume
       false,                  // no boolean operation
       0,                      // copy number
       checkOverlaps);         // checking overlaps 

   G4VisAttributes* foil_VisAtt = new G4VisAttributes(G4Colour(0,1.0,0));
   foil_logical->SetVisAttributes(foil_VisAtt);

return physWorld;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4RunManager.hh"

void DetectorConstruction::UpdateGeometry() {
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

