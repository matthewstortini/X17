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

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "TLorentzVector.h"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction* PrimaryGeneratorAction::fgInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// static access function via G4RunManager
const PrimaryGeneratorAction* PrimaryGeneratorAction::Instance() {   
   return fgInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(), fParticleGun(0) {

   fParticleGun  = new G4ParticleGun();
   particle_definition.push_back(G4Electron::ElectronDefinition());
   particle_definition.push_back(G4Positron::PositronDefinition());

   for (int i = 0;i<particle_definition.size();i++) {
      particle_mass.push_back(particle_definition[i]->GetPDGMass()/(1.*GeV));
   }

   particle_lorentz.reserve(particle_definition.size());
 
   fX17MassCmd = new G4UIcmdWithADouble("/generator/setX17Mass", this);
   fX17MassCmd->SetGuidance("Set the mass of the hypothetical Boson in units of GeV");
   X17mass = 0.01670;

   fResonanceEnergyCmd = new G4UIcmdWithADouble("/generator/setResonanceEnergy", this);
   fResonanceEnergyCmd->SetGuidance("Set resonance energy in units of GeV");
   resonanceenergy = 0.01815;
 
   fDecayModeCmd = new G4UIcmdWithAString("/generator/setDecayMode", this);
   std::string candidates = "gamma x17 gun capture";
   fDecayModeCmd->SetCandidates(candidates.c_str());
   fDecayModeCmd->SetGuidance("Set decay mode of interest");
   fMode = "kX17";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooO

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
   
   delete fParticleGun;
   delete fX17MassCmd;
   delete fResonanceEnergyCmd;
   delete fDecayModeCmd;
   fgInstance = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

   G4double cosTheta = 2*G4UniformRand() - 1.;
   G4double phi = 2*CLHEP::pi*G4UniformRand();
   G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
   G4double ux = sinTheta*std::cos(phi),
            uy = sinTheta*std::sin(phi),
            uz = cosTheta;
  
   // if not running in gun mode, we load in the positions of our excited states
   // these positions are to be obtained through previous runs in gun mode
   if ( fMode != "kGun" && fMode != "kCapture" ) {
      std::ifstream xPositionsFile("X17_xpositions.txt");
      if (xPositionsFile.is_open()) {
         while (xPositionsFile.good()) {
            getline(xPositionsFile, xPositionString);
            double xPositionDouble = atof(xPositionString.c_str());
            xPositionsVec.push_back(xPositionDouble);
         }
         xPositionsVec.pop_back();
      }

      std::ifstream yPositionsFile("X17_ypositions.txt");
      if (yPositionsFile.is_open()) {
         while (yPositionsFile.good()) {
            getline(yPositionsFile, yPositionString);
            double yPositionDouble = atof(yPositionString.c_str());
            yPositionsVec.push_back(yPositionDouble);
         }
         yPositionsVec.pop_back();
      }

      std::ifstream zPositionsFile("X17_zpositions.txt");
      if (zPositionsFile.is_open()) {
         while (zPositionsFile.good()) {
            getline(zPositionsFile, zPositionString);
            double zPositionDouble = atof(zPositionString.c_str());
            zPositionsVec.push_back(zPositionDouble);
         }
         zPositionsVec.pop_back();
      }
   }
   
   if ( fMode == "kX17" ) {
      FourVector.SetPx(ux*std::sqrt(resonanceenergy*resonanceenergy-X17mass*X17mass));
      FourVector.SetPy(uy*std::sqrt(resonanceenergy*resonanceenergy-X17mass*X17mass));
      FourVector.SetPz(uz*std::sqrt(resonanceenergy*resonanceenergy-X17mass*X17mass));
      FourVector.SetE(resonanceenergy); //in GeV

      PhaseSpaceEvent.SetDecay(FourVector, particle_definition.size(),particle_mass.data());

      // random number between 0 and xPositionsVec.size()-1
      positionElement = rand() % xPositionsVec.size();

      double Weight = PhaseSpaceEvent.Generate();
      for ( int i = 0;i<particle_definition.size();i++ ) {
         G4PrimaryVertex* vertex = new G4PrimaryVertex(G4ThreeVector(xPositionsVec[positionElement], yPositionsVec[positionElement], zPositionsVec[positionElement]),0.0*s);
         particle_lorentz[i] = PhaseSpaceEvent.GetDecay(i);
         G4PrimaryParticle* thePrimaryParticle = new G4PrimaryParticle(particle_definition[i], particle_lorentz[i]->Px()*GeV,particle_lorentz[i]->Py()*GeV,particle_lorentz[i]->Pz()*GeV);
         vertex->SetPrimary(thePrimaryParticle);
         vertex->SetWeight(Weight);
         anEvent->AddPrimaryVertex(vertex);
      }
   }
   
   if ( fMode == "kGamma" ) {
      // random number between 0 and xPositionsVec.size()-1
      positionElement = rand() % xPositionsVec.size();
      G4PrimaryVertex* vertex = new G4PrimaryVertex(G4ThreeVector(xPositionsVec[positionElement], yPositionsVec[positionElement], zPositionsVec[positionElement]),0.0*s);
      G4PrimaryParticle* thePrimaryParticle = new G4PrimaryParticle(G4Gamma::GammaDefinition(),ux*resonanceenergy*GeV,uy*resonanceenergy*GeV,uz*resonanceenergy*GeV);
      vertex->SetPrimary(thePrimaryParticle);
      anEvent->AddPrimaryVertex(vertex);
   }
   
   if ( fMode == "kGun" || fMode == "kCapture" ) {;}

   fParticleGun->GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetNewValue(G4UIcommand *command, G4String newValues) {
   
   if(command == fDecayModeCmd) {   
      if(newValues == "x17") fMode = "kX17";
      if(newValues == "gamma") fMode = "kGamma";
      if(newValues == "gun") fMode = "kGun";
      if(newValues == "capture") fMode = "kCapture";
   }

   if (command == fX17MassCmd) X17mass = fX17MassCmd->GetNewDoubleValue(newValues);
   if (command == fResonanceEnergyCmd) resonanceenergy = fResonanceEnergyCmd->GetNewDoubleValue(newValues);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::string PrimaryGeneratorAction::GetfMode() {
   return fMode;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
