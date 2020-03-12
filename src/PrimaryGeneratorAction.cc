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

#include <iostream>
using namespace std;

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
 
   fX17MassCmd = new G4UIcmdWithADouble("/generator/setX17Mass", this);
   fX17MassCmd->SetGuidance("Set the mass of the hypothetical Boson in units of GeV");
   X17mass = 0.01670;

   fResonanceEnergyCmd = new G4UIcmdWithADouble("/generator/setResonanceEnergy", this);
   fResonanceEnergyCmd->SetGuidance("Set resonance energy in units of GeV");
   resonanceenergy = 0.01815;
 
   fDecayModeCmd = new G4UIcmdWithAString("/generator/setDecayMode", this);
   std::string candidates = "x17 gamma e+e- gun capture";
   fDecayModeCmd->SetCandidates(candidates.c_str());
   fDecayModeCmd->SetGuidance("Set decay mode of interest");
   fMode = "kX17";

   // if not running in gun mode, we load in the positions/momenta of our excited states
   // these positions/momenta are to be obtained through previous runs in capture mode
   if ( fMode != "kGun" && fMode != "kCapture" ) {
      std::ifstream xPositionsFile("capture_xpositions.txt");
      if (xPositionsFile.is_open()) {
         while (xPositionsFile.good()) {
            getline(xPositionsFile, xPositionString);
            double xPositionDouble = atof(xPositionString.c_str());
            xPositionsVec.push_back(xPositionDouble);
         }
         xPositionsVec.pop_back();
      }

      std::ifstream yPositionsFile("capture_ypositions.txt");
      if (yPositionsFile.is_open()) {
         while (yPositionsFile.good()) {
            getline(yPositionsFile, yPositionString);
            double yPositionDouble = atof(yPositionString.c_str());
            yPositionsVec.push_back(yPositionDouble);
         }
         yPositionsVec.pop_back();
      }

      std::ifstream zPositionsFile("capture_zpositions.txt");
      if (zPositionsFile.is_open()) {
         while (zPositionsFile.good()) {
            getline(zPositionsFile, zPositionString);
            double zPositionDouble = atof(zPositionString.c_str());
            zPositionsVec.push_back(zPositionDouble);
         }
         zPositionsVec.pop_back();
      }

      std::ifstream xMomentaFile("capture_xmomenta.txt");
      if (xMomentaFile.is_open()) {
         while (xMomentaFile.good()) {
            getline(xMomentaFile, xMomentaString);
            double xMomentaDouble = atof(xMomentaString.c_str());
            xMomentaDouble /= 1000;   // to convert MeV/c to GeV/c
            xMomentaVec.push_back(xMomentaDouble);
         }
         xMomentaVec.pop_back();
      }

      std::ifstream yMomentaFile("capture_ymomenta.txt");
      if (yMomentaFile.is_open()) {
         while (yMomentaFile.good()) {
            getline(yMomentaFile, yMomentaString);
            double yMomentaDouble = atof(yMomentaString.c_str());
            yMomentaDouble /= 1000;   // to convert MeV/c to GeV/c
            yMomentaVec.push_back(yMomentaDouble);
         }
         yMomentaVec.pop_back();
      }

      std::ifstream zMomentaFile("capture_zmomenta.txt");
      if (zMomentaFile.is_open()) {
         while (zMomentaFile.good()) {
            getline(zMomentaFile, zMomentaString);
            double zMomentaDouble = atof(zMomentaString.c_str());
            zMomentaDouble /= 1000;   // to convert MeV/c to GeV/c
            zMomentaVec.push_back(zMomentaDouble);
         }
         zMomentaVec.pop_back();
      }
   }
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
  
   if ( fMode == "kX17" ) {
      // input necessary masses for this decay mode
      G4double ExcitedStateDecayMasses [2];
      G4double X17DecayMasses [2];
      ExcitedStateDecayMasses[0] = MassBeryllium8;
      ExcitedStateDecayMasses[1] = X17mass;
      X17DecayMasses[0] = G4Electron::ElectronDefinition()->GetPDGMass()/(1.*GeV);
      X17DecayMasses[1] = G4Positron::PositronDefinition()->GetPDGMass()/(1*GeV);

      // random number between 0 and xPositionsVec.size()-1 to choose random capture from capture data loaded in constructor
      positionElement = rand() % xPositionsVec.size();

      // create 4-vector for excited state beryllium based on momentum of captured proton and known energy of excited state
      FourVectorExcitedBeryllium8.SetPx(xMomentaVec[positionElement]);
      FourVectorExcitedBeryllium8.SetPy(yMomentaVec[positionElement]);
      FourVectorExcitedBeryllium8.SetPz(zMomentaVec[positionElement]);
      FourVectorExcitedBeryllium8.SetE(MassBeryllium8+resonanceenergy);

      // generate excited state decay and create pointers to daughter product's 4-vectors
      ExcitedStateDecay.SetDecay(FourVectorExcitedBeryllium8, 2, ExcitedStateDecayMasses);
      ExcitedStateDecay.Generate();
      FourVectorBeryllium8Ptr = ExcitedStateDecay.GetDecay(0);
      FourVectorX17Ptr = ExcitedStateDecay.GetDecay(1);

      // create a 4-vector for the X17 to get ready for next decay
      FourVectorX17.SetPx(FourVectorX17Ptr->Px());
      FourVectorX17.SetPy(FourVectorX17Ptr->Py());
      FourVectorX17.SetPz(FourVectorX17Ptr->Pz());
      FourVectorX17.SetE(FourVectorX17Ptr->E());

      // generate X17 decay and create pointers to the e+/e- pair's 4-vectors
      X17Decay.SetDecay(FourVectorX17, 2, X17DecayMasses);
      X17Decay.Generate();
      FourVectorElectronPtr = X17Decay.GetDecay(0);
      FourVectorPositronPtr = X17Decay.GetDecay(1);

      // create e- and e+ primaries using parameters given from their 4-vectors, and place primaries at position of capture
      G4PrimaryVertex* ElectronVertex = new G4PrimaryVertex(G4ThreeVector(xPositionsVec[positionElement], yPositionsVec[positionElement], zPositionsVec[positionElement]),0.0*s);
      G4PrimaryParticle* ElectronPrimary = new G4PrimaryParticle(G4Electron::ElectronDefinition(), FourVectorElectronPtr->Px()*GeV,FourVectorElectronPtr->Py()*GeV,FourVectorElectronPtr->Pz()*GeV);
      ElectronVertex->SetPrimary(ElectronPrimary);
      anEvent->AddPrimaryVertex(ElectronVertex);
      G4PrimaryVertex* PositronVertex = new G4PrimaryVertex(G4ThreeVector(xPositionsVec[positionElement], yPositionsVec[positionElement], zPositionsVec[positionElement]),0.0*s);
      G4PrimaryParticle* PositronPrimary = new G4PrimaryParticle(G4Positron::PositronDefinition(), FourVectorPositronPtr->Px()*GeV,FourVectorPositronPtr->Py()*GeV,FourVectorPositronPtr->Pz()*GeV);
      PositronVertex->SetPrimary(PositronPrimary);
      anEvent->AddPrimaryVertex(PositronVertex);
   }
   
   if ( fMode == "kGamma" ) {
      // random number between 0 and xPositionsVec.size()-1
      positionElement = rand() % xPositionsVec.size();
      G4PrimaryVertex* vertex = new G4PrimaryVertex(G4ThreeVector(xPositionsVec[positionElement], yPositionsVec[positionElement], zPositionsVec[positionElement]),0.0*s);
      G4PrimaryParticle* thePrimaryParticle = new G4PrimaryParticle(G4Gamma::GammaDefinition(),ux*resonanceenergy*GeV,uy*resonanceenergy*GeV,uz*resonanceenergy*GeV);
      vertex->SetPrimary(thePrimaryParticle);
      anEvent->AddPrimaryVertex(vertex);
   }

   if ( fMode == "kEplusEminus" ) {
      // input necessary masses for this decay mode
      G4double ExcitedStateDecayMasses [3];
      ExcitedStateDecayMasses[0] = MassBeryllium8;
      ExcitedStateDecayMasses[1] = G4Electron::ElectronDefinition()->GetPDGMass()/(1.*GeV);
      ExcitedStateDecayMasses[2] = G4Positron::PositronDefinition()->GetPDGMass()/(1*GeV);

      // random number between 0 and xPositionsVec.size()-1 to choose random capture from capture data loaded in constructor
      positionElement = rand() % xPositionsVec.size();

      // create 4-vector for excited state beryllium based on momentum of captured proton and known energy of excited state
      FourVectorExcitedBeryllium8.SetPx(xMomentaVec[positionElement]);
      FourVectorExcitedBeryllium8.SetPy(yMomentaVec[positionElement]);
      FourVectorExcitedBeryllium8.SetPz(zMomentaVec[positionElement]);
      FourVectorExcitedBeryllium8.SetE(MassBeryllium8+resonanceenergy);

      // generate excited state decay and create pointers to daughter product's 4-vectors
      ExcitedStateDecay.SetDecay(FourVectorExcitedBeryllium8, 3, ExcitedStateDecayMasses);
      ExcitedStateDecay.Generate();
      FourVectorBeryllium8Ptr = ExcitedStateDecay.GetDecay(0);
      FourVectorElectronPtr = ExcitedStateDecay.GetDecay(1);
      FourVectorPositronPtr = ExcitedStateDecay.GetDecay(2);

      // create e- and e+ primaries using parameters given from their 4-vectors, and place primaries at position of capture
      G4PrimaryVertex* ElectronVertex = new G4PrimaryVertex(G4ThreeVector(xPositionsVec[positionElement], yPositionsVec[positionElement], zPositionsVec[positionElement]),0.0*s);
      G4PrimaryParticle* ElectronPrimary = new G4PrimaryParticle(G4Electron::ElectronDefinition(), FourVectorElectronPtr->Px()*GeV,FourVectorElectronPtr->Py()*GeV,FourVectorElectronPtr->Pz()*GeV);
      ElectronVertex->SetPrimary(ElectronPrimary);
      anEvent->AddPrimaryVertex(ElectronVertex);
      G4PrimaryVertex* PositronVertex = new G4PrimaryVertex(G4ThreeVector(xPositionsVec[positionElement], yPositionsVec[positionElement], zPositionsVec[positionElement]),0.0*s);
      G4PrimaryParticle* PositronPrimary = new G4PrimaryParticle(G4Positron::PositronDefinition(), FourVectorPositronPtr->Px()*GeV,FourVectorPositronPtr->Py()*GeV,FourVectorPositronPtr->Pz()*GeV);
      PositronVertex->SetPrimary(PositronPrimary);
      anEvent->AddPrimaryVertex(PositronVertex);
   }
   
   if ( fMode == "kGun" || fMode == "kCapture" ) {;}

   fParticleGun->GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetNewValue(G4UIcommand *command, G4String newValues) {
   
   if(command == fDecayModeCmd) {   
      if(newValues == "x17") fMode = "kX17";
      if(newValues == "gamma") fMode = "kGamma";
      if(newValues == "e+e-") fMode = "kEplusEminus";   
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
