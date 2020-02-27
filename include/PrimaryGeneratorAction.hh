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

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "globals.hh"
#include <vector>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ParticleGun;
class G4Event;
class DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// The primary generator action class with particle gum.
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction, public G4UImessenger {

  public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();

    // static access method
    static const PrimaryGeneratorAction* Instance();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);

    //G4GeneralParticleSource* fParticleGun;

    // method to access particle gun
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
    G4ParticleGun*  fParticleGun; // pointer a to G4 gun class
    double eventtype;

    // for setting macro commands (i.e. decay mode, boson mass, resonance energy)
    void SetNewValue(G4UIcommand *command, G4String newValues);

  private:
    static PrimaryGeneratorAction* fgInstance;

    G4double X17mass;
    G4double resonanceenergy;

    TGenPhaseSpace PhaseSpaceEvent;
    TLorentzVector X17;
    TLorentzVector direction;
    TLorentzVector combined;

    std::vector<double> particle_mass;
    std::vector<TLorentzVector*> particle_lorentz;
    std::vector<G4ParticleDefinition*> particle_definition;

    G4UIcmdWithAString* fDecayModeCmd;
    enum EMode { kGamma, kX17, kGun };
    EMode fMode;

    G4UIcmdWithADouble* fX17MassCmd;
    G4UIcmdWithADouble* fResonanceEnergyCmd;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
