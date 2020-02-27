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

#include <regex>

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4UserSteppingAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4PhysicalVolumeStore.hh"

#include "g4root.hh"
#include "g4xml.hh"
#include "g4csv.hh"
#include "g4hdf5.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SteppingAction : public G4UserSteppingAction, public G4UImessenger {

  protected:
    G4UIcommand* fVolIDCmd;
    G4UIcmdWithAString* fOutputFormatCmd;
    G4UIcmdWithAString* fOutputOptionCmd;
    G4UIcmdWithABool* fRecordAllStepsCmd;
    G4UIcmdWithAString* fListVolsCmd;

    enum EFormat { kCsv, kXml, kRoot, kHdf5 };
    EFormat fFormat;
    enum EOption { kStepWise, kEventWise };
    EOption fOption;
    bool fRecordAllSteps;

    vector< pair<regex,string> > fPatternPairs;

    G4int fNEvents;
    G4int fEventNumber;
    vector<G4int> fPID;
    vector<G4int> fTrackID;
    vector<G4int> fParentID;
    vector<G4int> fStepNumber;
    vector<G4double> fKE;
    vector<G4double> fEDep;
    vector<G4double> fX;
    vector<G4double> fY;
    vector<G4double> fZ;
    vector<G4double> fLX;
    vector<G4double> fLY;
    vector<G4double> fLZ;
    vector<G4double> fPdX;
    vector<G4double> fPdY;
    vector<G4double> fPdZ;
    vector<G4double> fT;
    vector<G4int> fVolID;
    vector<G4int> fIRep;

    map<G4VPhysicalVolume*, int> fVolIDMap;

  public:
    SteppingAction();
    G4VAnalysisManager* GetAnalysisManager();
    ~SteppingAction();

    void SetNewValue(G4UIcommand *command, G4String newValues);

    void ResetVars();

    void WriteRow(G4VAnalysisManager* man);

    void UserSteppingAction(const G4Step *step);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
