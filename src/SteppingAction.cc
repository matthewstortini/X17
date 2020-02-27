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

#include "SteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::ResetVars() {
   
   fPID.clear();
   fTrackID.clear();
   fParentID.clear();
   fStepNumber.clear();
   fKE.clear();
   fEDep.clear();
   fX.clear();
   fY.clear();
   fZ.clear();
   fLX.clear();
   fLY.clear();
   fLZ.clear();
   fPdX.clear();
   fPdY.clear();
   fPdZ.clear();
   fT.clear();
   fVolID.clear();
   fIRep.clear();
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction() : fNEvents(0), fEventNumber(0) {

   ResetVars();

   fVolIDCmd = new G4UIcommand("/geometry/setVolID", this);
   fVolIDCmd->SetParameter(new G4UIparameter("pattern", 's', false));
   fVolIDCmd->SetParameter(new G4UIparameter("replacement", 's', false));
   fVolIDCmd->SetGuidance("Volumes with name matching [pattern] will be given volume ID "
                          "based on the [replacement] rule. Replacement rule must produce an integer."
                          " Patterns which replace to 0 or -1 are forbidden and will be omitted.");

   fOutputFormatCmd = new G4UIcmdWithAString("/analysis/setOutputFormat", this);
   string candidates = "csv xml root hdf5";
   fOutputFormatCmd->SetCandidates(candidates.c_str());
   fOutputFormatCmd->SetGuidance("Set output format");
   fFormat = kCsv;

   fOutputOptionCmd = new G4UIcmdWithAString("/analysis/setOutputOption", this);
   candidates = "stepwise eventwise";
   fOutputOptionCmd->SetCandidates(candidates.c_str());
   fOutputOptionCmd->SetGuidance("Set output option:");
   fOutputOptionCmd->SetGuidance("  stepwise: one row per step");
   fOutputOptionCmd->SetGuidance("  eventwise: one row per event");
   fOption = kStepWise;

   fRecordAllStepsCmd = new G4UIcmdWithABool("/tracking/recordAllSteps", this);
   fRecordAllStepsCmd->SetParameterName("recordAllSteps", true);
   fRecordAllStepsCmd->SetDefaultValue(true);
   fRecordAllStepsCmd->SetGuidance("Write out every single step, not just those in sensitive volumes.");
   fRecordAllSteps = false;

   fListVolsCmd = new G4UIcmdWithAString("/geometry/listPhysVols", this);
   fListVolsCmd->SetParameterName("pattern", true);
   fListVolsCmd->SetGuidance("List name of all instantiated physical volumes");
   fListVolsCmd->SetGuidance("Optionally supply a regex pattern to only list matching volume names");
   fListVolsCmd->AvailableForStates(G4State_Idle, G4State_GeomClosed, G4State_EventProc);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VAnalysisManager* SteppingAction::GetAnalysisManager() {

   if(fFormat == kCsv) return G4Csv::G4AnalysisManager::Instance();
   if(fFormat == kXml) return G4Xml::G4AnalysisManager::Instance();
   if(fFormat == kRoot) return G4Root::G4AnalysisManager::Instance();
   if(fFormat == kHdf5) return G4Hdf5::G4AnalysisManager::Instance();
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction() {
      
   G4VAnalysisManager* man = GetAnalysisManager();
   
   if(man->IsOpenFile()) {
      if(fOption == kEventWise && fPID.size()>0) WriteRow(man);
      man->Write();
      man->CloseFile();
   }

   delete man;
   delete fVolIDCmd;
   delete fOutputFormatCmd;
   delete fOutputOptionCmd;
   delete fRecordAllStepsCmd;
   delete fListVolsCmd;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::SetNewValue(G4UIcommand *command, G4String newValues) {

   if(command == fVolIDCmd) {
      istringstream iss(newValues);
      string pattern;
      string replacement;
      iss >> pattern >> replacement;
      cout << "in: " << pattern << ' ' << replacement << endl;
      fPatternPairs.push_back(pair<regex,string>(regex(pattern),replacement));
   }
   
   if(command == fOutputFormatCmd) {
      // also set recommended options.
      // override option by subsequent call to /analysis/setOutputOption
      if(newValues == "csv") {
         fFormat = kCsv;
         fOption = kStepWise;
      }
      if(newValues == "xml") {
         fFormat = kXml;
         fOption = kEventWise;
      }
      if(newValues == "root") {
         fFormat = kRoot;
         fOption = kEventWise;
      }
      if(newValues == "hdf5") {
         fFormat = kHdf5;
         fOption = kStepWise;
      }
      GetAnalysisManager(); // call once to make all of the /analysis commands available
   }

   if(command == fOutputOptionCmd) {
      if(newValues == "stepwise") fOption = kStepWise;
      if(newValues == "eventwise") fOption = kEventWise;
   }
      
   if(command == fRecordAllStepsCmd) fRecordAllSteps = fRecordAllStepsCmd->GetNewBoolValue(newValues);
   
   if(command == fListVolsCmd) {
      regex pattern(newValues);
      bool doMatching = (newValues != "");
      G4PhysicalVolumeStore* volumeStore = G4PhysicalVolumeStore::GetInstance();
      cout << "Physical volumes";
      if(doMatching) {cout << " matching pattern " << newValues;}
      cout << ":" << endl;
      for(size_t i=0; i<volumeStore->size(); i++) {
         string name = volumeStore->at(i)->GetName();
         if(!doMatching || regex_match(name, pattern)) {cout << name << endl;}
      }
   }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::WriteRow(G4VAnalysisManager* man) {
      
   man->FillNtupleIColumn(0, fNEvents);
   man->FillNtupleIColumn(1, fEventNumber);
   int row = 2;
   
   if(fOption == kStepWise) {
      size_t i = fPID.size()-1;
      man->FillNtupleIColumn(row++, fPID[i]);
      man->FillNtupleIColumn(row++, fTrackID[i]);
      man->FillNtupleIColumn(row++, fParentID[i]);
      man->FillNtupleIColumn(row++, fStepNumber[i]);
      man->FillNtupleDColumn(row++, fKE[i]);
      man->FillNtupleDColumn(row++, fEDep[i]);
      man->FillNtupleDColumn(row++, fX[i]);
      man->FillNtupleDColumn(row++, fY[i]);
      man->FillNtupleDColumn(row++, fZ[i]);
      man->FillNtupleDColumn(row++, fLX[i]);
      man->FillNtupleDColumn(row++, fLY[i]);
      man->FillNtupleDColumn(row++, fLZ[i]);
      man->FillNtupleDColumn(row++, fPdX[i]);
      man->FillNtupleDColumn(row++, fPdY[i]);
      man->FillNtupleDColumn(row++, fPdZ[i]);
      man->FillNtupleDColumn(row++, fT[i]);
      man->FillNtupleIColumn(row++, fVolID[i]);
      man->FillNtupleIColumn(row++, fIRep[i]);
   }
      
   // for event-wise, manager copies data from vectors over
   // automatically in the next line
   man->AddNtupleRow();
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step *step) {

   G4VAnalysisManager* man = GetAnalysisManager();

   if(!man->IsOpenFile()) {
      // need to create the ntuple before opening the file in order to avoid
      // writing error in csv, xml, and hdf5
      man->CreateNtuple("g4sntuple", "steps data");
      man->CreateNtupleIColumn("nEvents");
      man->CreateNtupleIColumn("event");
      if(fOption == kEventWise) {
         man->CreateNtupleIColumn("pid", fPID);
         man->CreateNtupleIColumn("trackID", fTrackID);
         man->CreateNtupleIColumn("parentID", fParentID);
         man->CreateNtupleIColumn("step", fStepNumber);
         man->CreateNtupleDColumn("KE", fKE);
         man->CreateNtupleDColumn("Edep", fEDep);
         man->CreateNtupleDColumn("x", fX);
         man->CreateNtupleDColumn("y", fY);
         man->CreateNtupleDColumn("z", fZ);
         man->CreateNtupleDColumn("lx", fLX);
         man->CreateNtupleDColumn("ly", fLY);
         man->CreateNtupleDColumn("lz", fLZ);
         man->CreateNtupleDColumn("pdx", fPdX);
         man->CreateNtupleDColumn("pdy", fPdY);
         man->CreateNtupleDColumn("pdz", fPdZ);
         man->CreateNtupleDColumn("t", fT);
         man->CreateNtupleIColumn("volID", fVolID);
         man->CreateNtupleIColumn("iRep", fIRep);
      }
      else if(fOption == kStepWise) {
         man->CreateNtupleIColumn("pid");
         man->CreateNtupleIColumn("trackID");
         man->CreateNtupleIColumn("parentID");
         man->CreateNtupleIColumn("step");
         man->CreateNtupleDColumn("KE");
         man->CreateNtupleDColumn("Edep");
         man->CreateNtupleDColumn("x");
         man->CreateNtupleDColumn("y");
         man->CreateNtupleDColumn("z");
         man->CreateNtupleDColumn("lx");
         man->CreateNtupleDColumn("ly");
         man->CreateNtupleDColumn("lz");
         man->CreateNtupleDColumn("pdx");
         man->CreateNtupleDColumn("pdy");
         man->CreateNtupleDColumn("pdz");
         man->CreateNtupleDColumn("t");
         man->CreateNtupleIColumn("volID");
         man->CreateNtupleIColumn("iRep");
      }
      else {
         cout << "ERROR: Unknown output option " << fOption << endl;
         return;
      }
        
      man->FinishNtuple();

      // look for filename set by macro command: /analysis/setFileName [name]
      if(man->GetFileName() == "") man->SetFileName("g4simpleout");
      cout << "Opening file " << man->GetFileName() << endl;
      man->OpenFile();

      ResetVars();
      fNEvents = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
      fVolIDMap.clear();
   }

   fEventNumber = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
   static G4int lastEventID = fEventNumber;
   
   if(fEventNumber != lastEventID) {
      if(fOption == kEventWise && fPID.size()>0) WriteRow(man);
      ResetVars();
      lastEventID = fEventNumber;
   }
      
   // post-step point will always work: only need to use the pre-step point
   // on the first step, for which the pre-step volume is always the same as
   // the post-step volume
   G4VPhysicalVolume* vpv = step->GetPostStepPoint()->GetPhysicalVolume();
   G4int id = fVolIDMap[vpv];
   
   if(id == 0 && fPatternPairs.size() > 0) {
      string name = (vpv == NULL) ? "NULL" : vpv->GetName();
      for(auto& pp : fPatternPairs) {
         if(regex_match(name, pp.first)) {
            string replaced = regex_replace(name,pp.first,pp.second);
            cout << "match: " << name << ' ' << /*pp.first.str() << ' ' <<*/ pp.second << ' ' << replaced << endl;
            //int id_new = stoi(regex_replace(name,pp.first,pp.second));
            int id_new = stoi(replaced);
            if (id_new == 0 || id_new == -1) {
               cout << "Volume " << name << ": Can't use ID = " << id << endl;
            }
            else {
               id = id_new;
            }
            break;
         }
      }
      if(id == 0 && !fRecordAllSteps) id = -1;
         fVolIDMap[vpv] = id;
   }

   // always record primary event info from pre-step of first step
   // if recording all steps, do this block to record prestep info
   if(fVolID.size() == 0 || (fRecordAllSteps && step->GetTrack()->GetCurrentStepNumber() == 1)) {
      fVolID.push_back(id == -1 ? 0 : id);
      fPID.push_back(step->GetTrack()->GetParticleDefinition()->GetPDGEncoding());
      fTrackID.push_back(step->GetTrack()->GetTrackID());
      fParentID.push_back(step->GetTrack()->GetParentID());
      fStepNumber.push_back(0); // call this step "0"
      fKE.push_back(step->GetPreStepPoint()->GetKineticEnergy());
      fEDep.push_back(0);
      G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();
      fX.push_back(pos.x());
      fY.push_back(pos.y());
      fZ.push_back(pos.z());
      G4TouchableHandle vol = step->GetPreStepPoint()->GetTouchableHandle();
      G4ThreeVector lPos = vol->GetHistory()->GetTopTransform().TransformPoint(pos);
      fLX.push_back(lPos.x());
      fLY.push_back(lPos.y());
      fLZ.push_back(lPos.z());
      G4ThreeVector momDir = step->GetPreStepPoint()->GetMomentumDirection();
      fPdX.push_back(momDir.x());
      fPdY.push_back(momDir.y());
      fPdZ.push_back(momDir.z());
      fT.push_back(step->GetPreStepPoint()->GetGlobalTime());
      fIRep.push_back(vol->GetReplicaNumber());
      if(fOption == kStepWise) WriteRow(man);
   }

   // If not in a sensitive volume, get out of here.
   if(id == -1) return;

   // Don't write Edep=0 steps (unless desired)
   if(!fRecordAllSteps && step->GetTotalEnergyDeposit() == 0) return;

   // Now record post-step info
   fVolID.push_back(id);
   fPID.push_back(step->GetTrack()->GetParticleDefinition()->GetPDGEncoding());
   fTrackID.push_back(step->GetTrack()->GetTrackID());
   fParentID.push_back(step->GetTrack()->GetParentID());
   fStepNumber.push_back(step->GetTrack()->GetCurrentStepNumber());
   fKE.push_back(step->GetTrack()->GetKineticEnergy());
   fEDep.push_back(step->GetTotalEnergyDeposit());
   G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();
   fX.push_back(pos.x());
   fY.push_back(pos.y());
   fZ.push_back(pos.z());
   G4TouchableHandle vol = step->GetPostStepPoint()->GetTouchableHandle();
   G4ThreeVector lPos = vol->GetHistory()->GetTopTransform().TransformPoint(pos);
   fLX.push_back(lPos.x());
   fLY.push_back(lPos.y());
   fLZ.push_back(lPos.z());
   G4ThreeVector momDir = step->GetPostStepPoint()->GetMomentumDirection();
   fPdX.push_back(momDir.x());
   fPdY.push_back(momDir.y());
   fPdZ.push_back(momDir.z());
   fT.push_back(step->GetPostStepPoint()->GetGlobalTime());
   fIRep.push_back(vol->GetReplicaNumber());
   if(fOption == kStepWise) WriteRow(man);
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

