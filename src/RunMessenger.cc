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

#include "RunMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "SteppingAction.hh"

#include "G4RunManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4Tokenizer.hh"

class RunMessenger::Driver {
public:
  Driver(RunMessenger *messenger);
  ~Driver();
  void SetNewValue(G4UIcommand *, G4String);

private:
  PrimaryGeneratorAction *fPrimaryGeneratorAction;
  SteppingAction *fSteppingAction;
  G4UIdirectory *fScatterDir;
  G4UIcommand *fSetMupTargetEnToEECmd;
  G4UIcommand *fSetMupTargetEnToMuMuCmd;
  G4UIcmdWithADoubleAndUnit *fSetTotalEnergyCmd;
};

RunMessenger::RunMessenger()
{
  fDriver = new Driver(this);
}

RunMessenger::~RunMessenger()
{
  delete fDriver;
}

void RunMessenger::SetNewValue(G4UIcommand *cmd, G4String val)
{
  fDriver->SetNewValue(cmd, val);
}

RunMessenger::Driver::Driver(RunMessenger *messenger)
{
  fPrimaryGeneratorAction = (PrimaryGeneratorAction *)
    G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  fSteppingAction = (SteppingAction *)
    G4RunManager::GetRunManager()->GetUserSteppingAction();

  fScatterDir = new G4UIdirectory("/scatter");
  fScatterDir->SetGuidance("Control scattering processes.");

  fSetMupTargetEnToEECmd = new G4UIcommand("/scatter/mupTargetEnToEE", messenger);
  fSetMupTargetEnToEECmd->SetParameter(new G4UIparameter("probability", 'd', false));
  fSetMupTargetEnToEECmd->GetParameter(0)->SetParameterRange("probability >= 0.0 && probability <= 1.0");
  fSetMupTargetEnToEECmd->SetParameter(new G4UIparameter("points_file", 's', false));
  fSetMupTargetEnToEECmd->SetGuidance("Configure MupTargetEnToEE process.");
  fSetMupTargetEnToEECmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fSetMupTargetEnToMuMuCmd = new G4UIcommand("/scatter/mupTargetEnToMuMu", messenger);
  fSetMupTargetEnToMuMuCmd->SetParameter(new G4UIparameter("probability", 'd', false));
  fSetMupTargetEnToMuMuCmd->GetParameter(0)->SetParameterRange("probability >= 0.0 && probability <= 1.0");
  fSetMupTargetEnToMuMuCmd->SetParameter(new G4UIparameter("points_file", 's', false));
  fSetMupTargetEnToMuMuCmd->SetGuidance("Configure MupTargetEnToMuMu process.");
  fSetMupTargetEnToMuMuCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fSetTotalEnergyCmd = new G4UIcmdWithADoubleAndUnit("/gun/totalEnergy", messenger);
  fSetTotalEnergyCmd->SetGuidance("Set total energy.");
  fSetTotalEnergyCmd->SetParameterName("TotalEnergy", false);
  fSetTotalEnergyCmd->SetUnitCategory("Energy");
  fSetTotalEnergyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RunMessenger::Driver::~Driver()
{
  delete fSetTotalEnergyCmd;
  delete fSetMupTargetEnToMuMuCmd;
  delete fSetMupTargetEnToEECmd;
  delete fScatterDir;
}

void RunMessenger::Driver::SetNewValue(G4UIcommand *cmd, G4String val)
{
  if(cmd == fSetTotalEnergyCmd) {
    fPrimaryGeneratorAction->SetTotalEnergy(fSetTotalEnergyCmd->GetNewDoubleValue(val));
  } else if(cmd == fSetMupTargetEnToEECmd || cmd == fSetMupTargetEnToMuMuCmd) {
    G4Tokenizer next(val);
    G4String probability_s = next();
    G4String points_file = next();
    G4String trailing = next();
    if(probability_s.empty() || points_file.empty() || !trailing.empty()) throw std::runtime_error("expect 2 arguments");
    G4double probability = std::stod(probability_s);
    if(cmd == fSetMupTargetEnToEECmd) {
      fSteppingAction->SetMupTargetEnToEE(probability, points_file);
    } else {
      fSteppingAction->SetMupTargetEnToMuMu(probability, points_file);
    }
  }
}
