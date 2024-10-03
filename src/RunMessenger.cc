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

#include "G4RunManager.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "MupTargetEnToLLPhysics.hh"
#include "PrimaryGeneratorAction.hh"
#include "SteppingAction.hh"

class RunMessenger::Driver {
public:
  Driver(RunMessenger *messenger);
  ~Driver();
  void SetNewValue(G4UIcommand *, G4String);

private:
  PrimaryGeneratorAction *fPrimaryGeneratorAction;
  SteppingAction *fSteppingAction;

  G4UIdirectory *fScatterDir;
  G4UIcommand *fSetMupTargetEnToLLCmd;
  G4UIcmdWithADoubleAndUnit *fSetTotalEnergyCmd;
};

RunMessenger::RunMessenger() { fDriver = new Driver(this); }

RunMessenger::~RunMessenger() { delete fDriver; }

void RunMessenger::SetNewValue(G4UIcommand *cmd, G4String val) { fDriver->SetNewValue(cmd, val); }

RunMessenger::Driver::Driver(RunMessenger *messenger)
{
  fPrimaryGeneratorAction = (PrimaryGeneratorAction *)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  fSteppingAction = (SteppingAction *)G4RunManager::GetRunManager()->GetUserSteppingAction();

  fScatterDir = new G4UIdirectory("/scatter/");
  fScatterDir->SetGuidance("Control scattering processes.");

  fSetMupTargetEnToLLCmd = new G4UIcommand("/scatter/mupTargetEnToLL", messenger);
  fSetMupTargetEnToLLCmd->SetParameter(new G4UIparameter("pid", 'i', false));
  fSetMupTargetEnToLLCmd->SetParameter(new G4UIparameter("points_file", 's', false));
  fSetMupTargetEnToLLCmd->SetParameter(new G4UIparameter("xssf", 'd', false));
  fSetMupTargetEnToLLCmd->SetGuidance("Configure MupTargetEnToLL process.");
  fSetMupTargetEnToLLCmd->AvailableForStates(G4State_Idle);

  fSetTotalEnergyCmd = new G4UIcmdWithADoubleAndUnit("/gun/totalEnergy", messenger);
  fSetTotalEnergyCmd->SetGuidance("Set total energy.");
  fSetTotalEnergyCmd->SetParameterName("TotalEnergy", false);
  fSetTotalEnergyCmd->SetUnitCategory("Energy");
  fSetTotalEnergyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RunMessenger::Driver::~Driver()
{
  delete fSetTotalEnergyCmd;
  delete fSetMupTargetEnToLLCmd;
  delete fScatterDir;
}

void RunMessenger::Driver::SetNewValue(G4UIcommand *cmd, G4String val)
{
  if(cmd == fSetTotalEnergyCmd) {
    fPrimaryGeneratorAction->SetTotalEnergy(fSetTotalEnergyCmd->GetNewDoubleValue(val));
  } else if(cmd == fSetMupTargetEnToLLCmd) {
    G4Tokenizer next(val);
    G4String pid_s = next();
    G4String points_file = next();
    G4String xssf_s = next();
    G4String trailing = next();
    if(pid_s.empty() || points_file.empty() || xssf_s.empty() || !trailing.empty()) {
      throw std::runtime_error("expect 3 arguments");
    }
    G4int pid = stoi(pid_s);
    G4double xssf = stod(xssf_s);
    MupTargetEnToLLPhysics::GetInstance()->Configure(pid, points_file, xssf);
  }
}
