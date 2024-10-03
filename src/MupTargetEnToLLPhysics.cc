#include "MupTargetEnToLLPhysics.hh"

#include <G4Electron.hh>
#include <G4MuonMinus.hh>
#include <G4MuonPlus.hh>
#include <G4Positron.hh>
#include <G4ProcessManager.hh>
#include <G4TauMinus.hh>
#include <G4TauPlus.hh>

#include "MupTargetEnToLLProcess.hh"

MupTargetEnToLLPhysics *MupTargetEnToLLPhysics::fInstance = new MupTargetEnToLLPhysics;
thread_local std::map<G4int, MupTargetEnToLLProcess *> MupTargetEnToLLPhysics::fProcesses;

MupTargetEnToLLPhysics::MupTargetEnToLLPhysics() { }

MupTargetEnToLLPhysics::~MupTargetEnToLLPhysics() { }

void MupTargetEnToLLPhysics::ConstructParticle()
{
  G4Electron::Definition();
  G4Positron::Definition();
  G4MuonPlus::Definition();
  G4MuonMinus::Definition();
  G4TauPlus::Definition();
  G4TauMinus::Definition();
}

void MupTargetEnToLLPhysics::ConstructProcess()
{
  G4ProcessManager *processManager = G4MuonPlus::Definition()->GetProcessManager();

  for(G4int lPid : { 11, 13, 15 }) {
    auto process = new MupTargetEnToLLProcess;
    processManager->AddDiscreteProcess(process);
    fProcesses.emplace(lPid, process);
  }
}

void MupTargetEnToLLPhysics::Configure(G4int lPid, G4String pointsFile, G4double xssf)
{
  MupTargetEnToLLProcess *process = fProcesses.at(lPid);
  process->Configure(lPid, pointsFile, xssf);
}
