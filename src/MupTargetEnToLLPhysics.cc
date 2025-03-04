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
thread_local MupTargetEnToLLProcess *MupTargetEnToLLPhysics::fProcess;

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

  fProcess = new MupTargetEnToLLProcess;
  processManager->AddDiscreteProcess(fProcess);
}

void MupTargetEnToLLPhysics::Configure(const std::vector<G4String> &pointsFile, G4double xssf)
{
  fProcess->Configure(pointsFile, xssf);
}
