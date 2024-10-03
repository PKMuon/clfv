#pragma once
#include <map>

#include "G4VPhysicsConstructor.hh"

class MupTargetEnToLLProcess;

class MupTargetEnToLLPhysics : public G4VPhysicsConstructor {
public:
  void ConstructParticle() override;
  void ConstructProcess() override;

  void Configure(G4int lPid, G4String pointsFile, G4double xssf = 1.0);
  static MupTargetEnToLLPhysics *GetInstance() { return fInstance; }

private:
  MupTargetEnToLLPhysics();
  ~MupTargetEnToLLPhysics() override;

  static MupTargetEnToLLPhysics *fInstance;
  static thread_local std::map<G4int, MupTargetEnToLLProcess *> fProcesses;
};
