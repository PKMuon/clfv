#pragma once
#include <vector>

#include "G4VPhysicsConstructor.hh"

class MupTargetEnToLLProcess;

class MupTargetEnToLLPhysics : public G4VPhysicsConstructor {
public:
  void ConstructParticle() override;
  void ConstructProcess() override;

  void Configure(const std::vector<G4String> &rootfiles, G4double xssf = 1.0);
  static MupTargetEnToLLPhysics *GetInstance() { return fInstance; }

private:
  MupTargetEnToLLPhysics();
  ~MupTargetEnToLLPhysics() override;

  static MupTargetEnToLLPhysics *fInstance;
  static thread_local MupTargetEnToLLProcess *fProcess;
};
