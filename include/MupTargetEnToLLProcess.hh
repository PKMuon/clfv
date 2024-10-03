#pragma once
#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <G4VDiscreteProcess.hh>

class MupTargetEnToLL;
class G4ParticleDefinition;
class Run;

class MupTargetEnToLLProcess : public G4VDiscreteProcess {
public:
  MupTargetEnToLLProcess();
  ~MupTargetEnToLLProcess() override;

  G4double PostStepGetPhysicalInteractionLength(
      const G4Track &track, G4double previousStepSize, G4ForceCondition *condition) override;
  G4VParticleChange *PostStepDoIt(const G4Track &, const G4Step &) override;
  G4double GetCrossSection(const G4double, const G4MaterialCutsCouple *) override;
  G4double MinPrimaryEnergy(const G4ParticleDefinition *, const G4Material *) override;

  void Configure(G4int lPid, G4String pointsFile, G4double xssf = 1.0);

protected:
  G4double GetMeanFreePath(const G4Track &aTrack, G4double previousStepSize, G4ForceCondition *condition) override;

private:
  Run *fRun;
  MupTargetEnToLL *fMupTargetEnToLL;
  G4ParticleDefinition *fLpDefinition;
  G4ParticleDefinition *fLnDefinition;
  G4double fXSSF;  // cross section scale factor
  G4double fMuonMass;
};
