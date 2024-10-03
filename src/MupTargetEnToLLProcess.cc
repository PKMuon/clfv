#include "MupTargetEnToLLProcess.hh"

#include <math.h>

#include <G4DynamicParticle.hh>
#include <G4Material.hh>
#include <G4MaterialCutsCouple.hh>
#include <G4ParticleChange.hh>
#include <G4ParticleTable.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>

#include "MupTargetEnToLL.hh"
#include "Run.hh"
#include "RunAction.hh"

MupTargetEnToLLProcess::MupTargetEnToLLProcess()
    : G4VDiscreteProcess("MupTargetEnToLLProcess", fUserDefined),
      fMupTargetEnToLL(nullptr),
      fLpDefinition(nullptr),
      fLnDefinition(nullptr),
      fXSSF(1.0),
      fMuonMass(G4ParticleTable::GetParticleTable()->GetParticle(13)->GetPDGMass())
{
  // empty
}

MupTargetEnToLLProcess::~MupTargetEnToLLProcess() { delete fMupTargetEnToLL; }

G4double MupTargetEnToLLProcess::PostStepGetPhysicalInteractionLength(
    const G4Track &track, G4double previousStepSize, G4ForceCondition *condition)
{
  if(!fMupTargetEnToLL) {
    *condition = InActivated;
    return INFINITY;
  }
  return GetMeanFreePath(track, previousStepSize, condition) * 0.001;
}

G4VParticleChange *MupTargetEnToLLProcess::PostStepDoIt(const G4Track &track, [[maybe_unused]] const G4Step &step)
{
  if(!fMupTargetEnToLL) return nullptr;
  if(track.GetParticleDefinition()->GetPDGEncoding() != -13) return nullptr;
  G4ThreeVector lpMomentum = track.GetMomentum(), lnMomentum;
  G4double xs = fMupTargetEnToLL->Scatter(lpMomentum, lnMomentum);
  if(!(xs > 0)) return nullptr;

  auto lp = new G4DynamicParticle(fLpDefinition, lpMomentum);
  auto ln = new G4DynamicParticle(fLnDefinition, lnMomentum);
  fRun->AddScatter(&track, lp, ln);

  auto change = new G4ParticleChange;
  change->Initialize(track);
  change->ProposeEnergy(0);
  change->ProposeVelocity(0);
  change->ProposeTrackStatus(fStopAndKill);
  change->AddSecondary(lp, true);
  change->AddSecondary(ln, true);
  return change;
}

G4double MupTargetEnToLLProcess::GetCrossSection(const G4double energy, const G4MaterialCutsCouple *couple)
{
  if(!fMupTargetEnToLL) return 0.0;
  const G4Material *material = couple->GetMaterial();
  G4double electronDensity = material->GetElectronDensity();
  G4double xs = fMupTargetEnToLL->CrossSection(energy + fMuonMass);
  return xs * electronDensity * fXSSF;
}

G4double MupTargetEnToLLProcess::MinPrimaryEnergy(const G4ParticleDefinition *definition, const G4Material *)
{
  if(!fMupTargetEnToLL) return INFINITY;
  if(definition->GetPDGEncoding() != -13) return INFINITY;
  return fMupTargetEnToLL->MinPrimaryEnergy() - definition->GetPDGMass();
}

void MupTargetEnToLLProcess::Configure(G4int lPid, G4String pointsFile, G4double xssf)
{
  fRun = ((RunAction *)G4RunManager::GetRunManager()->GetUserRunAction())->GetRun();
  fMupTargetEnToLL = new MupTargetEnToLL(lPid, pointsFile);
  fLpDefinition = G4ParticleTable::GetParticleTable()->FindParticle(-abs(lPid));
  fLnDefinition = G4ParticleTable::GetParticleTable()->FindParticle(+abs(lPid));
  fXSSF = xssf;
}

G4double MupTargetEnToLLProcess::GetMeanFreePath(
    const G4Track &track, [[maybe_unused]] G4double previousStepSize, G4ForceCondition *condition)
{
  if(track.GetParticleDefinition()->GetPDGEncoding() != -13) return INFINITY;
  G4double xs = GetCrossSection(track.GetKineticEnergy(), track.GetMaterialCutsCouple());
  *condition = xs > 0 ? Forced : NotForced;
  return 1 / xs;
}
