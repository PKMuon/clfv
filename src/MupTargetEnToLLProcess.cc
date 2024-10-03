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
  //if(GetLPid()) G4cout << __FUNCTION__ << "(" << GetLPid() << ", " << track.GetTrackID() << ")" << G4endl;
  G4double stepLength = INFINITY;
  if(!fMupTargetEnToLL) {
    *condition = InActivated;
  } else {
    stepLength = GetMeanFreePath(track, previousStepSize, condition) * 0.001;
  }
  //if(GetLPid()) {
  //  G4cout << __FUNCTION__ << "(" << GetLPid() << ", " << track.GetTrackID() << ") -> " << stepLength << G4endl;
  //}
  return stepLength;
}

G4VParticleChange *MupTargetEnToLLProcess::PostStepDoIt(const G4Track &track, const G4Step &step)
{
  //if(GetLPid()) G4cout << __FUNCTION__ << "(" << GetLPid() << ", " << track.GetTrackID() << ")" << G4endl;
  //bool changed = false;
  G4ParticleChange *change = new G4ParticleChange;
  change->Initialize(track);
  do {
    if(!fMupTargetEnToLL) break;
    if(track.GetParticleDefinition()->GetPDGEncoding() != -13) break;
    G4ThreeVector lpMomentum = track.GetMomentum(), lnMomentum;
    G4double xs = fMupTargetEnToLL->Scatter(lpMomentum, lnMomentum);
    if(!(xs > 0)) break;

    xs *= track.GetMaterial()->GetElectronDensity() * fXSSF;
    G4double stepLength = step.GetStepLength();
    G4double logProbKeep = -xs * stepLength;
    G4double logRandom = log(G4UniformRand());
    //G4cout << __FUNCTION__ << "(" << GetLPid() << ", " << track.GetTrackID() << "): stepLength=" << stepLength
    //       << " logProbKeep=" << std::fixed << std::setprecision(8) << logProbKeep << " logRandom=" << logRandom
    //       << std::defaultfloat << G4endl;
    if(logRandom < logProbKeep) break;

    //changed = true;
    auto lp = new G4DynamicParticle(fLpDefinition, lpMomentum);
    auto ln = new G4DynamicParticle(fLnDefinition, lnMomentum);
    fRun->AddScatter(&track, lp, ln);

    change->ProposeEnergy(0);
    change->ProposeVelocity(0);
    change->ProposeTrackStatus(fStopAndKill);
    change->AddSecondary(lp, track.GetGlobalTime(), true);
    change->AddSecondary(ln, track.GetGlobalTime(), true);
  } while(0);
  //if(GetLPid()) {
  //  G4cout << __FUNCTION__ << "(" << GetLPid() << ", " << track.GetTrackID() << ") -> " << std::boolalpha << changed
  //         << std::noboolalpha << G4endl;
  //}
  return change;
}

G4double MupTargetEnToLLProcess::GetCrossSection(const G4double energy, const G4MaterialCutsCouple *couple)
{
  //if(GetLPid()) G4cout << __FUNCTION__ << "(" << GetLPid() << ", " << energy << ")" << G4endl;
  G4double xs = 0.0;
  if(fMupTargetEnToLL) {
    const G4Material *material = couple->GetMaterial();
    G4double electronDensity = material->GetElectronDensity();
    xs = fMupTargetEnToLL->CrossSection(energy + fMuonMass);
    xs *= electronDensity * fXSSF;
  }
  //if(GetLPid()) G4cout << __FUNCTION__ << "(" << GetLPid() << ", " << energy << ") -> " << xs << G4endl;
  return xs;
}

G4double MupTargetEnToLLProcess::MinPrimaryEnergy(const G4ParticleDefinition *definition, const G4Material *)
{
  //if(GetLPid()) G4cout << __FUNCTION__ << "(" << GetLPid() << ", " << definition->GetPDGEncoding() << ")" << G4endl;
  G4double minPrimaryEnergy = INFINITY;
  if(fMupTargetEnToLL && definition->GetPDGEncoding() == -13) {
    minPrimaryEnergy = fMupTargetEnToLL->MinPrimaryEnergy() - definition->GetPDGMass();
  }
  //if(GetLPid()) {
  //  G4cout << __FUNCTION__ << "(" << GetLPid() << ", " << definition->GetPDGEncoding() << ") -> " << minPrimaryEnergy
  //         << G4endl;
  //}
  return minPrimaryEnergy;
}

void MupTargetEnToLLProcess::Configure(G4int lPid, G4String pointsFile, G4double xssf)
{
  fRun = ((RunAction *)G4RunManager::GetRunManager()->GetUserRunAction())->GetRun();
  fMupTargetEnToLL = new MupTargetEnToLL(lPid, pointsFile);
  fLpDefinition = G4ParticleTable::GetParticleTable()->FindParticle(-abs(lPid));
  fLnDefinition = G4ParticleTable::GetParticleTable()->FindParticle(+abs(lPid));
  fXSSF = xssf;
}

G4int MupTargetEnToLLProcess::GetLPid() const { return fLnDefinition ? fLnDefinition->GetPDGEncoding() : 0; }

G4double MupTargetEnToLLProcess::GetMeanFreePath(
    const G4Track &track, [[maybe_unused]] G4double previousStepSize, G4ForceCondition *condition)
{
  //if(GetLPid()) G4cout << __FUNCTION__ << "(" << GetLPid() << ", " << track.GetTrackID() << ")" << G4endl;
  G4double mfp = INFINITY;
  if(track.GetParticleDefinition()->GetPDGEncoding() == -13) {
    G4double xs = GetCrossSection(track.GetKineticEnergy(), track.GetMaterialCutsCouple());
    *condition = xs > 0 ? Forced : NotForced;
    mfp = 1 / xs;
  }
  //if(GetLPid()) G4cout << __FUNCTION__ << "(" << GetLPid() << ", " << track.GetTrackID() << ") -> " << mfp << G4endl;
  return mfp;
}
