#include "MupTargetEnToLLProcess.hh"

#include <math.h>

#include <G4DynamicParticle.hh>
#include <G4Electron.hh>
#include <G4Material.hh>
#include <G4MaterialCutsCouple.hh>
#include <G4MuonPlus.hh>
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
      fRun(nullptr),
      fMupTargetEnToLL(nullptr),
      fMuonPlus(nullptr),
      fElectron(nullptr),
      fXSSF(0.0),
      fMuonMass(0.0)
{
  // empty
}

MupTargetEnToLLProcess::~MupTargetEnToLLProcess() { delete fMupTargetEnToLL; }

G4double MupTargetEnToLLProcess::PostStepGetPhysicalInteractionLength(
    const G4Track &track, G4double previousStepSize, G4ForceCondition *condition)
{
  //G4cout << __FUNCTION__ << "(" << track.GetTrackID() << ")" << G4endl;
  G4double stepLength = GetMeanFreePath(track, previousStepSize, condition) * 0.001;
  //G4cout << __FUNCTION__ << "(" << track.GetTrackID() << ") -> " << stepLength << G4endl;
  return stepLength;
}

G4VParticleChange *MupTargetEnToLLProcess::PostStepDoIt(const G4Track &track, const G4Step &step)
{
  //G4cout << __FUNCTION__ << "(" << track.GetTrackID() << ")" << G4endl;
  //bool changed = false;
  thread_local G4ParticleChange change;
  change.Initialize(track);
  do {
    if(!fMupTargetEnToLL) break;
    if(track.GetParticleDefinition()->GetPDGEncoding() != -13) break;
    G4double xs = GetCrossSection(track.GetKineticEnergy(), track.GetMaterialCutsCouple());
    if(!(xs > 0)) break;

    G4double stepLength = step.GetStepLength();
    G4double logProbKeep = -xs * stepLength;
    G4double logRandom = log(G4UniformRand());
    //G4cout << __FUNCTION__ << "(" << track.GetTrackID() << "): stepLength=" << stepLength
    //       << " logProbKeep=" << std::fixed << std::setprecision(8) << logProbKeep << " logRandom=" << logRandom
    //       << std::defaultfloat << G4endl;
    if(logRandom < logProbKeep) break;

    //changed = true;
    G4ThreeVector lpMomentum = track.GetMomentum(), lnMomentum;
    fMupTargetEnToLL->Scatter(lpMomentum, lnMomentum);
    auto lp = new G4DynamicParticle(fMuonPlus, lpMomentum);
    auto ln = new G4DynamicParticle(fElectron, lnMomentum);
    fRun->AddScatter(&track, lp, ln);

    change.ProposeEnergy(0);
    change.ProposeVelocity(0);
    change.ProposeTrackStatus(fStopAndKill);
    change.AddSecondary(lp, track.GetGlobalTime(), true);
    change.AddSecondary(ln, track.GetGlobalTime(), true);
  } while(0);
  //G4cout << __FUNCTION__ << "(" << track.GetTrackID() << ") -> " << std::boolalpha << changed
  //       << std::noboolalpha << G4endl;
  return &change;
}

G4double MupTargetEnToLLProcess::GetCrossSection(const G4double energy, const G4MaterialCutsCouple *couple)
{
  //G4cout << __FUNCTION__ << "(" << energy << ")" << G4endl;
  G4double xs = 0.0;
  if(fMupTargetEnToLL) {
    const G4Material *material = couple->GetMaterial();
    G4double electronDensity = material->GetElectronDensity();
    xs = fMupTargetEnToLL->CrossSection(energy + fMuonMass);
    xs *= electronDensity * fXSSF;
  }
  //G4cout << __FUNCTION__ << "(" << energy << ") -> " << xs << G4endl;
  return xs;
}

G4double MupTargetEnToLLProcess::MinPrimaryEnergy(const G4ParticleDefinition *definition, const G4Material *)
{
  //G4cout << __FUNCTION__ << "(" << definition->GetPDGEncoding() << ")" << G4endl;
  G4double minPrimaryEnergy = INFINITY;
  if(fMupTargetEnToLL && definition->GetPDGEncoding() == -13) {
    minPrimaryEnergy = fMupTargetEnToLL->MinPrimaryEnergy() - definition->GetPDGMass();
  }
  //G4cout << __FUNCTION__ << "(" << definition->GetPDGEncoding() << ") -> " << minPrimaryEnergy
  //       << G4endl;
  return minPrimaryEnergy;
}

void MupTargetEnToLLProcess::Configure(const std::vector<G4String> &rootfiles, G4double xssf)
{
  fRun = ((RunAction *)G4RunManager::GetRunManager()->GetUserRunAction())->GetRun();
  delete fMupTargetEnToLL;
  fMupTargetEnToLL = new MupTargetEnToLL(rootfiles);
  fMuonPlus = G4MuonPlus::Definition();
  fElectron = G4Electron::Definition();
  fXSSF = xssf;
  fMuonMass = fMuonPlus->GetPDGMass();
}

G4double MupTargetEnToLLProcess::GetMeanFreePath(
    const G4Track &track, [[maybe_unused]] G4double previousStepSize, G4ForceCondition *condition)
{
  //G4cout << __FUNCTION__ << "(" << track.GetTrackID() << ")" << G4endl;
  G4double mfp = INFINITY;
  if(!fMupTargetEnToLL) {
    *condition = InActivated;
  } else {
    if(track.GetParticleDefinition()->GetPDGEncoding() == -13) {
      G4double xs = GetCrossSection(track.GetKineticEnergy(), track.GetMaterialCutsCouple());
      mfp = 1 / xs;
    }
    *condition = mfp == INFINITY ? NotForced : Forced;
  }
  //G4cout << __FUNCTION__ << "(" << track.GetTrackID() << ") -> " << mfp << G4endl;
  return mfp;
}
