#include "ScatterProcess.hh"
#include <G4ios.hh>
#include <G4ParticleTable.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4DynamicParticle.hh>
#include <G4EventManager.hh>
#include <G4StackManager.hh>
#include <Randomize.hh>
#include <algorithm>
#include <math.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

#ifndef PDG_PATH
#define PDG_PATH  ""
#endif  /* PDG_PATH */

MupTargetEnToLL::MupTargetEnToLL(int l_pid, const char *points_file)
{
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  lp_def = particleTable->FindParticle(-abs(l_pid));
  ln_def = particleTable->FindParticle(+abs(l_pid));
  e_mass = particleTable->FindParticle(11)->GetPDGMass();
  mu_mass = particleTable->FindParticle(13)->GetPDGMass();
  l_mass = ln_def->GetPDGMass();
  LoadPoints(points_file);
}

void MupTargetEnToLL::Scatter(double mup_energy, double lp_out_alpha, G4ThreeVector &lp_out_p, G4ThreeVector &ln_out_p) const
{
  // Compute COM energy
  double e2_com = mu_mass*mu_mass + 2*mup_energy*e_mass + e_mass*e_mass, e_com = sqrt(e2_com);

  // Compute Lorentz boost.
  double e = mup_energy + e_mass, e2 = e*e;
  double p2 = e2 - e2_com, p = sqrt(p2);
  double gamma = e / e_com, beta = p / e;

  // Compute l+ momentum in COM frame.
  double lp_p_com = sqrt(e2_com / 4 - l_mass*l_mass);

  // Write results via alpha.
  double lp_out_pt = lp_p_com * sin(lp_out_alpha);
  double lp_out_phi = G4UniformRand() * (2 * M_PI);
  lp_out_p.setX(lp_out_pt * cos(lp_out_phi));
  lp_out_p.setY(lp_out_pt * sin(lp_out_phi));
  lp_out_p.setZ(gamma * (lp_p_com * cos(lp_out_alpha) + beta * e_com / 2));
  ln_out_p.setX(-lp_out_p.getX());
  ln_out_p.setY(-lp_out_p.getY());
  ln_out_p.setZ(p - lp_out_p.getZ());
}

void MupTargetEnToLL::Scatter(G4ThreeVector &lp_p, double lp_out_alpha, G4ThreeVector &ln_out_p) const
{
  // Save original direction.
  double theta = lp_p.getTheta();
  double phi = lp_p.getPhi();

  // Calculate results in the beam-z frame.
  double mup_energy = sqrt(lp_p.mag2() + mu_mass*mu_mass);
  Scatter(mup_energy, lp_out_alpha, lp_p, ln_out_p);

  // Rotate back to the original frame.
  lp_p.rotateY(theta);
  lp_p.rotateZ(phi);
  ln_out_p.rotateY(theta);
  ln_out_p.rotateZ(phi);
}

double MupTargetEnToLL::Scatter(G4ThreeVector &lp_p, G4ThreeVector &ln_out_p) const
{
  double mup_energy = sqrt(lp_p.mag2() + mu_mass*mu_mass);
  auto [xs, lp_out_alpha] = Sample(mup_energy);
  if(xs == 0) {
    ln_out_p = {NAN, NAN, NAN};
    return 0;
  }
  Scatter(lp_p, lp_out_alpha, ln_out_p);
  return xs;
}

void MupTargetEnToLL::LoadPoints(const char *points_file)
{
  points.clear();

  G4cout << "Loading points from " << points_file << G4endl;
  auto file = TFile::Open(points_file);
  if(!file->IsOpen()) G4Exception("MupTargetEnToLL::LoadPoints", "FileOpenError", FatalException, "Failed to open points file.");
  auto tree = (TTree *)file->Get("ki_points");
  if(!tree) G4Exception("MupTargetEnToLL::LoadPoints", "TreeGetError", FatalException, "Failed to get points tree.");

  double mup_energy, xs;
  TH1 *lp_out_alpha = nullptr;
  tree->SetBranchAddress("muon_energy", &mup_energy);
  tree->SetBranchAddress("xs", &xs);
  tree->SetBranchAddress("alpha_pos", &lp_out_alpha);

  Long64_t ientry;
  for(ientry = 0; tree->GetEntry(ientry); ++ientry) {
    auto lp_out_alpha_clone = (TH1 *)lp_out_alpha->Clone();
    lp_out_alpha_clone->SetDirectory(nullptr);
    points.emplace_back(mup_energy * GeV, xs, lp_out_alpha_clone);
  }
  G4cout << "Loaded " << ientry << " points" << G4endl;
  std::sort(points.begin(), points.end());

  file->Close();
}

static std::pair<double, double> linear_interp_weights(double x1, double x2, double x)
{
  return {
    (x - x2) / (x1 - x2),
    (x1 - x) / (x1 - x2),
  };
}

std::pair<double, double> MupTargetEnToLL::Sample(double mup_energy) const
{
  // Locate end points.
  auto right = std::upper_bound(points.begin(), points.end(), mup_energy, [](double e, auto &p) { return e < std::get<0>(p); });
  if(right == points.begin() || right == points.end()) return {0, NAN};
  auto left = prev(right);
  auto &[l_mup_energy, l_xs, l_lp_out_alpha] = *left;
  auto &[r_mup_energy, r_xs, r_rp_out_arpha] = *right;
  auto [l_weight, r_weight] = linear_interp_weights(l_mup_energy, r_mup_energy, mup_energy);

  // Interpolate for xs.
  double xs = l_weight * l_xs + r_weight * r_xs;

  // Use weighted sampling in place of histogram interpolation.
  auto &lp_out_alpha_hist = G4UniformRand() < l_weight ? l_lp_out_alpha : r_rp_out_arpha;
  double lp_out_alpha = lp_out_alpha_hist->GetRandom();

  return {xs, lp_out_alpha};
}

double MupTargetEnToLL::Scatter(G4Track *&lp_track, G4Track *&ln_out_track) const
{
  // Compute xs and scattered momenta.
  G4ThreeVector lp_p = lp_track->GetMomentum(), ln_out_p;
  double xs = Scatter(lp_p, ln_out_p);
  if(xs == 0) return 0;

  // Collect track information and kill the original track.
  const G4ThreeVector &position = lp_track->GetPosition();
  const G4TouchableHandle &handle = lp_track->GetTouchableHandle();
  G4double globalTime = lp_track->GetGlobalTime();
  G4int parent_id = lp_track->GetTrackID();
  lp_track->SetTrackStatus(G4TrackStatus::fStopAndKill);

  // Create scattered tracks.
  auto lp_particle = new G4DynamicParticle(lp_def, lp_p);
  auto ln_particle = new G4DynamicParticle(ln_def, ln_out_p);
  lp_track = new G4Track(lp_particle, globalTime, position);
  ln_out_track = new G4Track(ln_particle, globalTime, position);
  lp_track->SetTrackID(0x10001);
  ln_out_track->SetTrackID(0x10002);

  // Register scattered tracks.
  G4StackManager *stackManager = G4EventManager::GetEventManager()->GetStackManager();
  for(G4Track *track : {lp_track, ln_out_track}) {
    track->SetParentID(parent_id);
    track->SetTouchableHandle(handle);
    stackManager->PushOneTrack(track);
  }

  return xs;
}
