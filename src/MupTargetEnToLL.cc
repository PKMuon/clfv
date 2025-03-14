#include "MupTargetEnToLL.hh"

#include <ExRootAnalysis/ExRootClasses.h>
#include <math.h>

#include <G4ParticleTable.hh>
#include <G4SystemOfUnits.hh>
#include <G4ios.hh>
#include <Randomize.hh>
#include <algorithm>

static G4double pb = 1e-12 * 1e-24 * cm2;

MupTargetEnToLL::MupTargetEnToLL(const std::vector<G4String> &rootfiles)
{
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  e_mass = particleTable->FindParticle(11)->GetPDGMass();
  mu_mass = particleTable->FindParticle(13)->GetPDGMass();
  Events = Particles = nullptr;

  size_t nfile = rootfiles.size();
  files.resize(nfile), points.resize(nfile);
  for(size_t i = 0; i < nfile; ++i) {
    files[i].reset(new TFile(rootfiles[i]));
    if(!files[i]->IsOpen()) throw std::runtime_error("Error opening file " + rootfiles[i]);
    auto &[energy, xs, tree] = points[i];
    energy = 0 * GeV;  // [TODO] Parse energy from filename.
    tree = (TTree *)files[i]->Get("LHEF");
    tree->SetBranchAddress("Events", &Events);
    tree->SetBranchAddress("Particles", &Particles);
    if(!tree->GetEntry(0)) throw std::runtime_error("Error reading file " + rootfiles[i]);
    auto event = (TRootLHEFEvent *)Events->At(0);
    xs = event->Weight * pb;
  }
  sort(points.begin(), points.end());  // [NOTE] Points do not match the order of rootfiles then.
}

MupTargetEnToLL::~MupTargetEnToLL()
{
  for(auto [energy, xs, tree] : points) {
    tree->SetBranchAddress("Events", NULL);
    tree->SetBranchAddress("Particles", NULL);
  }
  delete Events;
  delete Particles;
}

void MupTargetEnToLL::Scatter(double mup_energy, double lp_out_alpha, double lp_out_phi, TLorentzVector p4_miss,
    G4ThreeVector &lp_out_p, G4ThreeVector &ln_out_p) const
{
  // Compute COM energy
  double e2_com = mu_mass * mu_mass + 2 * mup_energy * e_mass + e_mass * e_mass, e_com = sqrt(e2_com);

  // Compute Lorentz boost.
  double e = mup_energy + e_mass, e2 = e * e;
  double p2 = e2 - e2_com, p = sqrt(p2);
  double gamma = e / e_com, beta = p / e;

  // Compute momenta in COM frame.  // [TODO]
  assert(lp_out_alpha || lp_out_phi || p4_miss.E());  // [TODO] Remove this line after implementation.

  // Boost to LAB frame. [TODO]
  lp_out_p.setX(0.0);
  lp_out_p.setY(0.0);
  lp_out_p.setZ(gamma * (0.0 + beta * 0.0));
  ln_out_p.setX(0.0);
  ln_out_p.setY(0.0);
  ln_out_p.setZ(gamma * (0.0 + beta * 0.0));
}

void MupTargetEnToLL::Scatter(
    G4ThreeVector &lp_p, double lp_out_alpha, double lp_out_phi, TLorentzVector p4_miss, G4ThreeVector &ln_out_p) const
{
  // Save original direction.
  double theta = lp_p.getTheta();
  double phi = lp_p.getPhi();

  // Calculate results in the beam-z frame.
  double mup_energy = sqrt(lp_p.mag2() + mu_mass * mu_mass);
  Scatter(mup_energy, lp_out_alpha, lp_out_phi, p4_miss, lp_p, ln_out_p);

  // Rotate back to the original frame.
  lp_p.rotateY(theta);
  lp_p.rotateZ(phi);
  ln_out_p.rotateY(theta);
  ln_out_p.rotateZ(phi);
}

double MupTargetEnToLL::Scatter(G4ThreeVector &lp_p, G4ThreeVector &ln_out_p) const
{
  double mup_energy = sqrt(lp_p.mag2() + mu_mass * mu_mass);
  auto [xs, lp_out_alpha, lp_out_phi, e_miss] = Sample(mup_energy);
  if(xs == 0) {
    ln_out_p = { NAN, NAN, NAN };
    return 0;
  }
  Scatter(lp_p, lp_out_alpha, lp_out_phi, e_miss, ln_out_p);
  return xs;
}

static std::pair<double, double> linear_interp_weights(double x1, double x2, double x)
{
  return {
    (x - x2) / (x1 - x2),
    (x1 - x) / (x1 - x2),
  };
}

double MupTargetEnToLL::CrossSection(double mup_energy) const
{
  // Locate end points.
  auto right =
      std::upper_bound(points.begin(), points.end(), mup_energy, [](double e, auto &p) { return e < std::get<0>(p); });
  if(right == points.begin() || right == points.end()) return 0.0;
  auto left = prev(right);
  auto &[l_mup_energy, l_xs, l_tree] = *left;
  auto &[r_mup_energy, r_xs, r_tree] = *right;
  auto [l_weight, r_weight] = linear_interp_weights(l_mup_energy, r_mup_energy, mup_energy);

  return l_weight * l_xs + r_weight * r_xs;
}

double MupTargetEnToLL::MinPrimaryEnergy() const { return points.empty() ? INFINITY : std::get<0>(points[0]); }

std::tuple<double, double, double, TLorentzVector> MupTargetEnToLL::Sample(double mup_energy) const
{
  // Locate end points.
  auto right =
      std::upper_bound(points.begin(), points.end(), mup_energy, [](double e, auto &p) { return e < std::get<0>(p); });
  if(right == points.begin() || right == points.end()) return { 0, NAN, NAN, { NAN, NAN, NAN, NAN } };
  auto left = prev(right);
  auto &[l_mup_energy, l_xs, l_tree] = *left;
  auto &[r_mup_energy, r_xs, r_tree] = *right;
  auto [l_weight, r_weight] = linear_interp_weights(l_mup_energy, r_mup_energy, mup_energy);

  // Interpolate for xs.
  double xs = l_weight * l_xs + r_weight * r_xs;

  // Use weighted sampling in place of histogram interpolation.
  auto &tree = G4UniformRand() < l_weight ? l_tree : r_tree;
  auto [lp_out_alpha, lp_out_phi, e_miss] = Draw(tree);

  return { xs, lp_out_alpha, lp_out_phi, e_miss };
}

std::tuple<double, double, TLorentzVector> MupTargetEnToLL::Draw(TTree *tree) const
{
  size_t n = tree->GetEntries();
  size_t i = G4UniformRand() * n;
  tree->GetEntry(i);

  // Process: mu+ e- > mu+ e- zp
  int n_mu = 0, n_e = 0, n_nu = 0;
  TLorentzVector p4_mu, p4_e, p4_nu1, p4_nu2;
  for(int j = 0; j < Particles->GetEntries(); ++j) {
    auto particle = (TRootLHEFParticle *)Particles->UncheckedAt(j);
    if(particle->Status != 1) continue;
    if(particle->PID == -13) {
      p4_mu.SetPxPyPzE(particle->Px, particle->Py, particle->Pz, particle->E);
      ++n_mu;
    } else if(particle->PID == 11) {
      p4_e.SetPxPyPzE(particle->Px, particle->Py, particle->Pz, particle->E);
      ++n_e;
    } else {
      (n_nu ? p4_nu2 : p4_nu1).SetPxPyPzE(particle->Px, particle->Py, particle->Pz, particle->E);
      ++n_nu;
    }
  }
  assert(n_mu == 1 && n_e == 1 && n_nu == 2);
  TLorentzVector p4_zp = p4_nu1 + p4_nu2;

  TVector3 b = (p4_mu + p4_e + p4_zp).BoostVector();
  p4_mu.Boost(-b), p4_e.Boost(-b), p4_zp.Boost(-b);
  return { p4_mu.Theta(), p4_mu.Phi(), p4_zp };
}
