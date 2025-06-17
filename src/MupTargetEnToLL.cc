#include "MupTargetEnToLL.hh"

#include <ExRootAnalysis/ExRootClasses.h>
#include <math.h>

#include <G4ParticleTable.hh>
#include <G4SystemOfUnits.hh>
#include <G4ios.hh>
#include <Randomize.hh>
#include <algorithm>
#include <filesystem>
#include <regex>

namespace fs = std::filesystem;

static G4double pb = 1e-12 * 1e-24 * cm2;

#define FPRE R"(.*)"  // [XXX]

MupTargetEnToLL::MupTargetEnToLL(const std::vector<G4String> &rootfiles)
{
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  e_mass = particleTable->FindParticle(11)->GetPDGMass();
  mu_mass = particleTable->FindParticle(13)->GetPDGMass();

  size_t nfile = rootfiles.size();
  G4cout << "Number of rootfiles: " << rootfiles.size() << G4endl;
  files.resize(nfile), points.resize(nfile);
  std::regex energy_regex(R"(MueZp_()" FPRE R"()GeV_Zp_()" FPRE R"()GeV.root$)");
  std::smatch match;
  for(size_t i = 0; i < nfile; ++i) {
    files[i] = new TFile(rootfiles[i]);
    if(!files[i]->IsOpen()) throw std::runtime_error("failed opening ROOT file: " + rootfiles[i]);
    auto &[energy, xs, tree, events, particles] = points[i];
    std::string basename = fs::path(rootfiles[i].c_str()).filename();
    if(std::regex_search(basename, match, energy_regex)) {
      energy = std::stod(match[1].str()) * GeV;
      G4cout << "Extracted energy: " << energy / MeV << " MeV" << G4endl;
    } else {
      throw std::runtime_error("unrecognized ROOT path: " + rootfiles[i]);
    }
    tree = (TTree *)files[i]->Get("LHEF");
    tree->SetBranchAddress("Event", &events);
    tree->SetBranchAddress("Particle", &particles);
    tree->GetEntry(0);
    static TRootLHEFEvent dummy;
    auto event = (TRootLHEFEvent *)events->At(0);
    xs = event->Weight * pb;
    G4cout << "Extracted cross section: " << xs / pb << " pb" << G4endl;
  }
  sort(points.begin(), points.end());  // [NOTE] Points do not match the order of rootfiles then.
}

MupTargetEnToLL::~MupTargetEnToLL()
{
  // [IMPORTANT NOTE]
  // Do NOT close opened ROOT files here.
  // This class may be destructed after TROOT.
}

void MupTargetEnToLL::Scatter(double mup_energy, double lp_out_alpha, double lp_out_phi, TLorentzVector p4_miss,
    G4ThreeVector &lp_out_p, G4ThreeVector &ln_out_p) const
{
  G4cout << "Entering Scatter function..." << G4endl;
  G4cout << "Computed mup_energy: " << mup_energy << G4endl;

  // Compute COM energy
  double e2_com = mu_mass * mu_mass + 2 * mup_energy * e_mass + e_mass * e_mass, e_com = sqrt(e2_com);

  // Compute Lorentz boost.
  double e = mup_energy + e_mass, e2 = e * e;
  double p2 = e2 - e2_com, p = sqrt(p2);
  double gamma = e / e_com, beta = p / e;

  // Compute momenta in COM frame.
  //double M = (-p4_miss.M() * p4_miss.M() + e_mass * e_mass - mu_mass * mu_mass) / (2 * p4_miss.E());
  //double A = sin(lp_out_alpha) * sin(lp_out_alpha) * cos(lp_out_phi) * cos(lp_out_phi) * p4_miss.Px() * p4_miss.Px()
  //    + sin(lp_out_alpha) * sin(lp_out_alpha) * sin(lp_out_phi) * sin(lp_out_phi) * p4_miss.Py() * p4_miss.Py()
  //    + cos(lp_out_alpha) * cos(lp_out_alpha) * p4_miss.Pz() * p4_miss.Pz()
  //    + 2 * sin(lp_out_alpha) * sin(lp_out_alpha) * cos(lp_out_phi) * sin(lp_out_phi) * p4_miss.Px() * p4_miss.Py()
  //    + 2 * sin(lp_out_alpha) * cos(lp_out_alpha) * cos(lp_out_phi) * p4_miss.Px() * p4_miss.Pz()
  //    + 2 * sin(lp_out_alpha) * cos(lp_out_alpha) * sin(lp_out_phi) * p4_miss.Py() * p4_miss.Pz() - 1;
  //double B = 2 * M * sin(lp_out_alpha) * cos(lp_out_phi) * p4_miss.Px()
  //    + 2 * M * sin(lp_out_alpha) * sin(lp_out_phi) * p4_miss.Py() + 2 * M * cos(lp_out_alpha) * p4_miss.Pz();
  //double C = M * M - mu_mass * mu_mass;
  //double disc = B * B - 4 * A * C;
  //G4cout << "Computed M: " << M << G4endl;
  //G4cout << "Computed A: " << A << G4endl;
  //G4cout << "Computed B: " << B << G4endl;
  //G4cout << "Computed C: " << C << G4endl;
  TVector3 dp_mu(0, 0, 1);
  dp_mu.SetTheta(lp_out_alpha), dp_mu.SetPhi(lp_out_phi);
  double A = p4_miss.P() * p4_miss.Vect().Unit().Dot(dp_mu);
  double B = (e2_com + mu_mass * mu_mass + p4_miss.M() * p4_miss.M() - e_mass * e_mass) / 2 - e_com * p4_miss.M();
  double C = p4_miss.E() - e_com;
  double disc = (A * B) * (A * B) - (A * A - C * C) * (B * B - C * C * mu_mass * mu_mass);
  double p_mu_1 = ((A * B) + sqrt(disc)) / (A * A - C * C);
  double p_mu_2 = ((A * B) - sqrt(disc)) / (A * A - C * C);
  G4cout << "Computed p_mu_1: " << p_mu_1 << G4endl;
  G4cout << "Computed p_mu_2: " << p_mu_2 << G4endl;

  //double lp_out_com_plus = (-B + sqrt(disc)) / (2 * A);
  //double lp_out_com_minus = (-B - sqrt(disc)) / (2 * A);
  double lp_out = (p_mu_1 >= 0) ? p_mu_1 : (p_mu_2 >= 0) ? p_mu_2 : -1;
  double lp_out_pt = lp_out * sin(lp_out_alpha);
  double lp_out_px = lp_out_pt * cos(lp_out_phi);
  double lp_out_py = lp_out_pt * sin(lp_out_phi);
  double lp_out_pz = lp_out * cos(lp_out_alpha);
  double lp_out_e = hypot(lp_out, mu_mass);
  double ln_out_px = -lp_out_px - p4_miss.Px();
  double ln_out_py = -lp_out_py - p4_miss.Py();
  double ln_out_pz = -lp_out_pz - p4_miss.Pz();
  double ln_out = hypot(hypot(ln_out_px, ln_out_py), ln_out_pz);
  double ln_out_e = hypot(ln_out, e_mass);
  double residual = lp_out_e + ln_out_e + p4_miss.E() - e_com;
  G4cout << "Computed lp_out: " << lp_out << G4endl;
  G4cout << "Computed ln_out: " << ln_out << G4endl;
  G4cout << "Computed residual: " << residual << G4endl;

  // Boost to LAB frame.
  lp_out_p.setX(lp_out_px);
  lp_out_p.setY(lp_out_py);
  lp_out_p.setZ(gamma * (lp_out_pz + beta * lp_out_e));
  ln_out_p.setX(ln_out_px);
  ln_out_p.setY(ln_out_py);
  ln_out_p.setZ(gamma * (ln_out_pz + beta * ln_out_e));

  G4cout << "-------------------------------------------" << G4endl;
  G4cout << "Lab frame outgoing muon momentum (lp_out_p):" << G4endl;
  G4cout << "  Px: " << lp_out_p.x() << "  Py: " << lp_out_p.y() << "  Pz: " << lp_out_p.z() << G4endl;
  G4cout << "Lab frame outgoing electron momentum (ln_out_p):" << G4endl;
  G4cout << "  Px: " << ln_out_p.x() << "  Py: " << ln_out_p.y() << "  Pz: " << ln_out_p.z() << G4endl;
  G4cout << "-------------------------------------------" << G4endl;
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
  if(right == points.end() && right != points.begin()) right = prev(right);  // [XXX] Temporarily enables extrapolation.
  if(right == points.begin() || right == points.end()) return 0.0;
  auto left = prev(right);
  auto &[l_mup_energy, l_xs, l_tree, l_events, l_particles] = *left;
  auto &[r_mup_energy, r_xs, r_tree, r_events, r_particles] = *right;
  auto [l_weight, r_weight] = linear_interp_weights(l_mup_energy, r_mup_energy, mup_energy);

  return l_weight * l_xs + r_weight * r_xs;
}

double MupTargetEnToLL::MinPrimaryEnergy() const { return points.empty() ? INFINITY : std::get<0>(points[0]); }

std::tuple<double, double, double, TLorentzVector> MupTargetEnToLL::Sample(double mup_energy) const
{
  // Locate end points.
  auto right =
      std::upper_bound(points.begin(), points.end(), mup_energy, [](double e, auto &p) { return e < std::get<0>(p); });
  if(right == points.end() && right != points.begin()) right = prev(right);  // [XXX] Temporarily enables extrapolation.
  if(right == points.begin() || right == points.end()) return { 0, NAN, NAN, { NAN, NAN, NAN, NAN } };
  auto left = prev(right);
  auto &[l_mup_energy, l_xs, l_tree, l_events, l_particles] = *left;
  auto &[r_mup_energy, r_xs, r_tree, r_events, r_particles] = *right;
  auto [l_weight, r_weight] = linear_interp_weights(l_mup_energy, r_mup_energy, mup_energy);

  // Interpolate for xs.
  double xs = l_weight * l_xs + r_weight * r_xs;

  // Use weighted sampling in place of histogram interpolation.
  auto &tree = G4UniformRand() < l_weight ? l_tree : r_tree;
  auto &particles = tree == l_tree ? l_particles : r_particles;
  auto [lp_out_alpha, lp_out_phi, e_miss] = Draw(tree, particles);

  return { xs, lp_out_alpha, lp_out_phi, e_miss };
}

std::tuple<double, double, TLorentzVector> MupTargetEnToLL::Draw(TTree *tree, TClonesArray *Particles) const
{
  size_t n = tree->GetEntries();
  size_t i = G4UniformRand() * n;
  tree->GetEntry(i);

  // Process: mu+ e- > mu+ e- zp, zp > v v
  int n_mu = 0, n_e = 0, n_nu = 0;
  TLorentzVector p4_mu, p4_e, p4_nu1, p4_nu2;
  for(int j = 0; j < Particles->GetEntries(); ++j) {
    static TRootLHEFParticle dummy;
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
