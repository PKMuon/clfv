#pragma once
#include <G4ThreeVector.hh>
#include <vector>
#include <utility>
#include <tuple>
#include <memory>
#include <TH1.h>

class G4Track;
class G4ParticleDefinition;

class MupTargetEnToLL {
public:
  MupTargetEnToLL(int l_pid, const char *points_file);

  // For incoming mu+ beams towards z+.
  void Scatter(double mup_energy, double lp_out_alpha, G4ThreeVector &lp_out_p, G4ThreeVector &ln_out_p) const;

  // For arbitrary incoming mu+ beams.
  void Scatter(G4ThreeVector &lp_p, double lp_out_alpha, G4ThreeVector &ln_out_p) const;

  // Request lp_out_alpha sampling.
  double Scatter(G4ThreeVector &lp_p, G4ThreeVector &ln_out_p) const;

  // Kill original lp_track, and produce lp_track and ln_out_track.
  double Scatter(G4Track *&lp_track, G4Track *&ln_out_track) const;

private:
  G4ParticleDefinition *lp_def, *ln_def;
  double e_mass, mu_mass, l_mass;
  std::vector<std::tuple<double, double, std::unique_ptr<TH1>>> points;  // (mup_energy, xs, lp_out_alpha)

  void LoadPoints(const char *points_file);
  std::pair<double, double> Sample(double mup_energy) const;  // (xs, lp_out_alpha)
};
