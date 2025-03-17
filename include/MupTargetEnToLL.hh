#pragma once
#include <TClonesArray.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>

#include <G4ThreeVector.hh>
#include <memory>
#include <tuple>
#include <vector>
#include <regex>

class MupTargetEnToLL {
public:
  MupTargetEnToLL(const std::vector<G4String> &rootfiles);
  ~MupTargetEnToLL();

  // For incoming mu+ beams towards z+.
  void Scatter(double mup_energy, double lp_out_alpha, double lp_out_phi, TLorentzVector p4_miss,
      G4ThreeVector &lp_out_p, G4ThreeVector &ln_out_p) const;

  // For arbitrary incoming mu+ beams.
  void Scatter(G4ThreeVector &lp_p, double lp_out_alpha, double lp_out_phi, TLorentzVector p4_miss,
      G4ThreeVector &ln_out_p) const;

  // Request sampling.
  double Scatter(G4ThreeVector &lp_p, G4ThreeVector &ln_out_p) const;

  double CrossSection(double mup_energy) const;
  double MinPrimaryEnergy() const;

private:
  double e_mass, mu_mass;
  std::vector<TFile *> files;

  // (mup_energy, xs, tree, events, particles)
  std::vector<std::tuple<double, double, TTree *, TClonesArray *, TClonesArray *>> points;

  std::regex energy_regex;

  // (xs, lp_out_alpha, ln_out_phi, p4_miss)
  std::tuple<double, double, double, TLorentzVector> Sample(double mup_energy) const;
  std::tuple<double, double, TLorentzVector> Draw(TTree *, TClonesArray *) const;
};
