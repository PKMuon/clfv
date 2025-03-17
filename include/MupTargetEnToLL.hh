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
  std::vector<std::unique_ptr<TFile>> files;
  std::vector<std::tuple<double, double, TTree *>> points;  // (mup_energy, xs, tree)
  TClonesArray *Events, *Particles;
  
  std::regex energy_regex;

  // (xs, lp_out_alpha, ln_out_phi, p4_miss)
  std::tuple<double, double, double, TLorentzVector> Sample(double mup_energy) const;
  std::tuple<double, double, TLorentzVector> Draw(TTree *) const;
  
  //int Event_Nparticles;
  //double Event_Weight;

  //static  const int MAX_PARTICLES = 1000000; // 设定最大粒子数
  //Int_t Particle_PID[MAX_PARTICLES], Particle_Status[MAX_PARTICLES];
  //Double_t Particle_Px[MAX_PARTICLES], Particle_Py[MAX_PARTICLES], Particle_Pz[MAX_PARTICLES], Particle_E[MAX_PARTICLES];

  //Int_t Particle_size = 0; // 粒子数
};
