#pragma once
#include <G4ThreeVector.hh>

class MupTargetEnToLL {
public:
  MupTargetEnToLL(int l_pid);

  // For incoming mu+ beams towards z+.
  void Scatter(double mup_energy, double mup_out_alpha, G4ThreeVector &mup_out_p, G4ThreeVector &mun_out_p);

  // For arbitrary incoming mu+ beams.
  void Scatter(G4ThreeVector &mup_p, double mup_out_alpha, G4ThreeVector &mun_out_p);

private:
  double e_mass, mu_mass, l_mass;
};
