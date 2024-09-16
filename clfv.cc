#include "ScatterProcess.hh"
#include <Randomize.hh>
#include <G4SystemOfUnits.hh>
#include <G4ios.hh>
#include <particle.h>

#ifndef DATA_DIR
#define DATA_DIR  "."
#endif  /* DATA_DIR */

#ifndef PDG_PATH
#define PDG_PATH  ""
#endif  /* PDG_PATH */

int main()
{
  ParticleDatabase db(PDG_PATH);
  double muon_mass = (double)db.query(13, "mass") * MeV;

  MupTargetEnToLL process(13, DATA_DIR "/muemumu_ki_2.200000000000000011e-01.root");

  for(int i = 0; i < 100; ++i) {
    double muon_energy = (G4UniformRand() * 99 + 1) * GeV;
    G4ThreeVector lp_p(0, 0, sqrt(muon_energy*muon_energy - muon_mass*muon_mass)), ln_p;
    G4cout << lp_p << G4endl;
    process.Scatter(lp_p, ln_p);
    G4cout << lp_p << "\t" << ln_p << G4endl;
  }

  return 0;
}
