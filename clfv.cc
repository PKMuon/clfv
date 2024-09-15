#include "ScatterProcess.hh"
#include <G4ios.hh>
#include <G4SystemOfUnits.hh>

int main()
{
  MupTargetEnToLL scatter(13);
  G4ThreeVector mup_out_p, mun_out_p;
  scatter.Scatter(40.0 * GeV, 2e-3, mup_out_p, mun_out_p);
  G4cout << mup_out_p << G4endl;
  G4cout << mun_out_p << G4endl;

  return 0;
}
