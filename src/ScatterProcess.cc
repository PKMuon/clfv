#include "ScatterProcess.hh"
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>
#include <particle.h>
#include <math.h>
//#include <iostream>
//#include <iomanip>

#ifndef PDG_PATH
#define PDG_PATH ""
#endif  /* PDG_PATH */

MupTargetEnToLL::MupTargetEnToLL(int l_pid)
{
  ParticleDatabase db(PDG_PATH);
  e_mass = (double)db.query(11, "mass") * MeV;
  mu_mass = (double)db.query(13, "mass") * MeV;
  l_mass = (double)db.query(l_pid, "mass") * MeV;
}

static double compute_alpha(double beta, double gamma, double e_com, double mup_p_com, double mup_out_theta)
{
  double A = tan(mup_out_theta) * gamma;
  double B = beta * e_com / (2 * mup_p_com);
  //std::cout << std::scientific << std::setprecision(15) << __func__ << ": A=" << A << "  B=" << B << std::endl;
  double alpha = 0;
  for(;;) {
    double alpha_new = atan(A * (1 + B / cos(alpha)));
    //std::cout << __func__ << ": alpha=" << alpha_new << std::endl;
    if(fabs(alpha - alpha_new) < DBL_EPSILON * alpha_new) return alpha_new;
    alpha = alpha_new;
  }
}

void MupTargetEnToLL::Scatter(double mup_energy, double mup_out_theta, G4ThreeVector &mup_out_p, G4ThreeVector &mun_out_p)
{
  // Compute COM energy
  double e2_com = mu_mass*mu_mass + 2*mup_energy*e_mass + e_mass*e_mass, e_com = sqrt(e2_com);

  // Compute Lorentz boost.
  double e = mup_energy + e_mass, e2 = e*e;
  double p2 = e2 - e2_com, p = sqrt(p2);
  double gamma = e / e_com, beta = p / e;

  // Compute mu+ momentum in COM frame.
  double mup_p_com = sqrt(e2_com / 4 - mu_mass*mu_mass);

  // Solve equation for alpha.
  double alpha = compute_alpha(beta, gamma, e_com, mup_p_com, mup_out_theta);

  // Write results via alpha.
  double mup_out_pt = mup_p_com * sin(alpha);
  double mup_out_phi = G4UniformRand() * (2 * M_PI);
  mup_out_p.setX(mup_out_pt * cos(mup_out_phi));
  mup_out_p.setY(mup_out_pt * sin(mup_out_phi));
  mup_out_p.setZ(gamma * (mup_p_com * cos(alpha) + beta * e_com / 2));
  mun_out_p.setX(-mup_out_p.getX());
  mun_out_p.setY(-mup_out_p.getY());
  mun_out_p.setZ(p - mup_out_p.getZ());
}
