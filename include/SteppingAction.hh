//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class MupTargetEnToLL;
class EventAction;
class Run;

class SteppingAction : public G4UserSteppingAction {
public:
  SteppingAction(EventAction *eventAction);
  ~SteppingAction() override = default;

  void UserSteppingAction(const G4Step *) override;

  void SetMupTargetEnToEE(double probability, const char *points_file) { SetMupTargetEnToLL(0, probability, points_file); }
  void SetMupTargetEnToMuMu(double probability, const char *points_file) { SetMupTargetEnToLL(1, probability, points_file); }

private:
  [[maybe_unused]] EventAction *fEventAction;
  [[maybe_unused]] Run *fRun;
  G4double fScatterProbability[2];  // e, mu
  MupTargetEnToLL *fScatterProcess[2];  // e, mu

  void SetMupTargetEnToLL(size_t index, double probability, const char *points_file);
};

#endif