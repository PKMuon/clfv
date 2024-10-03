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

#ifndef Run_h
#define Run_h 1

#include <Rtypes.h>

#include "globals.hh"

class TFile;
class TTree;
class G4Track;
class G4Step;

class Run {
public:
  Run();
  virtual ~Run();

  void AutoSave();
  void FillAndReset();
  void AddTrack(const G4Track *);
  void AddStep(const G4Step *);
  void AddScatter(double probability, double xs);

private:
  class Manager;
  Manager *fManager;

protected:
  static G4String fDirName;
  static G4String fTreeName;
  static G4String fTreeTitle;
  TFile *fFile;
  TTree *fTree;
};

#endif
