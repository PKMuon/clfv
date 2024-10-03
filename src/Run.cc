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

#include "Run.hh"

#include <TClonesArray.h>
#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>
#include <unistd.h>

#include <filesystem>

#include "G4LogicalVolumeStore.hh"
#include "G4Step.hh"
#include "G4Threading.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "Object.hh"

namespace fs = std::filesystem;

namespace {

struct ROOTInitializer {
  ROOTInitializer() { ROOT::EnableThreadSafety(); }
} rootInitializer [[maybe_unused]];

}  // namespace

G4String Run::fDirName = "tree/" + std::to_string(getpid());
G4String Run::fTreeName = "tree";
G4String Run::fTreeTitle = "tree";

class Run::Manager {
public:
  Manager();
  ~Manager();

  void Branch(TTree *tree);
  void PreFill();
  void Reset();
  void AddTrack(const G4Track *track);
  void AddStep(const G4Step *track);
  void AddScatter(double probability, double xs);
  void SaveCuts();

private:
  TFile *fFile;
  TTree *fCuts;
  TClonesArray Tracks;
  TClonesArray Cuts;
  Double_t EnergyDeposit, NonIonizingEnergyDeposit;
  Double_t ScatterProbability, ScatterXS;
};

Run::Run()
{
  G4String filename = fDirName + "/" + std::to_string(G4Threading::G4GetThreadId()) + ".root";
  auto parent = fs::path(filename.c_str()).parent_path();
  if(!parent.empty()) fs::create_directories(parent);

  fFile = TFile::Open(filename, "RECREATE");
  if(!fFile->IsOpen()) {
    G4cerr << "Error opening file " << filename << G4endl;
    exit(1);
  }
  fTree = new TTree(fTreeName, fTreeTitle);
  fManager = new Manager;
  fManager->Branch(fTree);
}

Run::~Run()
{
  delete fManager;
  fFile->cd();
  fTree->Write(NULL, fTree->kOverwrite);
  fFile->Close();
}

void Run::AutoSave()
{
  fTree->AutoSave("SaveSelf, Overwrite");
  fManager->SaveCuts();
}

void Run::FillAndReset()
{
  fManager->PreFill();
  fTree->Fill();
  fManager->Reset();
}

void Run::AddTrack(const G4Track *track) { fManager->AddTrack(track); }

void Run::AddStep(const G4Step *step) { fManager->AddStep(step); }

void Run::AddScatter(double probability, double xs) { fManager->AddScatter(probability, xs); }

Run::Manager::Manager()
    : Tracks("Track"), Cuts("Cuts"), EnergyDeposit(0), NonIonizingEnergyDeposit(0), ScatterProbability(0), ScatterXS(0)
{
  fFile = NULL;
  fCuts = NULL;
}

Run::Manager::~Manager()
{
  if(fFile == NULL) return;
  fFile->cd();
  fCuts->Write(NULL, fCuts->kOverwrite);
}

void Run::Manager::Branch(TTree *tree)
{
  tree->Branch("Tracks", &Tracks);
  tree->Branch("EnergyDeposit", &EnergyDeposit);
  tree->Branch("NonIonizingEnergyDeposit", &NonIonizingEnergyDeposit);
  tree->Branch("ScatterProbability", &ScatterProbability);
  tree->Branch("ScatterXS", &ScatterXS);

  fFile = tree->GetCurrentFile();
  fFile->cd();
  fCuts = new TTree("cuts", "cuts");
  fCuts->Branch("Cuts", &Cuts);
}

void Run::Manager::PreFill()
{
  //// Inplace index sort.
  //Int_t n = Tracks.GetEntries();
  //for(Int_t i = 0; i < n; ++i) {
  //  auto track = (Track *)Tracks[i];
  //  while(track->Id - 1 != i) {
  //    if(track->Id <= 0 || track->Id > n) {
  //      throw std::runtime_error("invalid track ID: " + std::to_string(track->Id));
  //    }
  //    if(((Track *)Tracks[track->Id - 1])->Id == track->Id) {
  //      throw std::runtime_error("duplicate track ID: " + std::to_string(track->Id));
  //    }
  //    TObject *object = track;
  //    std::swap(object, Tracks[track->Id - 1]);
  //    track = (Track *)object;
  //  }
  //  Tracks[i] = track;
  //}
}

void Run::Manager::Reset()
{
  Tracks.Clear();
  EnergyDeposit = 0;
  NonIonizingEnergyDeposit = 0;
  ScatterProbability = 0;
  ScatterXS = 0;
}

void Run::Manager::AddTrack(const G4Track *track) { *(Track *)Tracks.ConstructedAt(Tracks.GetEntries()) = *track; }

void Run::Manager::AddStep(const G4Step *step)
{
  EnergyDeposit += step->GetTotalEnergyDeposit();
  NonIonizingEnergyDeposit += step->GetNonIonizingEnergyDeposit();
}

void Run::Manager::AddScatter(double probability, double xs)
{
  ScatterProbability = probability;
  ScatterXS = xs;
}

void Run::Manager::SaveCuts()
{
  *(::Cuts *)Cuts.ConstructedAt(0) = *G4LogicalVolumeStore::GetInstance()->GetVolume("world");
  fCuts->Fill();
  fCuts->AutoSave("SaveSelf, Overwrite");
}
