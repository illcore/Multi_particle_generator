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
//
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <iostream>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//29-bin BESS-TeV spectrum - by AIT
void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
//Discrete probabilities for protons - by AIT
double A[]= {0, 0.139931881, 0.128737331, 0.119142002, 0.108347257, 0.096353095,
0.083159518, 0.069566135, 0.057971779, 0.04757684, 0.038061472, 0.029385695,
0.022509043, 0.017351553, 0.012593869, 0.008995621, 0.006356905, 0.00447782,
0.003082499, 0.002130963, 0.001451294, 0.000991517, 0.000647685, 0.000435788,
0.00028666, 0.000193506, 0.000125939, 8.27597E-05, 5.35739E-05}; 
//Discrete energies for protons
double B[]= {0.0, 1080, 1260, 1470, 1710, 2000,	2330, 2710, 3160,	
3690, 4300, 5010, 5840,	6810, 7930, 9250, 10800, 12600,	14700,	
17100, 19900, 23200, 27100, 31600, 36800, 42900, 50000, 58300,	
68000};
///////////////////////////////////////////////////////////////////////////////////////
//Discrete probabilities for alpha - by AIT
double C[]= {0, 0.180861048, 0.165616056, 0.139283796, 0.11121915, 0.090777001,
0.075185532, 0.062712356, 0.047467363, 0.033850813, 0.026574794, 0.019783843,
0.013789789, 0.009805302, 0.007172076, 0.005127861, 0.003534066, 0.002342185,
0.001631907, 0.00113298, 0.000737996, 0.000505857, 0.000335043, 0.000204075,
0.00014656, 9.8746E-05, 6.37518E-05, 3.25688E-05, 7.48391E-06}; 
//Discrete energies for alpha
double D[]= {0.0, 1080, 1260, 1470, 1710, 2000,	2330, 2710, 3160, 3690,	4290,	
4980, 5840, 6800, 7940,	9240, 10800, 12600, 14700, 17000, 19900, 23200,	
27100, 31400, 36700, 42900, 49900, 62500, 86100};
double Grid1[29];
double Grid2[29];
G4double protonprob=0.893;
G4double alphaprob=0.107;
double sum1=0;
for(int x=0; x < sizeof(Grid1)/sizeof(Grid1[0]); x++){
  sum1=sum1+A[x];
  Grid1[x]=sum1;
  std::ofstream GridFile1;
  GridFile1.open("Probability_grid_proton.txt", std::ios::app);
  GridFile1 <<  Grid1[x] << G4endl;
  GridFile1.close();
  }
double sum2=0;  
for(int x=0; x < sizeof(Grid2)/sizeof(Grid2[0]); x++){
  sum2=sum2+C[x];
  Grid2[x]=sum2;
  std::ofstream GridFile2;
  GridFile2.open("Probability_grid_alpha.txt", std::ios::app);
  GridFile2 <<  Grid2[x] << G4endl;
  GridFile2.close();
  }          
// default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  for (int n_particle = 1; n_particle < 100000; n_particle++){  
  G4double pseudonumberI=G4UniformRand();
  if (pseudonumberI <= protonprob){
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="proton");
  fParticleGun->SetParticleDefinition(particle);  
  G4double y0 = 85*cm;
  G4double z0 = 0.5*cm;
  G4double x0 = 0.5*cm; 
  x0 = -x0+2*x0*G4UniformRand();
  z0 = -z0+2*z0*G4UniformRand();
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  G4double Energy; //Just for initialization
  G4double pseudonumberII=G4UniformRand();
  for (int i=0; i < sizeof(Grid1)/sizeof(Grid1[0]); i++){
  if(pseudonumberII > Grid1[i] && pseudonumberII <= Grid1[i+1]){
  Energy=B[i+1];
  std::ofstream EnergyFile1;
  EnergyFile1.open("Energy_proton.txt", std::ios::app);
  EnergyFile1 <<  Energy << G4endl;
  EnergyFile1.close();
  fParticleGun->SetParticleEnergy(Energy);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,-1,0));
  }
  }
  }
  if (pseudonumberI > 1-alphaprob && pseudonumberI <= 1){
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="alpha");
    fParticleGun->SetParticleDefinition(particle);
  G4double y0 = 85*cm;
  G4double z0 = 0.5*cm;
  G4double x0 = 0.5*cm; 
  x0 = -x0+2*x0*G4UniformRand();
  z0 = -z0+2*z0*G4UniformRand();
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  G4double Energy; //Just for initialization
  G4double pseudonumberIII=G4UniformRand();
  for (int i=0; i < sizeof(Grid2)/sizeof(Grid2[0]); i++){
  if(pseudonumberIII > Grid2[i] && pseudonumberIII <= Grid2[i+1]){
  Energy=D[i+1];
  std::ofstream EnergyFile2;
  EnergyFile2.open("Energy_alpha.txt", std::ios::app);
  EnergyFile2 <<  Energy << G4endl;
  EnergyFile2.close();
  fParticleGun->SetParticleEnergy(Energy);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,-1,0));  
  }
  }
  }
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

