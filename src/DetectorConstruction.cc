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
// $Id: DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "SensitiveDetector.hh"
#include "Hit.hh"

#include "G4SDManager.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4AssemblyVolume.hh"
#include <fstream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{

  // ------------ some initial settings -----------
  // get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  // ------------ envelope and world size ------------
  G4double env_sizeXY = 2600*cm, env_sizeZ = 1111*cm;
  G4double world_sizeXY = 1.1*env_sizeXY;
  G4double world_sizeZ  = 1.1*env_sizeZ;

  // material: vacuum - for world and envelope
  G4Material* vacuum_material = new G4Material("Vacuum",
              1.0 , 1.01*g/mole, 1.0E-25*g/cm3,
              kStateGas, 2.73*kelvin, 3.0E-18*pascal );

  // ---------------- place world and envelope --------------
  // world
  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        vacuum_material,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  // envelope
  G4Box* solidEnv =
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size

  G4LogicalVolume* logicEnv =
    new G4LogicalVolume(solidEnv,            //its solid
                        vacuum_material,             //its material
                        "Envelope");         //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  // ---------------- get detector config --------------

  std::fstream configFile;
  configFile.open("../src/config.txt", std::ios_base::in);

  // Initialize variables
  G4double det1_thickness;

  // Load into variables
  configFile >> det1_thickness;

  configFile.close();

  // world offset
  G4double world_offset = env_sizeZ*.45;

  // Dimensions for detectors (detector 1 and 2 use the same planar dimensions)
  G4double detector_dimX = 6.3*cm;
  G4double detector_dimY = 6.3*cm;

  //G4double detector_dimX = 1.364*cm;
  //G4double detector_dimY = 1.364*cm;

  G4double detector1_thickness = det1_thickness*um;
  G4double detector2_thickness = 140.0*um;

  G4double distance_between_detectors = 30.0*mm;

  // Window dimensions
  G4double window_dimX    = 6.3*cm; 
  G4double window_dimY    = 6.3*cm;  
  G4double window_thickness = 100.0*um;
  G4double window_gap       = 1.0*mm;

  // ---------------- set materials for the detectors --------------

  // (Element name, symbol, atomic number, atomic mass) (as floats)
  G4Element* Si = new G4Element("Silicon","Si", 14., 28.0855*g/mole); // main wafer material for detector
  G4Element* B = new G4Element("Boron","B", 5., 10.811*g/mole);   // possible doping material

  // Final doped silicon material to be used in the electron detector
  G4Material* DopedSilicon = new G4Material("DopedSilicon", 5.8*g/cm3, 2); // last argument is number of components in material
  DopedSilicon->AddElement(Si, 99.9*perCent);
  DopedSilicon->AddElement(B, 0.1*perCent);

  // ---------------- set detector positions --------------
  G4ThreeVector detector1_pos  = G4ThreeVector(0, 0, world_offset);

  G4VSolid* detector1_solid = new G4Box("detector",
                   0.5*detector_dimX, 0.5*detector_dimY, 0.5*detector1_thickness);

  G4VSolid* detector2_solid = new G4Box("detector",
                  0.5*detector_dimX, 0.5*detector_dimY, 0.5*detector2_thickness);

  // ---------------- create detector 1 --------------

  detector1 =
  new G4LogicalVolume(detector1_solid,      //its solid
                      DopedSilicon,        //its material
                      "detector1");      //its name

  new G4PVPlacement(0,                     //no rotation
                  detector1_pos,            //at position
                  detector1,                //its logical volume
                  "detector1",           //its name
                  logicEnv,                //its mother  volume
                  false,                   //no boolean operation
                  0,                       //copy number
                  checkOverlaps);          //overlaps checking

  // ---------------- create detector 2 --------------
  G4ThreeVector detector2_pos  = G4ThreeVector(0, 0, detector1_thickness/2 + detector2_thickness/2 + distance_between_detectors + world_offset);

  detector2 =
  new G4LogicalVolume(detector2_solid,      //its solid
                      DopedSilicon,        //its material
                      "detector2");      //its name

  new G4PVPlacement(0,                     //no rotation
                  detector2_pos,            //at position
                  detector2,                //its logical volume
                  "detector2",           //its name
                  logicEnv,                //its mother  volume
                  false,                   //no boolean operation
                  0,                       //copy number
                  checkOverlaps);          //overlaps checking
  // ---------------- set window material --------------
  // comment out the window for scattering distributions
  
  G4Material* window_material = nist->FindOrBuildMaterial("G4_Be");
  G4VSolid*   window_solid = new G4Box("window", 0.5*window_dimX, 0.5*window_dimY, 0.5*window_thickness);

  // ---------------- set window position --------------

  G4ThreeVector window_pos;
  window_pos = G4ThreeVector(0, 0, -(detector1_thickness/2 + window_thickness/2 + window_gap)+world_offset);

  // ---------------- create window --------------
  
  G4LogicalVolume* window =
  new G4LogicalVolume(window_solid,         // its solid
                      window_material,      // its material
                      "window");    // its name
  
  new G4PVPlacement(0,                       //no rotation
                    window_pos,              //at position
                    window,                  //its logical volume
                    "window",       //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  // ---------------- create coded aperture --------------
  /*
  // comment this out for the 2 detector configuration
  std::fstream configFile;
  configFile.open("../src/config.txt", std::ios_base::in);

  // Initialize variables
  G4double nelements, ca_thickness_um, ca_gap_cm, hole_size_mm;
  G4String filename;

  // Load into variables
  configFile >> nelements >> ca_thickness_um >> ca_gap_cm >> hole_size_mm >> filename;

  configFile.close();

  // mask size  
  //G4double nelements = 11;
  //G4String filename = "../src/mask_designs/11mosaicMURA_matrix.txt";
  // set coded aperture parameters

  G4double ca_thickness = ca_thickness_um*um;
  G4double ca_gap = ca_gap_cm*cm; 
  G4Material* ca_material = nist->FindOrBuildMaterial("G4_W");

  G4double ca_pos = -(detector1_thickness/2 + ca_gap - ca_thickness/2)+world_offset; // defined so that the gap is from the front of the first detector to front of mask
  G4double hole_size = hole_size_mm*mm; //0.649*mm; // same as element size
  G4double ca_size = (nelements * hole_size);
  
  // for mosiac
  //G4double ca_size = (nelements * hole_size * 2) - hole_size; 
  
  G4double mask_offset = -(ca_size/2 - hole_size/2); // centering to correct for origin

  // create the mask element -- pad with 0.5um to ensure overlap
  G4VSolid* ca_element = new G4Box("hole",
                  0.5*hole_size + 0.5*um, 0.5*hole_size + 0.5*um, 0.5*ca_thickness);
  
  G4LogicalVolume* logicMask =
  new G4LogicalVolume(ca_element,         // its solid
                      ca_material,      // its material
                      "logicMask");    // its name

  // start adding the blocks of tungsten
  // name the variables
  G4String placementXY_str; 
  G4double placementX, placementY; 
  G4String token;
  
  std::ifstream placementFile(filename, std::ios_base::in);
  
  // Get number of lines in file
  unsigned int numberOfBoxes = 0;
  while(getline(placementFile, placementXY_str, '\n'))
    { numberOfBoxes++; }
  
  placementFile.close();

  // Reopen file to start from first line
  placementFile.open(filename, std::ios_base::in);
  getline(placementFile, placementXY_str, '\n');
  
  token = placementXY_str.substr(
		  0, 
  		  placementXY_str.find(',')); 
  
  placementX = std::stod(token);
  
  token = placementXY_str.substr(
		  placementXY_str.find(',')+1, 
		  placementXY_str.find('\n'));
  
  placementY = std::stod(token);

  // save first element -- convert from mm to cm
  G4double first_hole_x = 0.1*placementX*cm;
  G4double first_hole_y = 0.1*placementY*cm;

  // place first element
  G4VPhysicalVolume *physMask = new G4PVPlacement(0,  
      G4ThreeVector(mask_offset + first_hole_x, -(mask_offset + first_hole_y), ca_pos), logicMask, "physMask", logicEnv, false, 0, checkOverlaps);

  // start looping through the file ----------
  // starts at 1 since logicAp1 uses first line of file 
  for(unsigned int i=1; i<numberOfBoxes; i++)
  {
    getline(placementFile, placementXY_str, '\n');

    token = placementXY_str.substr(
		  0, 
  		  placementXY_str.find(',')); 
    
    placementX = std::stod(token); 
 
    token = placementXY_str.substr(
		  placementXY_str.find(',')+1, 
		  placementXY_str.find('\n'));
    
    placementY = std::stod(token); 

    // place the elements
    // for some dumb reason i flipped the y axis...
    G4VPhysicalVolume *physMask = new G4PVPlacement(0,  
         G4ThreeVector(mask_offset + (0.1*placementX*cm), -(mask_offset + (0.1*placementY*cm)), ca_pos), 
         logicMask, "physmask", logicEnv, false, 0, checkOverlaps);
  }

  placementFile.close();

  // finally add an outline of the 2*hole size around it
  G4VSolid* xoutline = new G4Box("hole",
                  0.5*(ca_size+2*hole_size),0.5*hole_size,  0.5*ca_thickness);
  G4VSolid* youtline = new G4Box("hole",
                   0.5*hole_size, 0.5*(ca_size+2*hole_size), 0.5*ca_thickness);
  
  G4LogicalVolume* xlogicOutline =
  new G4LogicalVolume(xoutline,         // its solid
                      ca_material,      // its material
                      "xlogicOutline");    // its name
  G4LogicalVolume* ylogicOutline =
  new G4LogicalVolume(youtline,         // its solid
                      ca_material,      // its material
                      "ylogicOutline");    // its name

  // place the outline
  G4VPhysicalVolume *xphysOutline1 = new G4PVPlacement(0,  
      G4ThreeVector(0, -0.5*(ca_size +2*hole_size - hole_size),ca_pos), 
      xlogicOutline, "xphysOutline", logicEnv, false, 0, checkOverlaps);
  // place the outline
  G4VPhysicalVolume *xphysOutline2 = new G4PVPlacement(0,  
      G4ThreeVector(0, 0.5*(ca_size +2*hole_size - hole_size),ca_pos), 
      xlogicOutline, "xphysOutline", logicEnv, false, 0, checkOverlaps);
  // place the outline
  G4VPhysicalVolume *yphysOutline1 = new G4PVPlacement(0,  
      G4ThreeVector(-0.5*(ca_size +2*hole_size - hole_size),0, ca_pos), 
      ylogicOutline, "yphysOutline", logicEnv, false, 0, checkOverlaps);
  // place the outline
  G4VPhysicalVolume *yphysOutline2 = new G4PVPlacement(0,  
      G4ThreeVector(0.5*(ca_size +2*hole_size - hole_size), 0,ca_pos), 
      ylogicOutline, "yphysOutline", logicEnv, false, 0, checkOverlaps);
  // always return the physical World
  */

  return physWorld;
}

void DetectorConstruction::ConstructSDandField(){
  G4String SDname;
  // register the SD 
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  MySensitiveDetector* sd1 = new MySensitiveDetector(SDname = "/SD1");
  SDman->AddNewDetector(sd1);
  detector1->SetSensitiveDetector(sd1);

  MySensitiveDetector* sd2 = new MySensitiveDetector(SDname = "/SD2");
  SDman->AddNewDetector(sd2);
  detector2->SetSensitiveDetector(sd2);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
