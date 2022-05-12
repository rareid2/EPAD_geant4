<h3 align="center">EPAD</h3>

<h3 align="center"> Electron Pitch Angle Detector -- instrument name is a work in progress </h3>

## About
GEANT4 is a toolkit for simulating the passage of radiation through matter. GEANT4 is required to compile this repo. 
Follow this [link](https://geant4.web.cern.ch/support/download) to install and build GEANT4. 

The purpose of this repo is to simulate an energetic electron detector for pitch angle measurements in the radiation belts. 

In this repo, you'll find:
- code to build a solid state silicon detector
- code to build a coded aperture with various designs into GEANT4, included Modified Uniformly Redundant Arrays (MURA)
- GEANT4 macro files to simulate radiation belt fluxes with various pitch angle distributions 

## Getting Started

Once you have GEANT4 installed and built succesffully, clone this repo and navigate to the main directory. 

Run:

`mkdir build` 

`cd build`

`cmake ..`

Finally, build with the following command, replace N with the number of cores you'd like to compile with: 

`make -j N `

## Runnning

To run in interactive mode, run the following line from the main directory: 

`build/main`

To try running a pitch angle distribution in the visualization environment, run the following in the command line in the visualization environment:

`Idle> /control/execute macros/run_example.mac`
