# $Id: CMakeLists.txt 86065 2014-11-07 08:51:15Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(detector)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#

option(CONSTRUCT "which construction to use - coded aperture, two detectors, pinhole" CA) #coded aperture by default
option(PARTICLE_SOURCE "configure with particle gun (pt src) or distribution" DS) #distribution by default 

if((${CONSTRUCT} STREQUAL CA) AND (${PARTICLE_SOURCE} STREQUAL DS)) # default state

  file(GLOB sources ${PROJECT_SOURCE_DIR}/src/ActionInitialization_DS.cc ${PROJECT_SOURCE_DIR}/src/EventAction_CA.cc ${PROJECT_SOURCE_DIR}/src/DetectorConstruction_CA.cc ${PROJECT_SOURCE_DIR}/src/Hit.cc ${PROJECT_SOURCE_DIR}/src/PrimaryGeneratorAction_DS.cc ${PROJECT_SOURCE_DIR}/src/PrimaryGeneratorMessenger.cc ${PROJECT_SOURCE_DIR}/src/RunAction_DS.cc ${PROJECT_SOURCE_DIR}/src/SensitiveDetector_CA.cc ${PROJECT_SOURCE_DIR}/src/SteppingAction.cc)  
  file(GLOB headers ${PROJECT_SOURCE_DIR}/include/ActionInitialization.hh ${PROJECT_SOURCE_DIR}/include/DetectorConstruction.hh ${PROJECT_SOURCE_DIR}/include/EventAction.hh ${PROJECT_SOURCE_DIR}/include/Hit.hh ${PROJECT_SOURCE_DIR}/include/PrimaryGeneratorAction_DS.hh ${PROJECT_SOURCE_DIR}/include/PrimaryGeneratorMessenger.hh ${PROJECT_SOURCE_DIR}/include/RunAction.hh ${PROJECT_SOURCE_DIR}/include/SensitiveDetector.hh ${PROJECT_SOURCE_DIR}/include/SteppingAction.hh)

elseif((${CONSTRUCT} STREQUAL CA) AND (${PARTICLE_SOURCE} STREQUAL PS))

  file(GLOB sources ${PROJECT_SOURCE_DIR}/src/ActionInitialization_PS.cc ${PROJECT_SOURCE_DIR}/src/EventAction_CA.cc ${PROJECT_SOURCE_DIR}/src/DetectorConstruction_CA.cc ${PROJECT_SOURCE_DIR}/src/Hit.cc ${PROJECT_SOURCE_DIR}/src/PrimaryGeneratorAction_PS.cc ${PROJECT_SOURCE_DIR}/src/RunAction_PS.cc ${PROJECT_SOURCE_DIR}/src/SensitiveDetector_CA.cc ${PROJECT_SOURCE_DIR}/src/SteppingAction.cc)  
  file(GLOB headers ${PROJECT_SOURCE_DIR}/include/ActionInitialization.hh ${PROJECT_SOURCE_DIR}/include/DetectorConstruction.hh ${PROJECT_SOURCE_DIR}/include/EventAction.hh ${PROJECT_SOURCE_DIR}/include/Hit.hh ${PROJECT_SOURCE_DIR}/include/PrimaryGeneratorAction_PS.hh ${PROJECT_SOURCE_DIR}/include/RunAction.hh ${PROJECT_SOURCE_DIR}/include/SensitiveDetector.hh ${PROJECT_SOURCE_DIR}/include/SteppingAction.hh)

elseif((${CONSTRUCT} STREQUAL TD) AND (${PARTICLE_SOURCE} STREQUAL DS))

  file(GLOB sources ${PROJECT_SOURCE_DIR}/src/ActionInitialization_DS.cc ${PROJECT_SOURCE_DIR}/src/EventAction_TD.cc ${PROJECT_SOURCE_DIR}/src/DetectorConstruction_TD.cc ${PROJECT_SOURCE_DIR}/src/Hit.cc ${PROJECT_SOURCE_DIR}/src/PrimaryGeneratorAction_DS.cc ${PROJECT_SOURCE_DIR}/src/PrimaryGeneratorMessenger.cc ${PROJECT_SOURCE_DIR}/src/RunAction_DS.cc ${PROJECT_SOURCE_DIR}/src/SensitiveDetector_TD.cc ${PROJECT_SOURCE_DIR}/src/SteppingAction.cc)  
  file(GLOB headers ${PROJECT_SOURCE_DIR}/include/ActionInitialization.hh ${PROJECT_SOURCE_DIR}/include/DetectorConstruction.hh ${PROJECT_SOURCE_DIR}/include/EventAction.hh ${PROJECT_SOURCE_DIR}/include/Hit.hh ${PROJECT_SOURCE_DIR}/include/PrimaryGeneratorAction_DS.hh ${PROJECT_SOURCE_DIR}/include/PrimaryGeneratorMessenger.hh ${PROJECT_SOURCE_DIR}/include/RunAction.hh ${PROJECT_SOURCE_DIR}/include/SensitiveDetector.hh ${PROJECT_SOURCE_DIR}/include/SteppingAction.hh)

elseif((${CONSTRUCT} STREQUAL TD) AND (${PARTICLE_SOURCE} STREQUAL PS))

  file(GLOB sources ${PROJECT_SOURCE_DIR}/src/ActionInitialization_PS.cc ${PROJECT_SOURCE_DIR}/src/EventAction_TD.cc ${PROJECT_SOURCE_DIR}/src/DetectorConstruction_TD.cc ${PROJECT_SOURCE_DIR}/src/Hit.cc ${PROJECT_SOURCE_DIR}/src/PrimaryGeneratorAction_PS.cc ${PROJECT_SOURCE_DIR}/src/RunAction_PS.cc ${PROJECT_SOURCE_DIR}/src/SensitiveDetector_TD.cc ${PROJECT_SOURCE_DIR}/src/SteppingAction.cc)  
  file(GLOB headers ${PROJECT_SOURCE_DIR}/include/ActionInitialization.hh ${PROJECT_SOURCE_DIR}/include/DetectorConstruction.hh ${PROJECT_SOURCE_DIR}/include/EventAction.hh ${PROJECT_SOURCE_DIR}/include/Hit.hh ${PROJECT_SOURCE_DIR}/include/PrimaryGeneratorAction_PS.hh ${PROJECT_SOURCE_DIR}/include/RunAction.hh ${PROJECT_SOURCE_DIR}/include/SensitiveDetector.hh ${PROJECT_SOURCE_DIR}/include/SteppingAction.hh)

else() # default state

  file(GLOB sources ${PROJECT_SOURCE_DIR}/src/ActionInitialization_DS.cc ${PROJECT_SOURCE_DIR}/src/EventAction_CA.cc ${PROJECT_SOURCE_DIR}/src/DetectorConstruction_CA.cc ${PROJECT_SOURCE_DIR}/src/Hit.cc ${PROJECT_SOURCE_DIR}/src/PrimaryGeneratorAction_DS.cc ${PROJECT_SOURCE_DIR}/src/PrimaryGeneratorMessenger.cc ${PROJECT_SOURCE_DIR}/src/RunAction_DS.cc ${PROJECT_SOURCE_DIR}/src/SensitiveDetector_CA.cc ${PROJECT_SOURCE_DIR}/src/SteppingAction.cc)  
  file(GLOB headers ${PROJECT_SOURCE_DIR}/include/ActionInitialization.hh ${PROJECT_SOURCE_DIR}/include/DetectorConstruction.hh ${PROJECT_SOURCE_DIR}/include/EventAction.hh ${PROJECT_SOURCE_DIR}/include/Hit.hh ${PROJECT_SOURCE_DIR}/include/PrimaryGeneratorAction_DS.hh ${PROJECT_SOURCE_DIR}/include/PrimaryGeneratorMessenger.hh ${PROJECT_SOURCE_DIR}/include/RunAction.hh ${PROJECT_SOURCE_DIR}/include/SensitiveDetector.hh ${PROJECT_SOURCE_DIR}/include/SteppingAction.hh)

endif()
unset(CONSTRUCT CACHE) # <---- this is the important!!
unset(PARTICLE_SOURCE CACHE) # <---- this is the important!!

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(main electron_detector_main.cc ${sources} ${headers})
target_link_libraries(main ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
# Grant N.B.: macros are generally run as follows: ./main ../macros/run___.mac
# so not all macro run files need to be copied
set(EXAMPLEB1_SCRIPTS
  init_vis.mac
  #macros/run_y_axis.mac
  #macros/run_adjoint_simulation_electron.mac
  vis.mac
  )

foreach(_script ${EXAMPLEB1_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(B1 DEPENDS main)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS main DESTINATION bin)


