# 09: URDF
This repository contains all the URDF files required for modeling hardware.

## Dependencies

- DART
 [Dart Homepage](https://dartsim.github.io)

## Build and Run

 1: Enter the project directory

 2: Build the project

	mkdir build
    cd build
    cmake ..
    make

 3: Run the project

    ./{project_name}

   Press [Space] to start Dart simulation from pop-up window

## Installation

 1: Enter the project directory

 2: Install the project

    mkdir build
    cd build
    cmake ..
    sudo make install

## Uninstall

 To remove system files created by the installation of this repo.

    sudo make uninstall

## 3DOF-WIP
Contains files for a 3DOF wheel inverted pendulum model.

## 7DOFArm
Contains files for a single 7DOF Arm model.

## Krang
Contains files for full body Krang model.

## KrangFixedWheels
Contains files for fixed wheels Krang model.

## KrangWaist
Contains files for fixed base Krang model.

## scenes
Contains files for different object and scene models.
