# Whole Body Control of Krang Simulation

## Dependencies

- DART
 [Dart Homepage](https://dartsim.github.io)

- config4cpp
 [Source Code Download](http://www.config4star.org/#main-source-code)

  How to add config4cpp to your system (Linux/BSD/OSX)

  1: Download and extract the source code of config4cpp

  2: cd config4cpp

  3: Follow the README.txt in config4cpp/

  4: Run the following commands to add the make'd files into your local system

        sudo cp bin/{config2cpp,config4cpp} /usr/local/bin &&
        sudo cp lib/libconfig4cpp.a /usr/local/lib &&
        sudo cp -r include/config4cpp /usr/local/include &&
        sudo chmod g+rx /usr/local/include/config4cpp &&
        sudo chmod o+rx /usr/local/include/config4cpp

  \*Note: You can just copy paste the above block of commands

- [09-URDF](https://github.gatech.edu/WholeBodyControlAttempt1/09-URDF)
 Install the repository.

## Use
In `examples/3dofddp/controlParams.cfg`, there are many configuration parameters that the user can play with. For now, their meaning will need to be understood by looking at their role in the source code itself.

To compile:

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make

To run:

    cd examples/3dofddp
    ./3dofddp
