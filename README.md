#This simple program is provided as a reproductivity sample of the numerical methods. It is not fully encapsulated and does not make use of advanced C++ features.

complie command:    bash build.sh
run command:       ./build/exe
Results are written to the 'result' folder. Example: ./result/SodComp300.dat

Simulation parameters are in the `in.in` file.  The sample program will automatically read this file at runtime.
'in.in'：six items.
#-------------------
   1.whatflow Lax           #Problem name (./src/InitialCDs.cpp): please choose one from the list below
                            #Lax
                            #Sod
                            #ShuOsher
                            #Blast
                            #LeBlanc

   2.epsilon 1e-6           #used in the WENOJS scheme(./src/WENO5JS)

   3.Time 0                 # Time = 0 → the integration time will use the program's default setting (./src/InitialCDs.cpp).  
                            # Time ≠ 0 → the integration time will advance according to the user-defined value.

   4.nx 200                 # Grid size

   5.cfl 0.5

   6.ReconstructType Comp   #   Reconstruction options (./src/FourMethodsForTests.cpp), please choose one:
                            #   Comp    → component-wise reconstruction from the Global Lax-Friedrichs split fluxes(./src/GLF_FVS.cpp)
                            #   Char    → characteristic-wise reconstruction
                            #   PS      → Pressure and entropy replacement treatment
                            #   Co      → improved common-weights method
#-------------------

