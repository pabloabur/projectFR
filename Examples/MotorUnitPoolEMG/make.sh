gfortran-7 -c ../../ogpf.f90
gfortran-7 -c ../../DynamicalArrays.f90
gfortran-7 -c ../../randomSeedInitialize.f90
gfortran-7 -c ../../Configuration.f90
gfortran-7 -c ../../PulseConductanceState.f90
gfortran-7 -c ../../ChannelConductance.f90
gfortran-7 -c ../../Compartment.f90
gfortran-7 -c ../../AxonDelay.f90
gfortran-7 -c ../../MotorUnit.f90
gfortran-7 -c ../../MuscularActivation.f90
gfortran-7 -c ../../MotorUnitPool.f90
gfortran-7 -c ../../MuscleNoHill.f90

gfortran-7 MotorUnitPoolEMG.f90 -o MotorUnitPoolEMG ../../MuscleNoHill.f90 ../../MuscularActivation.f90 ../../MotorUnitPool.f90 ../../AxonDelay.f90 ../../MotorUnit.f90   ../../DynamicalArrays.f90 ../../randomSeedInitialize.f90 ../../ogpf.f90 ../../Configuration.f90 ../../PulseConductanceState.f90 ../../ChannelConductance.f90 ../../Compartment.f90


