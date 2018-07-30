ifort -c ../../ogpf.f90
ifort -c ../../DynamicalArrays.f90
ifort -c ../../randomSeedInitialize.f90
ifort -c ../../Configuration.f90
ifort -c ../../PulseConductanceState.f90
ifort -c ../../ChannelConductance.f90
ifort -c ../../Compartment.f90
ifort -c ../../AxonDelay.f90
ifort -c ../../MotorUnit.f90
ifort -c ../../MuscularActivation.f90
ifort -c ../../MotorUnitPool.f90
ifort -c ../../MuscleNoHill.f90

ifort MotorUnitPoolEMG.f90 -o MotorUnitPoolEMG ../../MuscularActivation.f90 ../../MuscleNoHill.f90  ../../MotorUnitPool.f90 ../../AxonDelay.f90 ../../MotorUnit.f90   ../../DynamicalArrays.f90 ../../randomSeedInitialize.f90 ../../ogpf.f90 ../../Configuration.f90 ../../PulseConductanceState.f90 ../../ChannelConductance.f90 ../../Compartment.f90



