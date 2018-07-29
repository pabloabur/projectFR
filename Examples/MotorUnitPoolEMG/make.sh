gfortran -c ../../ogpf.f90
gfortran -c ../../DynamicalArrays.f90
gfortran -c ../../randomSeedInitialize.f90
gfortran -c ../../Configuration.f90
gfortran -c ../../PulseConductanceState.f90
gfortran -c ../../ChannelConductance.f90
gfortran -c ../../Compartment.f90
gfortran -c ../../MotorUnit.f90
gfortran -c ../../AxonDelay.f90
gfortran -c ../../MotorUnitPool.f90

gfortran MotorUnitPoolEMG.f90 -o MotorUnitPoolEMG ../../MotorUnitPool.f90 ../../AxonDelay.f90 ../../MotorUnit.f90   ../../DynamicalArrays.f90 ../../randomSeedInitialize.f90 ../../ogpf.f90 ../../Configuration.f90 ../../PulseConductanceState.f90 ../../ChannelConductance.f90 ../../Compartment.f90


