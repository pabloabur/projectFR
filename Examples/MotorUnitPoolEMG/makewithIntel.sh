ifort -c ../../queue.f90
ifort -c ../../ogpf.f90
ifort -c ../../DynamicalArrays.f90
ifort -c ../../CharacterArray.f90
ifort -c ../../CharacterMatrix.f90
ifort -c ../../randomSeedInitialize.f90
ifort -c ../../Configuration.f90
ifort -c ../../PulseConductanceState.f90
ifort -c ../../ChannelConductance.f90
ifort -c ../../Synapse.f90
ifort -c ../../SynapsePointer.f90
ifort -c ../../Compartment.f90
ifort -c ../../AxonDelay.f90
ifort -c ../../MotorUnit.f90
ifort -c ../../MuscularActivation.f90
ifort -c ../../MotorUnitPool.f90
ifort -c ../../PointProcessGenerator.f90
ifort -c ../../NeuralTractUnit.f90
ifort -c ../../NeuralTract.f90
ifort -c ../../MuscleNoHill.f90
ifort -c ../../SynapsesFactory.f90


ifort MotorUnitPoolEMG.f90 -o MotorUnitPoolEMG  ../../SynapsePointer.f90 ../../SynapsesFactory.f90 ../../Synapse.f90 ../../CharacterMatrix.f90 ../../PointProcessGenerator.f90 ../../NeuralTractUnit.f90 ../../NeuralTract.f90  ../../queue.f90 ../../CharacterArray.f90 ../../MuscleNoHill.f90 ../../MuscularActivation.f90 ../../MotorUnitPool.f90 ../../AxonDelay.f90 ../../Compartment.f90 ../../MotorUnit.f90   ../../DynamicalArrays.f90 ../../randomSeedInitialize.f90 ../../ogpf.f90 ../../Configuration.f90 ../../PulseConductanceState.f90 ../../ChannelConductance.f90 