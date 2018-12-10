gfortran -c -g ../../queue.f90
gfortran -c -g ../../ogpf.f90
gfortran -c -g ../../DynamicalArrays.f90
gfortran -c -g ../../String.f90
gfortran -c -g ../../CharacterArray.f90
gfortran -c -g ../../CharacterMatrix.f90
gfortran -c -g ../../randomSeedInitialize.f90
gfortran -c -g ../../Configuration.f90
gfortran -c -g ../../PulseConductanceState.f90
gfortran -c -g ../../ChannelConductance.f90
gfortran -c -g ../../Synapse.f90
gfortran -c -g ../../SynapsePointer.f90
gfortran -c -g ../../Compartment.f90
gfortran -c -g ../../AxonDelay.f90
gfortran -c -g ../../MotorUnit.f90
gfortran -c -g ~/intel/mkl/include/mkl_spblas.f90
gfortran -c -g ../../MuscularActivation.f90
gfortran -c -g ../../MuscleNoHill.f90
gfortran -c -g ../../MuscleHill.f90
gfortran -c -g ../../MuscleSpindle.f90
gfortran -c -g ../../GolgiTendonOrgan.f90
gfortran -c -g ../../MotorUnitPool.f90
gfortran -c -g ../../PointProcessGenerator.f90
gfortran -c -g ../../NeuralTractUnit.f90
gfortran -c -g ../../NeuralTract.f90
gfortran -c -g ../../Interneuron.f90
gfortran -c -g ../../InterneuronPool.f90
gfortran -c -g ../../SynapticNoise.f90
gfortran -c -g ../../AfferentUnit.f90
gfortran -c -g ../../AfferentPool.f90
gfortran -c -g ../../SynapsesFactory.f90
gfortran -c -g ../../MusclePointer.f90
gfortran -c -g ../../jointAnkleForceTask.f90


gfortran AntidromicStimulationofMNandRC.f90 -o AntidromicStimulationofMNandRC -O0 ../../GolgiTendonOrgan.f90 ../../MusclePointer.f90 ../../jointAnkleForceTask.f90 ../../AfferentUnit.f90 ../../AfferentPool.f90 ../../String.f90 ../../SynapticNoise.f90 ../../MuscleSpindle.f90 ../../MuscleHill.f90 ../../Interneuron.f90 ../../InterneuronPool.f90 ../../SynapsePointer.f90 ../../SynapsesFactory.f90 ../../Synapse.f90 ../../CharacterMatrix.f90 ../../PointProcessGenerator.f90 ../../NeuralTractUnit.f90 ../../NeuralTract.f90  ../../queue.f90 ../../CharacterArray.f90 ../../MuscleNoHill.f90 ../../MuscularActivation.f90 ../../MotorUnitPool.f90 ../../AxonDelay.f90 ../../Compartment.f90 ../../MotorUnit.f90   ../../DynamicalArrays.f90 ../../randomSeedInitialize.f90 ../../ogpf.f90 ../../Configuration.f90 ../../PulseConductanceState.f90 ../../ChannelConductance.f90    -Wl,--start-group ~/intel/mkl/lib/intel64/libmkl_gf_lp64.a ~/intel/mkl/lib/intel64/libmkl_sequential.a ~/intel/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl -m64 -I ~/intel/mkl/include
