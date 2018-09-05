gfortran -c ../../queue.f90
gfortran -c ../../ogpf.f90
gfortran -c ../../DynamicalArrays.f90
gfortran -c ../../String.f90
gfortran -c ../../CharacterArray.f90
gfortran -c ../../CharacterMatrix.f90
gfortran -c ../../randomSeedInitialize.f90
gfortran -c ../../Configuration.f90
gfortran -c ../../PulseConductanceState.f90
gfortran -c ../../ChannelConductance.f90
gfortran -c ../../Synapse.f90
gfortran -c ../../SynapsePointer.f90
gfortran -c ../../Compartment.f90
gfortran -c ../../AxonDelay.f90
gfortran -c ../../MotorUnit.f90
gfortran -c ~/intel/mkl/include/mkl_spblas.f90
gfortran -c ../../MuscularActivation.f90
gfortran -c ../../MuscleNoHill.f90
gfortran -c ../../MuscleHill.f90
gfortran -c ../../MuscleSpindle.f90
gfortran -c ../../MotorUnitPool.f90
gfortran -c ../../PointProcessGenerator.f90
gfortran -c ../../NeuralTractUnit.f90
gfortran -c ../../NeuralTract.f90
gfortran -c ../../Interneuron.f90
gfortran -c ../../InterneuronPool.f90
gfortran -c ../../SynapticNoise.f90
gfortran -c ../../AfferentUnit.f90
gfortran -c ../../AfferentPool.f90
gfortran -c ../../SynapsesFactory.f90
gfortran -c ../../MusclePointer.f90
gfortran -c ../../jointAnkleForceTask.f90


gfortran -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -finit-real=nan -std=f2003 -O0 FrequencyAnalysis.f90 -o FrequencyAnalysis -O3 MusclePointer.o jointAnkleForceTask.o AfferentUnit.o AfferentPool.o String.o SynapticNoise.o MuscleSpindle.o MuscleHill.o Interneuron.o InterneuronPool.o SynapsePointer.o SynapsesFactory.o Synapse.o CharacterMatrix.o PointProcessGenerator.o NeuralTractUnit.o NeuralTract.o  queue.o CharacterArray.o MuscleNoHill.o MuscularActivation.o MotorUnitPool.o AxonDelay.o Compartment.o MotorUnit.o   DynamicalArrays.o randomSeedInitialize.o ogpf.o Configuration.o PulseConductanceState.o ChannelConductance.o -Wl,--start-group ~/intel/mkl/lib/intel64/libmkl_gf_lp64.a ~/intel/mkl/lib/intel64/libmkl_sequential.a ~/intel/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl -m64 -I ~/intel/mkl/include
#gfortran StaticDynamicProperties.f90 -o StaticDynamicProperties -O3 ../../MusclePointer.f90 ../../jointAnkleForceTask.f90 ../../AfferentUnit.f90 ../../AfferentPool.f90 ../../String.f90 ../../SynapticNoise.f90 ../../MuscleSpindle.f90 ../../MuscleHill.f90 ../../Interneuron.f90 ../../InterneuronPool.f90 ../../SynapsePointer.f90 ../../SynapsesFactory.f90 ../../Synapse.f90 ../../CharacterMatrix.f90 ../../PointProcessGenerator.f90 ../../NeuralTractUnit.f90 ../../NeuralTract.f90  ../../queue.f90 ../../CharacterArray.f90 ../../MuscleNoHill.f90 ../../MuscularActivation.f90 ../../MotorUnitPool.f90 ../../AxonDelay.f90 ../../Compartment.f90 ../../MotorUnit.f90   ../../DynamicalArrays.f90 ../../randomSeedInitialize.f90 ../../ogpf.f90 ../../Configuration.f90 ../../PulseConductanceState.f90 ../../ChannelConductance.f90 -Wl,--start-group ~/intel/mkl/lib/intel64/libmkl_gf_lp64.a ~/intel/mkl/lib/intel64/libmkl_sequential.a ~/intel/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl -m64 -I ~/intel/mkl/include

