gfortran -c ../../PointProcessGenerator.f90
gfortran -c ../../ogpf.f90
gfortran -c ../../NeuralTractUnit.f90
gfortran -c ../../DynamicalArrays.f90
gfortran -c ../../randomSeedInitialize.f90
gfortran -c ../../Configuration.f90

gfortran NeuralTractExample.f90 -o NeuralTractExample ../../NeuralTractUnit.f90 ../../PointProcessGenerator.f90 ../../DynamicalArrays.f90 ../../randomSeedInitialize.f90 ../../ogpf.f90 ../../Configuration.f90 


