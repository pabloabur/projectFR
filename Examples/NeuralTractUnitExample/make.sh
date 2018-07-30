gfortran -c ../../PointProcessGenerator.f90
gfortran -c ../../ogpf.f90
gfortran -c ../../NeuralTractUnit.f90
gfortran -c ../../DynamicalArrays.f90
gfortran -c ../../randomSeedInitialize.f90

gfortran NeuralTractUnitExample.f90 -o NeuralTractUnitExample ../../NeuralTractUnit.f90 ../../PointProcessGenerator.f90 ../../DynamicalArrays.f90 ../../randomSeedInitialize.f90 ../../ogpf.f90 


