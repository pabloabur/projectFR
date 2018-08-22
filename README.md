# projectFR

This is a Neuromusculoskletal simulator written in Fortran. 

This software uses MKL Intel Library, avaliable at [https://software.intel.com/en-us/mkl]. After installing it on your computer run the following command at a terminal:

    source $MKLDIR/bin/mklvars.sh [intel64/ia32]

where $MKLDIR is the directory where the MKL library was installed. Use the option intel64 for 64 bits archutecture and ia32 for 32 bits architecture. For example, on my computer the  command is: source /opt/intel/mkl/bin/mklvars.sh intel64 . Almost all of the recent (10 years old) computers are 64 bits.

After that, you can run any of the examples in the Examples directory. Each example is contained in one directory, where you can find a file make.sh and a file makewithintel.sh . The former compiles the program using the GNU Fortran 8.1 compiler and the latter compiles the program using the ifort (Intel Fortran Compiler).  Currently, the ifort compiler is not compiling this version of the program.

