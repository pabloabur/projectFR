program ForceVariability
    use MotorUnitPoolClass
    use NeuralTractClass
    use ConfigurationClass
    use ogpf 
    use randomSeedInitialize
    use MuscleNoHillClass
    use CharacterArrayClass
    use CharacterMatrixClass
    use QueueClass
    use SynapsesFactoryModule
    implicit none 
    !integer, parameter :: wp = kind(1.0d0)
    type(Configuration) :: conf
    real(wp), parameter :: PI = 4 * atan(1.0_wp)    
    real(wp) :: dt
    real(wp) :: tf
    integer :: timeLength
    integer :: i, j
    real(wp), dimension(:), allocatable :: t, MNv_mV
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp) :: FR
    integer :: GammaOrder 
    character(len = 80) :: pool = 'CMExt'
    character(len = 80) :: filename = 'confMNPoolWithDescendingCommand.rmto'
    type(CharacterArray) :: dynamics  
    Character(len = 80) :: synapseDyn
    type(Queue) :: spikesQueue
    integer :: newItem
    logical, dimension(10) :: logList
    integer, dimension(10) :: logCount
    type(CharacterMatrix) :: Synapses
    type(MotorUnitPool), dimension(1), target :: motorUnitPools
    type(NeuralTract), dimension(1) :: neuralTractPools    
    real(wp) :: alpha, beta
    integer, parameter ::  M = 1600, K = 1600, N = 1
    real(wp) :: A(M,K), B(K,N), C(M,N)

    print '(A)', "This example computes real matrix C=alpha*A*B+beta*C"
    print '(A)', "using Intel(R) MKL function dgemm, where A, B, and C"
    print '(A)', "are matrices and alpha and beta are double precision "
    print '(A)', "scalars"
    print '(A)', ""

    print '(A)', "Initializing data for matrix multiplication C=A*B for "
    print '(a,I5,a,I5,a,I5,a,I5,a)', " matrix A(",M," x",K, ") and matrix B(", K," x", N, ")"
    print *, ""
    alpha = 1.0 
    beta = 0.0
    call init_random_seed()
    print *, "Intializing matrix data"
    print *, ""
    call RANDOM_NUMBER(A)

    call RANDOM_NUMBER(B)

    C(:,:) = 0.0

    print '(A)', "Computing matrix product using Intel(R) MKL DGEMM "
    print '(A)', "subroutine"
    call cpu_time(tic)
        do i = 1, 300
            call dgemv('N',M,K,alpha,A,M,B,1,beta,C,1)
        end do
    call cpu_time(toc)
    print '(F15.6, A)', (toc - tic), ' seconds'
    print '(A)', "Computations completed."
    print '(A)', ""

    print '(A)', "Top left corner of matrix A:"
    print '(ES12.4,1x)', ((A(I,J), J = 1,MIN(K,6)), I = 1,MIN(M,6))
    print *, ""

    print '(A)', "Top left corner of matrix B:"
    print '(F12.0,1x)', ((B(i,j),j = 1,MIN(N,6)), i = 1,MIN(K,6))
    print '(A)', ""



    print '(A)', "Top left corner of matrix C:"
    print '(ES12.4,1x)', ((C(i,j), j = 1,MIN(N,1)), i = 1,MIN(M,1))
    print *, ""

    call cpu_time(tic)
    do i = 1, 300
        C = matmul(A,B)
    end do
    call cpu_time(toc)
    print '(F15.6, A)', (toc - tic), ' seconds'

    print '(A)', "Top left corner of matrix C:"
    print '(ES12.4,1x)', ((C(i,j), j = 1,MIN(N,1)), i = 1,MIN(M,1))
    print *, ""

    print '(A)', "Example completed."
    
    

    

    
    
    

    
    
end program ForceVariability