program NeuralTractUnitExample
    use NeuralTractUnitClass
    use ogpf 
    use randomSeedInitialize
    implicit none 
    !integer, parameter :: wp = kind(1.0d0)
    type(NeuralTractUnit) :: NT_unit
    real(wp), parameter :: ISI_ms = 4.0
    real(wp), parameter :: FR_Hz = 1000/ISI_ms
    integer, parameter :: GammaOrder = 1
    real(wp), parameter :: PI = 4 * atan(1.0_wp)    
    !use, intrinsic :: iso_fortran_env, only : wp => real64
    real(wp), parameter :: dt = 0.05
    real(wp), parameter :: tf = 10000
    integer, parameter :: timeLength = int(tf/dt)
    integer :: i
    real(wp), dimension(timeLength) :: t
    real(wp) :: tic, toc
    character(len = 6) :: pool = 'NT_ext'
    integer :: index = 1
    type(gpf) :: gp
    real(wp), dimension(:), allocatable :: indices
    real(wp), dimension(:), allocatable :: ISI
    real(wp), dimension(:), allocatable :: ISIhist, ISIhistLimits
    real(wp) :: histogramMax
    real(wp) :: M, SD, CV, MFR, SD_FR

    t = [(dt*i, i=1, timeLength)]

    
    NT_unit = NeuralTractUnit(pool, index)
    
    
    call init_random_seed()


    call cpu_time(tic)
    do i = 1, timeLength
        call NT_unit%atualizeNeuralTractUnit(t(i), FR_Hz*dt/1000.0, GammaOrder)
    enddo
    call cpu_time(toc)
    print "(E15.6, A)", toc - tic,  " seconds"

    

    allocate(indices(size(NT_unit%spikesGenerator%points)))
    
    do i = 1, size(NT_unit%spikesGenerator%points)
        indices(i) = index    
    enddo
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call gp%title('Neural tract spikes instants')
    call gp%xlabel('t (s))')
    call gp%ylabel('Descending command index')
    call gp%plot(NT_unit%spikesGenerator%points/1000, indices, 'with points pt 5 lc rgb "#0008B0"')
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(ISI(size(NT_unit%spikesGenerator%points)-1))
    do i = 1, size(ISI)
        ISI(i) = NT_unit%spikesGenerator%points(i+1) - NT_unit%spikesGenerator%points(i) 
    end do

    allocate(ISIhist(10))
    allocate(ISIhistLimits(10))
    histogramMax = maxval(ISI)
    
    do i = 1, size(ISIhist)
        ISIhist(i) = count(ISI.gt.((i-1)*histogramMax/10.0).and.ISI.lt.(i*histogramMax/10.0))
        ISIhistLimits(i) = (i)*histogramMax/10.0
    end do

    call gp%title('ISI histogram')
    call gp%xlabel('ISI')
    call gp%ylabel('count')

    call gp%plot(ISIhistLimits,ISIhist, 'with impulses lw 2.5', &
            x2=ISIhistLimits, y2=ISIhist,  ls2='with points pt 7')
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    M = sum(ISI)/size(ISI)
    SD = SQRT(sum((ISI-M)**2)/size(ISI))
    CV = SD / M
    MFR = 1000/M
    
    print "(A, F15.6, A)", 'Mean = ', M, ' ms'
    print "(A, F15.6, A)", 'STD = ', SD, ' ms'
    print "(A, F15.6)", 'CV = ', CV


end program NeuralTractUnitExample