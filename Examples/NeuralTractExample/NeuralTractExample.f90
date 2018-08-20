program NeuralTractExample
    use NeuralTractClass
    use ConfigurationClass
    use ogpf 
    use randomSeedInitialize
    implicit none 
    !integer, parameter :: wp = kind(1.0d0)
    type(Configuration) :: conf
    type(NeuralTract) :: neuralTractPool
    real(wp), parameter :: ISI_ms = 12.0
    real(wp), parameter :: FR_Hz = 1000.0/ISI_ms
    integer, parameter :: GammaOrder = 1
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    !use, intrinsic :: iso_fortran_env, only : wp => real64
    real(wp), parameter :: dt = 0.05
    real(wp), parameter :: tf = 10000
    integer, parameter :: timeLength = int(tf/dt)
    integer :: i
    real(wp), dimension(timeLength) :: t
    real(wp) :: tic, toc
    type(gpf) :: gp
    character(len = 6) :: pool = 'CMExt'
    character(len = 80) :: filename = 'confNeuralTractExample.rmto'
    character(len=30)::paramTag
    real(wp) :: parammReal
    integer :: parammInt
    character(len = 80) :: paramChar
    real(wp), dimension(:), allocatable :: spikesIndices
    real(wp), dimension(:), allocatable :: spikesInstants
    

    conf = Configuration(filename)
    neuralTractPool = NeuralTract(conf, pool)
    
    t = [(dt*i, i=1, timeLength)]
   
    call init_random_seed()


    call cpu_time(tic)
    do i = 1, timeLength
        call neuralTractPool%atualizePool(t(i), FR_Hz, GammaOrder)
    enddo
    call cpu_time(toc)
    print "(F15.6, A)", toc - tic,  " seconds"

    call neuralTractPool%listSpikes()

    
    allocate(spikesIndices(size(neuralTractPool%poolTerminalSpikes(:,2))))
    allocate(spikesInstants(size(neuralTractPool%poolTerminalSpikes(:,1))))

    spikesIndices = neuralTractPool%poolTerminalSpikes(:,2)
    spikesInstants = neuralTractPool%poolTerminalSpikes(:,1)

    call gp%title('Neural tract spikes instants')
    call gp%xlabel('t (s))')
    call gp%ylabel('Descending command index')
    call gp%plot(spikesInstants, spikesIndices, 'with points pt 5 lc rgb "#0008B0"')

end program NeuralTractExample