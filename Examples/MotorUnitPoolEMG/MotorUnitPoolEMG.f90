program MotorUnitPoolEMG
    use MotorUnitPoolClass
    use ConfigurationClass
    use ogpf 
    use randomSeedInitialize
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
    character(len = 6) :: pool = 'SOL'
    character(len = 80) :: filename = 'confMotorUnitPoolEMG.rmto'
    type(MotorUnitPool) :: poolSOL
    
    call init_random_seed()
    conf = Configuration(filename)
    

    tf = conf%simDuration_ms
    dt = conf%timeStep_ms
    timeLength = int(tf/dt)
    allocate(t(timeLength))
    allocate(MNv_mV(timeLength))

    poolSOL = MotorUnitPool(conf, pool)
      
    t = [(dt*i, i=1, timeLength)]
    
    call CPU_TIME(tic)
    do i = 1, size(t)
        do j = 1, poolSOL%MUnumber
            poolSOL%iInjected(2*j) = 10.0
        end do
        MNv_mV(i) = poolSol%v_mV(30)        
        call poolSOL%atualizeMotorUnitPool(t(i))
    end do
    call CPU_TIME(toc)

    print '(F15.6, A)', toc - tic, ' seconds'

    if (allocated(poolSOL%unit(50)%somaSpikeTrain)) print '(I3)', size(poolSOL%unit(50)%somaSpikeTrain)

    call poolSol%listSpikes()
    
    call gp%title('Membrane potential of the soma of the MN #15')
    call gp%xlabel('t (ms))')
    call gp%ylabel('Descending command index')
    call gp%plot(t, MNv_mV, 'with line lw 2 lc rgb "#0008B0"')
    
    
    call gp%title('MN spike instants at the soma')
    call gp%xlabel('t (s))')
    call gp%ylabel('Descending command index')
    call gp%plot(poolSOL%poolSomaSpikes(:,1), poolSOL%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    call gp%title('MN spike instants at the terminal')
    call gp%xlabel('t (s))')
    call gp%ylabel('Descending command index')
    call gp%plot(poolSOL%poolTerminalSpikes(:,1), poolSOL%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')
end program MotorUnitPoolEMG