program MotorUnitPoolEMG
    use MotorUnitPoolClass
    use NeuralTractClass
    use InterneuronPoolClass
    use ConfigurationClass
    use ogpf 
    use randomSeedInitialize
    use SynapsesFactoryModule
    use SynapticNoiseClass
    use MuscleNoHillClass
    implicit none 
    !integer, parameter :: wp = kind(1.0d0)
    type(Configuration) :: conf
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    real(wp) :: dt
    real(wp) :: tf
    integer :: timeLength
    integer :: i, j
    real(wp), dimension(:), allocatable :: t, MNv_mV
    real(wp) :: tic, toc
    type(gpf) :: gp
    character(len = 6) :: pool = 'SOL'
    character(len = 80) :: filename = 'confMotorUnitPoolEMG.rmto'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools    

    
    

    call init_random_seed()
    conf = Configuration(filename)
    

    tf = conf%simDuration_ms
    dt = conf%timeStep_ms
    timeLength = int(tf/dt)
    allocate(t(timeLength))
    allocate(MNv_mV(timeLength))

    allocate(neuralTractPools(0))
    
    allocate(interneuronPools(0))
    allocate(motorUnitPools(1))

    

    motorUnitPools(1) = MotorUnitPool(conf, pool)
      
    t = [(dt*i, i=1, timeLength)]
    
    call cpu_time(tic)
    do i = 1, size(t)
        do j = 1, motorUnitPools(1)%MUnumber
            motorUnitPools(1)%iInjected(2*j) = 40.0
        end do
        MNv_mV(i) = motorUnitPools(1)%v_mV(30)        
        call motorUnitPools(1)%atualizeMotorUnitPool(t(i))
    end do
    call cpu_time(toc)

    print '(F15.6, A)', toc - tic, ' seconds'

    

    call motorUnitPools(1)%listSpikes()
    call motorUnitPools(1)%getMotorUnitPoolEMG()
    
    call gp%title('Membrane potential of the soma of the MN #15')
    call gp%xlabel('t (ms))')
    call gp%ylabel('Descending command index')
    call gp%plot(t, MNv_mV, 'with line lw 2 lc rgb "#0008B0"')
    
    
    call gp%title('MN spike instants at the soma')
    call gp%xlabel('t (s))')
    call gp%ylabel('Descending command index')
    call gp%plot(motorUnitPools(1)%poolSomaSpikes(:,1), &
    motorUnitPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    call gp%title('MN spike instants at the terminal')
    call gp%xlabel('t (s))')
    call gp%ylabel('Descending command index')
    call gp%plot(motorUnitPools(1)%poolTerminalSpikes(:,1), &
    motorUnitPools(1)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    call gp%title('Muscle force')
    call gp%xlabel('t (ms))')
    call gp%ylabel('Force (N)')
    call gp%plot(t, motorUnitPools(1)%NoHillMuscle%force, 'with line lw 2 lc rgb "#0008B0"')
    
    call gp%title('Muscle EMG')
    call gp%xlabel('t (ms))')
    call gp%ylabel('EMG (mV)')
    call gp%plot(t, motorUnitPools(1)%emg, 'with line lw 2 lc rgb "#0008B0"')
end program MotorUnitPoolEMG