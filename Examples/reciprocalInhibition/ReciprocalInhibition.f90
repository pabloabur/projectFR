program ReciprocalInhibition
    use MotorUnitPoolClass
    use NeuralTractClass
    use InterneuronPoolClass
    use SynapticNoiseClass
    use ConfigurationClass
    use ogpf 
    use randomSeedInitialize
    use SynapsesFactoryModule
    use jointAnkleForceTaskClass
    use AfferentPoolClass
    implicit none 
    !integer, parameter :: wp = kind(1.0d0)
    type(Configuration) :: conf
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    real(wp) :: dt
    real(wp) :: tf
    integer :: timeLength
    integer :: i, j
    real(wp), dimension(:), allocatable :: t, MN_VmV
    real(wp) :: tic, toc
    type(gpf) :: gp
    character(len = 80) :: pool, muscle
    character(len = 80) :: filename = 'confReciprocalInhibition.rmto'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools   
    type(AfferentPool), dimension(:), allocatable:: afferentPools    
    type(jointAnkleForceTask) :: ankle
    real(wp) :: IaFR
    
    
    call init_random_seed()

    conf = Configuration(filename)
    allocate(afferentPools(3))
    pool = 'Ia'
    muscle = 'SOL'
    afferentPools(1) = AfferentPool(conf, pool, muscle)
    
    pool = 'Ia'
    muscle = 'LG'
    afferentPools(2) = AfferentPool(conf, pool, muscle)

    pool = 'Ia'
    muscle = 'TA'
    afferentPools(3) = AfferentPool(conf, pool, muscle)

    allocate(neuralTractPools(0))
    
    allocate(motorUnitPools(1))
    pool = 'SOL'
    motorUnitPools(1) = MotorUnitPool(conf, pool)    

    allocate(interneuronPools(1))
    pool = 'IaIn'
    muscle = 'ext'
    interneuronPools(1) = InterneuronPool(conf, pool, muscle)
   
    

    synapticNoisePools = synapseFactory(conf, neuralTractPools, &
                                        motorUnitPools, &
                                        interneuronPools, &
                                        afferentPools)
    
    tf = conf%simDuration_ms
    dt = conf%timeStep_ms
    timeLength = int(tf/dt)

    allocate(t(timeLength))
    allocate(MN_VmV(timeLength))
    
    t = [(dt*(i-1), i=1, timeLength)]
    
    IaFR = 0.0
    call cpu_time(tic)
    do i = 1, size(t)
        do j = 1, size(motorUnitPools)
            call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 30.0_wp, 30.0_wp)
        end do
        do j = 1, size(afferentPools)            
            call afferentPools(j)%atualizeAfferentPool(t(i), IaFR)
        end do
        call interneuronPools(1)%atualizeInterneuronPool(t(i))
        MN_VmV(i) = motorUnitPools(1)%v_mV(50)
    end do    
    call cpu_time(toc)

    print '(F15.6, A)', toc - tic, ' seconds'
    
    call interneuronPools(1)%listSpikes()
    call motorUnitPools(1)%listSpikes()
    call afferentPools(1)%listSpikes()
    call motorUnitPools(1)%getMotorUnitPoolEMG()
    
    call gp%title('EMG')
    call gp%xlabel('t (ms))')
    call gp%ylabel('EMG (mV)')
    call gp%plot(t, motorUnitPools(1)%emg, 'with line lw 2 lc rgb "#0008B0"')
    
    call gp%title('Membrane potential')
    call gp%xlabel('t (ms))')
    call gp%ylabel('V (mV)')
    call gp%plot(t, MN_VmV, 'with line lw 2 lc rgb "green"')
   
    
    call gp%title('IN spike instants at the soma')
    call gp%xlabel('t (ms))')
    call gp%ylabel('IN index')
    call gp%plot(interneuronPools(1)%poolSomaSpikes(:,1), &
    interneuronPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    call gp%title('MN spike instants at the terminal')
    call gp%xlabel('t (ms))')
    call gp%ylabel('MN index')
    call gp%plot(motorunitPools(1)%poolTerminalSpikes(:,1), &
    motorunitPools(1)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    call gp%title('IN spike instants at the soma')
    call gp%xlabel('t (ms))')
    call gp%ylabel('IN index')
    call gp%plot(x1=t, &
    y1=motorUnitPools(1)%unit(1)%nerveStimulus_mA, &
    ls1='with line lw 2 lc rgb "blue"', &
    x2 = t, y2 = afferentPools(3)%unit(1)%nerveStimulus_mA, &
    ls2='with line lw 2 lc rgb "red"')
    
    
end program ReciprocalInhibition