program ImpedanceAnkle
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
    real(wp), dimension(:), allocatable :: t, disturbance
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp) :: FR
    integer :: GammaOrder 
    character(len = 80) :: pool, muscle
    character(len = 80) :: filename = 'confImpedanceAnkle.rmto'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools   
    type(AfferentPool), dimension(:), allocatable:: afferentPools    
    type(jointAnkleForceTask) :: ankle
    real(wp) :: angle
    real(wp) , dimension(200):: Amp, phase
    
    
    call init_random_seed()

    conf = Configuration(filename)

    tf = conf%simDuration_ms
    dt = conf%timeStep_ms
    timeLength = int(tf/dt)

    allocate(t(timeLength))
    allocate(disturbance(timeLength))
    t = [(dt*(i-1), i=1, timeLength)]

    disturbance(:) = 0.0

    call random_number(phase)  
    do i = 1, timeLength
              
        do j = 1, 200            
            disturbance(i) = disturbance(i) + sin(2*pi*50.0/200.0*j*t(i)/1000 + 2*pi*phase(j))
        end do
    end do

    
    

    allocate(afferentPools(4))
    pool = 'Ia'
    muscle = 'SOL'
    afferentPools(1) = AfferentPool(conf, pool, muscle)

    pool = 'Ia'
    muscle = 'MG'
    afferentPools(2) = AfferentPool(conf, pool, muscle)

    pool = 'Ia'
    muscle = 'LG'
    afferentPools(3) = AfferentPool(conf, pool, muscle)

    pool = 'Ia'
    muscle = 'TA'
    afferentPools(4) = AfferentPool(conf, pool, muscle)

    allocate(neuralTractPools(1))
    pool = 'CMExt'
    neuralTractPools(1) = NeuralTract(conf, pool)
    
    allocate(motorUnitPools(4))
    pool = 'SOL'
    motorUnitPools(1) = MotorUnitPool(conf, pool)    

    pool = 'MG'
    motorUnitPools(2) = MotorUnitPool(conf, pool)    

    pool = 'LG'
    motorUnitPools(3) = MotorUnitPool(conf, pool)    

    pool = 'TA'
    motorUnitPools(4) = MotorUnitPool(conf, pool)    


    ankle = jointAnkleForceTask(conf, motorUnitPools)
    allocate(interneuronPools(0))

    synapticNoisePools = synapseFactory(conf, neuralTractPools, &
                                        motorUnitPools, &
                                        interneuronPools, &
                                        afferentPools)
    


    
    
    FR = 1000.0/12.0
    GammaOrder = 10
    
    call cpu_time(tic)
    do i = 2, size(t)
        angle = 0.004*disturbance(i)
        call ankle%atualizeAnkle(t(i), angle)      
        ! do j = 1, size(neuralTractPools)
        !     call neuralTractPools(j)%atualizePool(t(i), FR, GammaOrder)
        ! end do
        do j = 1, 4
            call motorUnitPools(j)%atualizeMotorUnitPool(t(i))
            ! call afferentPools(j)%atualizeAfferentPool(t(i), motorUnitPools(j)%spindle%IaFR_Hz)            
        end do
        call ankle%computeTorque(t(i))
    end do    
    call cpu_time(toc)

    print '(F15.6, A)', toc - tic, ' seconds'

    
    
    
    !call neuralTractPools(1)%listSpikes()
    call motorUnitPools(1)%listSpikes()
    !call afferentPools(1)%listSpikes()
    call motorUnitPools(1)%getMotorUnitPoolEMG()
    
    call gp%title('MN spike instants at the soma')
    call gp%xlabel('t (s))')
    call gp%ylabel('Motoneuron index')
    call gp%plot(motorUnitPools(1)%poolSomaSpikes(:,1), motorUnitPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    call gp%title('MN spike instants at the terminal')
    call gp%xlabel('t (s))')
    call gp%ylabel('Motoneuron index')
    call gp%plot(motorUnitPools(1)%poolTerminalSpikes(:,1), &
    motorUnitPools(1)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    ! call gp%title('Descending command spike instants at the terminal')
    ! call gp%xlabel('t (s))')
    ! call gp%ylabel('Descending command index')
    ! call gp%plot(neuralTractPools(1)%poolTerminalSpikes(:,1), &
    ! neuralTractPools(1)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    ! call gp%title('afferents spike instants at the terminal')
    ! call gp%xlabel('t (s))')
    ! call gp%ylabel('Ia index')
    ! call gp%plot(afferentPools(1)%poolTerminalSpikes(:,1), &
    ! afferentPools(1)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')


    ! call gp%title('Muscle force')
    ! call gp%xlabel('t (ms))')
    ! call gp%ylabel('Force (N)')
    ! call gp%plot(t, motorUnitPools(1)%HillMuscle%force, 'with line lw 2 lc rgb "#0008B0"')

    ! call gp%title('Muscle length')
    ! call gp%xlabel('t (ms))')
    ! call gp%ylabel('length (m)')
    ! call gp%plot(t, motorUnitPools(1)%HillMuscle%length_m, 'with line lw 2 lc rgb "#0008B0"')

    ! call gp%title('Muscle velocity')
    ! call gp%xlabel('t (ms))')
    ! call gp%ylabel('velocity (m/ms)')
    ! call gp%plot(t, motorUnitPools(1)%HillMuscle%velocity_m_ms, 'with line lw 2 lc rgb "#0008B0"')

    call gp%title('Ankle torque')
    call gp%xlabel('t (ms))')
    call gp%ylabel('torque (N.m)')
    call gp%plot(t, ankle%ankleTorque_Nm, 'with line lw 2 lc rgb "#0008B0"')

    call gp%title('Ankle angle')
    call gp%xlabel('t (ms))')
    call gp%ylabel('angle (degree)')
    call gp%plot(t, ankle%ankleAngle_rad*180/pi, 'with line lw 2 lc rgb "#0008B0"')
    
    
end program ImpedanceAnkle