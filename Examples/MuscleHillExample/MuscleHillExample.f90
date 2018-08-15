program MuscleHillExample
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
    real(wp), dimension(:), allocatable :: t
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp) :: FR
    integer :: GammaOrder 
    character(len = 80) :: pool
    character(len = 80) :: filename = 'confMuscleHillExample.rmto'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools   
    type(AfferentPool), dimension(:), allocatable:: afferentPools    
    type(jointAnkleForceTask) :: ankle
    real(wp) :: angle
    
    
    call init_random_seed()

    conf = Configuration(filename)
    allocate(afferentPools(0))
    allocate(neuralTractPools(1))
    pool = 'CMExt'
    neuralTractPools(1) = NeuralTract(conf, pool)
    allocate(motorUnitPools(1))
    pool = 'SOL'
    motorUnitPools(1) = MotorUnitPool(conf, pool)    
    ankle = jointAnkleForceTask(conf, motorUnitPools)
    allocate(interneuronPools(0))

    synapticNoisePools = synapseFactory(conf, neuralTractPools, &
                                        motorUnitPools, &
                                        interneuronPools, &
                                        afferentPools)
    
    tf = conf%simDuration_ms
    dt = conf%timeStep_ms
    timeLength = int(tf/dt)

    allocate(t(timeLength))
    
    t = [(dt*(i-1), i=1, timeLength)]
    
    FR = 1000.0/12.0
    GammaOrder = 10
    
    call cpu_time(tic)
    do i = 2, size(t)
        angle = 0.1*sin(2*pi*t(i)/1000.0)    
        call ankle%atualizeAnkle(t(i), angle)      
        do j = 1, size(neuralTractPools)
            call neuralTractPools(j)%atualizePool(t(i), FR, GammaOrder)
        end do
        do j = 1, size(motorUnitPools)
            call motorUnitPools(j)%atualizeMotorUnitPool(t(i))
        end do
    end do    
    call cpu_time(toc)

    print '(F15.6, A)', toc - tic, ' seconds'

    
    
    
    call neuralTractPools(1)%listSpikes()
    call motorUnitPools(1)%listSpikes()
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

    call gp%title('MN spike instants at the terminal')
    call gp%xlabel('t (s))')
    call gp%ylabel('Motoneuron index')
    call gp%plot(neuralTractPools(1)%poolTerminalSpikes(:,1), &
    neuralTractPools(1)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')


    call gp%title('Muscle force')
    call gp%xlabel('t (ms))')
    call gp%ylabel('Force (N)')
    call gp%plot(t, motorUnitPools(1)%HillMuscle%force, 'with line lw 2 lc rgb "#0008B0"')

    call gp%title('Muscle length')
    call gp%xlabel('t (ms))')
    call gp%ylabel('length (m)')
    call gp%plot(t, motorUnitPools(1)%HillMuscle%length_m, 'with line lw 2 lc rgb "#0008B0"')

    call gp%title('Muscle velocity')
    call gp%xlabel('t (ms))')
    call gp%ylabel('velocity (m/ms)')
    call gp%plot(t, motorUnitPools(1)%HillMuscle%velocity_m_ms, 'with line lw 2 lc rgb "#0008B0"')

    call gp%title('Ankle angle')
    call gp%xlabel('t (ms))')
    call gp%ylabel('angle (degree)')
    call gp%plot(t, ankle%ankleAngle_rad, 'with line lw 2 lc rgb "#0008B0"')
    
    
end program MuscleHillExample