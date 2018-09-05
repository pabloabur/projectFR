program PosturalControl
    use MotorUnitPoolClass
    use NeuralTractClass
    use InterneuronPoolClass
    use SynapticNoiseClass
    use ConfigurationClass
    use ogpf 
    use randomSeedInitialize
    use SynapsesFactoryModule
    use postureTaskClass
    use AfferentPoolClass
    implicit none 
    include 'mkl.fi'
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
    character(len = 80) :: pool, muscle
    character(len = 80) :: filename = 'confPosturalControl.rmto'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools   
    type(AfferentPool), dimension(:), allocatable:: afferentPools    
    type(postureTask) :: body
    real(wp) :: angle
    real(wp) , dimension(200):: Amp, phase
    logical :: continueFlag
    real(wp) :: dynGamma, statGamma
    integer :: subject, trial


    call init_random_seed()

    conf = Configuration(filename)

    tf = conf%simDuration_ms
    dt = conf%timeStep_ms
    timeLength = int(tf/dt)

    allocate(t(timeLength))
    
    t = [(dt*(i-1), i=1, timeLength)]

    
    FR = 73.0
    GammaOrder = 25
    
    filename = 'GammaFallTime.txt'
    open(1, file=filename, status = 'replace') 

    subject = 1
    do while (subject <= 1)
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


        body = postureTask(conf, motorUnitPools)
        allocate(interneuronPools(0))

        synapticNoisePools = synapseFactory(conf, neuralTractPools, &
                                            motorUnitPools, &
                                            interneuronPools, &
                                            afferentPools)
        call random_number(dynGamma)
        dynGamma = 30+4*dynGamma
        call random_number(statGamma)
        statGamma = 30 + 4*statGamma
        print '(A, F15.6)', 'Dynamic Gamma = ',  dynGamma
        print '(A, F15.6)', 'Static Gamma = ', statGamma
        continueFlag = .true.
        trial = 1
        do while (continueFlag .and. trial <= 2)
            continueFlag = .false.
            tic =  dsecnd()
            i = 1
            
            do while (i <= size(t))
                do j = 1, size(neuralTractPools)
                    call neuralTractPools(j)%atualizePool(t(i), FR, GammaOrder)
                end do
                do j = 1, 4
                    call motorUnitPools(j)%atualizeMotorUnitPool(t(i), dynGamma, statGamma)
                    call afferentPools(j)%atualizeAfferentPool(t(i), motorUnitPools(j)%spindle%IaFR_Hz)            
                end do
                call body%atualizeBody(t(i))
                if (abs(body%ankleAngle_rad(i)) > pi/6.0) then
                    print '(A, 1X, F15.6, 1X, A)', 'Body fell after ', t(i)/1000.0, ' seconds'
                    continueFlag = .true.                
                    write(1, '(F15.6, 1X, F15.6,1X,F15.6, 1X, I2, 1X, I2)') ([t(i), dynGamma, statGamma]), ([subject, trial])
                    i = size(t)
                    call random_number(dynGamma)
                    dynGamma = 30 + 4*dynGamma
                    call random_number(statGamma)
                    statGamma = 30 + 4*statGamma
                    print '(A, F15.6)', 'Dynamic Gamma = ',  dynGamma
                    print '(A, F15.6)', 'Static Gamma = ', statGamma
                end if 
                i = i + 1
            end do    
            toc =  dsecnd()

            print '(F15.6, A)', toc - tic, ' seconds'       
            
            ! call neuralTractPools(1)%listSpikes()
            ! call motorUnitPools(1)%listSpikes()
            ! call afferentPools(1)%listSpikes()
            ! !call motorUnitPools(1)%getMotorUnitPoolEMG()
            
            ! call gp%title('MN spike instants at the soma')
            ! call gp%xlabel('t (s))')
            ! call gp%ylabel('Motoneuron index')
            ! call gp%plot(motorUnitPools(1)%poolSomaSpikes(:,1), &
            ! motorUnitPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

            ! call gp%title('AF spike instants at the soma')
            ! call gp%xlabel('t (s))')
            ! call gp%ylabel('Afferent index')
            ! call gp%plot(afferentPools(1)%poolTerminalSpikes(:,1), &
            ! afferentPools(1)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

            ! call gp%title('MN spike instants at the terminal')
            ! call gp%xlabel('t (s))')
            ! call gp%ylabel('Motoneuron index')
            ! call gp%plot(motorUnitPools(1)%poolTerminalSpikes(:,1), &
            ! motorUnitPools(1)%poolTerminalSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

            ! call gp%title('Muscle force')
            ! call gp%xlabel('t (ms))')
            ! call gp%ylabel('Force (N)')
            ! call gp%plot(t, motorUnitPools(1)%HillMuscle%force, 'with line lw 2 lc rgb "#0008B0"')

            ! call gp%title('Ankle torque')
            ! call gp%xlabel('t (ms))')
            ! call gp%ylabel('torque (N.m)')
            ! call gp%plot(t, body%ankleTorque_Nm, 'with line lw 2 lc rgb "#0008B0"')

            ! call gp%title('Ankle angle')
            ! call gp%xlabel('t (ms))')
            ! call gp%ylabel('angle (degree)')
            ! call gp%plot(t, body%ankleAngle_rad*180.0/pi, 'with line lw 2 lc rgb "#0008B0"')

            do j = 1, size(neuralTractPools)
                call neuralTractPools(j)%reset()
            end do            
            do j = 1, 4
                call motorUnitPools(j)%reset()
            end do            
            do j = 1, 4
                call afferentPools(j)%reset()
            end do            
            call body%reset()
            i = 1
            trial = trial + 1
        end do
        
        write(1, '(F15.6, 1X, F15.6,1X,F15.6, 1X, I2, 1X, I2)') ([t(i), dynGamma, statGamma]), ([subject, trial])
        subject = subject + 1
        deallocate(afferentPools)
        deallocate(motorUnitPools)
        deallocate(neuralTractPools)
        deallocate(synapticNoisePools)
        deallocate(interneuronPools)
    end do
    close(1)
    
end program PosturalControl