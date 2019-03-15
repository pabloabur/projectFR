program IPSP
    use MotorUnitPoolClass
    use NeuralTractClass
    use InterneuronPoolClass
    use SynapticNoiseClass
    use ConfigurationClass
    use ogpf 
    use randomSeedInitialize
    use CharacterArrayClass
    use CharacterMatrixClass
    use QueueClass
    use AfferentPoolClass
    use SynapsesFactoryModule

    implicit none 
    !integer, parameter :: wp = kind(1.0d0)
    type(Configuration) :: conf
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    integer :: timeLength
    integer :: i, j
    real(wp), dimension(:), allocatable :: t, MNSv_mV, MNFRv_mV, MNFFv_mV
    real(wp), dimension(:), allocatable :: RC1, RC2, RC3
    type(gpf) :: gp
    character(len = 80) :: pool, group
    character(len = 80) :: filename = '../../conf.rmto'
    character(len = 45) :: path = '/home/pablo/osf/Master-Thesis-Data/cell/'
    character(len = 25) :: folderName = 'IPSP/'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools     
    type(AfferentPool), dimension(:), allocatable:: afferentPools     
    character(len=5) :: param
    character(len=80) :: paramTag
    character(len=80) :: value1, value2
    real(wp) :: dt
    real(wp) :: tf
    ! Input parameters
    logical, parameter :: probDecay = .false.
    character(len=3), parameter :: nS = '2', nFR = '2', &
        nFF = '2', nRC = '6', nCM = '0', nMN = '6' ! nS+nFR+nFF

    call init_random_seed()

    conf = Configuration(filename)
    conf%simDuration_ms = 60

    !Changing configuration file
    paramTag = 'Number_CMExt'
    value1 = nCM
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'MUnumber_MG-S'
    value1 = nS
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'MUnumber_MG-FR'
    value1 = nFR
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'MUnumber_MG-FF'
    value1 = nFF
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'Number_RC_ext'
    value1 = nRC
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'Number_MG'
    value1 = nMN
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    !!!!!!!!!!!!!!!! Done to ensure that neighboring neurons will not affect each other
    paramtag = 'position:MG-'
    value1 = '0'
    value2 = '1200'
    call conf%changeconfigurationparameter(paramtag, value1, value2)
    paramtag = 'position:RC_ext-'
    value1 = '0'
    value2 = '1200'
    call conf%changeconfigurationparameter(paramtag, value1, value2)

    !!!!!!!!!!!!!!!! RC
    ! Turning off spontaneous activity
    paramtag = 'NoiseFunction_RC_ext'
    value1 = '0'!'7'
    value2 = ''
    call conf%changeconfigurationparameter(paramtag, value1, value2)
    ! Dynamics of MN-RC synapse
    paramtag = 'dyn:MG-S>RC_ext-@soma|excitatory'
    value1 = 'None'
    value2 = ''
    call conf%changeconfigurationparameter(paramtag, value1, value2)
    paramTag = 'dyn:MG-FR>RC_ext-@soma|excitatory'
    value1 = 'None'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'dyn:MG-FF>RC_ext-@soma|excitatory'
    value1 = 'None'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    !!!!!!!!!!!!!!!! Independent noise
    ! TODO try high FR, low g, and vice versa
    paramTag = 'NoiseFunction_MG'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    ! Removing influence of stimulus (required)
    paramTag = 'stimIntensity_PTN'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    print *, 'Building neural elements'
    allocate(neuralTractPools(0))

    allocate(motorUnitPools(1))
    pool = 'MG'
    motorUnitPools(1) = MotorUnitPool(conf, pool)    

    allocate(interneuronPools(1))
    pool = 'RC'
    group = 'ext'    
    interneuronPools(1) = InterneuronPool(conf, pool, group)

    allocate(afferentPools(0))

    synapticNoisePools = synapseFactory(conf, neuralTractPools, &
                                        motorUnitPools, &
                                        interneuronPools, &
                                        afferentPools, probDecay)
    
    tf = conf%simDuration_ms
    dt = conf%timeStep_ms
    timeLength = nint(tf/dt)
    
    allocate(t(timeLength))
    allocate(MNFRv_mV(timeLength))
    allocate(MNFFv_mV(timeLength))
    allocate(MNSv_mV(timeLength))
    allocate(RC1(timeLength))
    allocate(RC2(timeLength))
    allocate(RC3(timeLength))

    t = [(dt*(i-1), i=1, timeLength)]

    do i = 1, size(t)
        if (t(i)>10.0.and.t(i)<10.5) then
            interneuronPools(1)%iInjected(1) = 6.0
            interneuronPools(1)%iInjected(3) = 6.0
            interneuronPools(1)%iInjected(5) = 6.0
        else
            interneuronPools(1)%iInjected(:) = 0.0
        end if
        do j = 1, size(interneuronPools)
            call interneuronPools(j)%atualizeInterneuronPool(t(i))
            RC1(i) = interneuronPools(j)%v_mV(1)
            RC2(i) = interneuronPools(j)%v_mV(3)
            RC3(i) = interneuronPools(j)%v_mV(5)
        end do
        do j = 1, size(motorUnitPools)
            call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
            MNSv_mV(i) = motorUnitPools(j)%v_mV(2*(1))
            MNFRv_mV(i) = motorUnitPools(j)%v_mV(2*(3))
            MNFFv_mV(i) = motorUnitPools(j)%v_mV(2*(5))
        end do
    end do

    print *, 'Position of RC 1'
    print *, interneuronPools(1)%unit(1)%position_mm
    print *, 'Position of RC 3'
    print *, interneuronPools(1)%unit(3)%position_mm
    print *, 'Position of RC 5'
    print *, interneuronPools(1)%unit(5)%position_mm
    print *, 'Position of MN 1'
    print *, motorUnitPools(1)%unit(1)%position_mm
    print *, 'Position of MN 3'
    print *, motorUnitPools(1)%unit(3)%position_mm
    print *, 'Position of MN 5'
    print *, motorUnitPools(1)%unit(5)%position_mm
    call interneuronPools(1)%listSpikes()
    call gp%title('RC spike instants at the soma')
    call gp%xlabel('t (s))')
    call gp%ylabel('Motoneuron index')
    call gp%plot(interneuronPools(1)%poolSomaSpikes(:,1), &
    interneuronPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')
    call gp%title('Membrane potential of the soma of the RC #1')
    call gp%xlabel('t (ms))')
    call gp%ylabel('Descending command index')
    call gp%plot(t, RC1, 'with line lw 2 lc rgb "#0008B0"')  
    call gp%title('Membrane potential of the soma of the RC #1')
    call gp%xlabel('t (ms))')
    call gp%ylabel('Descending command index')
    call gp%plot(t, RC2, 'with line lw 2 lc rgb "#0008B0"')  
    call gp%title('Membrane potential of the soma of the RC #1')
    call gp%xlabel('t (ms))')
    call gp%ylabel('Descending command index')
    call gp%plot(t, RC3, 'with line lw 2 lc rgb "#0008B0"')  

    !call gp%title('Membrane potential of the soma of the MN #1')
    !call gp%xlabel('t (ms))')
    !call gp%ylabel('Membrane potential (mV)')
    !call gp%plot(t, MNSv_mV, 'with line lw 2 lc rgb "#0008B0"')
    !
    !call gp%title('Membrane potential of the soma of the MN #1')
    !call gp%xlabel('t (ms))')
    !call gp%ylabel('Membrane potential (mV)')
    !call gp%plot(t, MNFRv_mV, 'with line lw 2 lc rgb "#0008B0"')

    !call gp%title('Membrane potential of the soma of the MN #1')
    !call gp%xlabel('t (ms))')
    !call gp%ylabel('Membrane potential (mV)')
    !call gp%plot(t, MNFFv_mV, 'with line lw 2 lc rgb "#0008B0"')

    filename = trim(path) // trim(folderName) // trim("IPSP.dat")
    open(1, file=filename, status = 'replace')
    do i = 1, size(t)
        write(1, '(F8.3, 1X, F15.10, 1X, F15.10, 1X, F15.10)') t(i), MNSv_mV(i), &
            MNFRv_mV(i), MNFFv_mV(i)
    end do
    close(1)
end program IPSP
