program EPSP
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
    real(wp), dimension(:), allocatable :: t, MNv_mV, RCv_mV
    type(gpf) :: gp
    character(len = 80) :: pool, group
    character(len = 80) :: filename = '../../conf.rmto'
    character(len = 80) :: path = '/home/pablo/osf/Master-Thesis-Data/cell/EPSP/'
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
    real(wp), parameter :: FFConducStrength = 0.023_wp, & 
        declineFactorMN = real(1, wp)/6, declineFactorRC = real(3.5, wp)/3
    character(len=3), parameter :: nS = '1', nFR = '0', &
        nFF = '0', nRC = '1', nCM = '0', nMN = '1' ! nS+nFR+nFF

    call init_random_seed()

    conf = Configuration(filename)
    conf%simDuration_ms = 100

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

        ! Columnar length
        paramTag = 'position:MG-'
        value1 = '0'
        value2 = '0'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'position:RC_ext-'
        value1 = '0'
        value2 = '0'
        call conf%changeConfigurationParameter(paramTag, value1, value2)

    !!!!!!!!!!!!!!!! RC
    ! Turning off spontaneous activity
    paramtag = 'NoiseFunction_RC_ext'
    value1 = '0'
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
    allocate(MNv_mV(timeLength))
    allocate(RCv_mV(timeLength))

    t = [(dt*(i-1), i=1, timeLength)]

    do i = 1, size(t)
        if (t(i)>10.and.t(i)<10.5) then
            motorUnitPools(1)%iInjected(2) = 45.0
        else
            motorUnitPools(1)%iInjected(2) = 0.0
        end if
        do j = 1, size(interneuronPools)
            call interneuronPools(j)%atualizeInterneuronPool(t(i))
            RCv_mV(i) = interneuronPools(j)%v_mV(1)        
        end do
        do j = 1, size(motorUnitPools)
            call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
            MNv_mV(i) = motorUnitPools(j)%v_mV(2)
        end do
    end do

    filename = trim(path) // 'EPSP.dat'
    open(1, file=filename, status = 'replace')
    do i = 1, size(t)
        write(1, '(F10.5, 1X, F10.5)') t(i), RCv_mV(i)
    end do
    close(1)
    !call gp%title('Membrane potential of the soma of the RC #1')
    !call gp%xlabel('t (ms))')
    !call gp%ylabel('Membrane potential (mV)')
    !call gp%plot(t, RCv_mV, 'with line lw 2 lc rgb "#0008B0"') 

    !call gp%title('Membrane potential of the soma of the MN #1')
    !call gp%xlabel('t (ms))')
    !call gp%ylabel('Membrane potential (mV)')
    !call gp%plot(t, MNv_mV, 'with line lw 2 lc rgb "#0008B0"')
    
end program EPSP
