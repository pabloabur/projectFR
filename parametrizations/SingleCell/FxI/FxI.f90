program FxI
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
    integer :: i, j, k, l
    real(wp), dimension(:), allocatable :: t, RCv_mV
    type(gpf) :: gp
    character(len = 80) :: pool, group
    character(len = 80) :: filename = '../../conf.rmto'
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
    real(wp) :: firingRate_pps(4)
    real(wp), allocatable :: steadyState(:)
    real(wp) :: auxSteadyState
    real(wp) :: st(30)
    real(wp) :: nd(30)
    real(wp) :: rd(30)
    real(wp) :: th(30)
    integer :: arrLen
    ! Input parameters
    logical, parameter :: probDecay = .false.
    real(wp), parameter :: FFConducStrength = 0.3_wp, & 
        declineFactorMN = real(1, wp)/6, declineFactorRC = real(3.5, wp)/3
    character(len=3), parameter :: nS = '0', nFR = '0', &
        nFF = '0', nRC = '1', nCM = '0', nMN = '0' ! nS+nFR+nFF
    real(wp) :: currents(30)
    do i = 1, size(currents)
        currents(i) = real(i, wp)/2
    end do

    call init_random_seed()

    conf = Configuration(filename)
    conf%simDuration_ms = 150

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

    call get_command_argument(1, param)
    if (param.eq.'old') then
        ! Connectivity
        paramTag = 'Con:RC_ext->MG-S@dendrite|inhibitory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:RC_ext->MG-FR@dendrite|inhibitory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:RC_ext->MG-FF@dendrite|inhibitory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:MG-S>RC_ext-@soma|excitatory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:MG-FR>RC_ext-@soma|excitatory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:MG-FF>RC_ext-@soma|excitatory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Conductances
        paramTag = 'gmax:RC_ext->MG-S@dendrite|inhibitory'
        value1 = '0.44'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:RC_ext->MG-FR@dendrite|inhibitory'
        value1 = '0.3'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:RC_ext->MG-FF@dendrite|inhibitory'
        value1 = '0.24'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:MG-S>RC_ext-@soma|excitatory'
        value1 = '0.15'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:MG-FR>RC_ext-@soma|excitatory'
        value1 = '0.17'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:MG-FF>RC_ext-@soma|excitatory'
        value1 = '0.3'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Morphology
        paramTag = 'd@soma:RC_ext-'
        value1 = '64.77885'
        value2 = '64.77885'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'l@soma:RC_ext-'
        value1 = '285'
        value2 = '285'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'res@soma:RC_ext-'
        value1 = '200'
        value2 = '200'
        call conf%changeConfigurationParameter(paramTag, value1, value2)

    else if (param.eq.'new') then
        ! Threshold
        paramTag = 'threshold:RC_ext-'
        value1 = '5'
        value2 = '15'
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Connectivity
        paramTag = 'Con:RC_ext->MG-S@dendrite|inhibitory'
        value1 = '4'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:RC_ext->MG-FR@dendrite|inhibitory'
        value1 = '4'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:RC_ext->MG-FF@dendrite|inhibitory'
        value1 = '4'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:MG-S>RC_ext-@soma|excitatory'
        value1 = '6'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:MG-FR>RC_ext-@soma|excitatory'
        value1 = '6'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:MG-FF>RC_ext-@soma|excitatory'
        value1 = '6'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Conductances
        paramTag = 'gmax:RC_ext->MG-S@dendrite|inhibitory'
        value1 = '0.44'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:RC_ext->MG-FR@dendrite|inhibitory'
        value1 = '0.44'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:RC_ext->MG-FF@dendrite|inhibitory'
        value1 = '0.44'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:MG-S>RC_ext-@soma|excitatory'
        value1 = '0.15'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:MG-FR>RC_ext-@soma|excitatory'
        value1 = '0.15'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:MG-FF>RC_ext-@soma|excitatory'
        value1 = '0.15'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Morphology
        paramTag = 'd@soma:RC_ext-'
        value1 = '25'
        value2 = '25'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'l@soma:RC_ext-'
        value1 = '242'
        value2 = '242'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'res@soma:RC_ext-'
        value1 = '760'
        value2 = '760'
        call conf%changeConfigurationParameter(paramTag, value1, value2)

    else if (param.eq.'final') then
        ! Threshold
        paramTag = 'threshold:RC_ext-'
        value1 = '18.9089'
        value2 = '18.9089'
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Connectivity
        paramTag = 'Con:RC_ext->MG-S@dendrite|inhibitory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:RC_ext->MG-FR@dendrite|inhibitory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:RC_ext->MG-FF@dendrite|inhibitory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:MG-S>RC_ext-@soma|excitatory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:MG-FR>RC_ext-@soma|excitatory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:MG-FF>RC_ext-@soma|excitatory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Conductances
        paramTag = 'gmax:RC_ext->MG-S@dendrite|inhibitory'
        value1 = '0.130'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:RC_ext->MG-FR@dendrite|inhibitory'
        value1 = '0.119'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:RC_ext->MG-FF@dendrite|inhibitory'
        value1 = '0.081'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:MG-S>RC_ext-@soma|excitatory'
        write(value1, '(F15.6)') FFConducStrength/2.2
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:MG-FR>RC_ext-@soma|excitatory'
        write(value1, '(F15.6)') FFConducStrength/1.8
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:MG-FF>RC_ext-@soma|excitatory'
        write(value1, '(F15.6)') FFConducStrength
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Morphology
        paramTag = 'd@soma:RC_ext-'
        value1 = '27'
        value2 = '27'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'l@soma:RC_ext-'
        value1 = '218.2168'
        value2 = '218.2168'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'res@soma:RC_ext-'
        value1 = '7000'
        value2 = '7000'
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Ks
        paramTag = 'gmax_Kf:RC_ext-@soma'
        value1 = '3300'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax_Ks:RC_ext-@soma'
        value1 = '2300000'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'beta_q:RC_ext-@soma'
        value1 = '0.02'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'alpha_q:RC_ext-@soma'
        value1 = '0.004'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'alpha_n:RC_ext-@soma'
        value1 = '6'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'beta_n:RC_ext-@soma'
        value1 = '0.5'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Decay factors
        paramTag = 'dec:MG-S>RC_ext-@soma|excitatory'
        write(value1, '(F15.6)') declineFactorMN
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'dec:MG-FR>RC_ext-@soma|excitatory'
        write(value1, '(F15.6)') declineFactorMN
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'dec:MG-FF>RC_ext-@soma|excitatory'
        write(value1, '(F15.6)') declineFactorMN
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'dec:RC_ext->MG-S@dendrite|inhibitory'
        write(value1, '(F15.6)') declineFactorRC
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'dec:RC_ext->MG-FR@dendrite|inhibitory'
        write(value1, '(F15.6)') declineFactorRC
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'dec:RC_ext->MG-FF@dendrite|inhibitory'
        write(value1, '(F15.6)') declineFactorRC
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Columnar length
        paramTag = 'position:MG-'
        value1 = '0'
        value2 = '7'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'position:RC_ext-'
        value1 = '0'
        value2 = '7'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
    else
        print *, 'Wrong parametrization option'
        stop (1)
    endif

    !!!!!!!!!!!!!!!! RC
    ! Turning off spontaneous activity
    paramtag = 'gmax:Noise>RC_ext-@soma|excitatory'
    value1 = '0.08'
    value2 = ''
    call conf%changeconfigurationparameter(paramtag, value1, value2)
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

    allocate(motorUnitPools(0))

    allocate(interneuronPools(1))
    pool = 'RC'
    group = 'ext'    
    interneuronPools(1) = InterneuronPool(conf, pool, group)

    allocate(afferentPools(0))

    synapticNoisePools = synapseFactory(conf, neuralTractPools, &
                                        motorUnitPools, &
                                        interneuronPools, &
                                        afferentPools)!, probDecay)
    
    tf = conf%simDuration_ms
    dt = conf%timeStep_ms
    timeLength = nint(tf/dt)
    
    allocate(t(timeLength))
    allocate(RCv_mV(timeLength))

    t = [(dt*(i-1), i=1, timeLength)]

    do k = 1, size(currents)
        ! Running simulation
        do i = 1, size(t)
            if (t(i)>10.0.and.t(i)<100.0) then
                interneuronPools(1)%iInjected(:) = currents(k)
            else
                interneuronPools(1)%iInjected(:) = 0.0
            end if
            do j = 1, size(interneuronPools)
                call interneuronPools(j)%atualizeInterneuronPool(t(i))
                RCv_mV(i) = interneuronPools(j)%v_mV(1)        
            end do
        end do

        call interneuronPools(1)%listSpikes()

        if (size(interneuronPools(1)%poolSomaSpikes(:,1))>0) then
            ! Getting first three intervals
            do l = 1, 3
                firingRate_pps(l) = 1000/(interneuronPools(1)%poolSomaSpikes(l+1,1)-&
                                    interneuronPools(1)%poolSomaSpikes(l,1))
            end do

            ! The arrLen count is done to avoid using a dynamic array
            arrLen = 0
            do l = 1, size(interneuronPools(1)%poolSomaSpikes(:,1))
                if (interneuronPools(1)%poolSomaSpikes(l,1)>80.0) then
                    arrLen = arrLen + 1
                end if
            end do

            ! When there are 1 spike in steady state, firing rate is calculated
            ! considering the last spike instant before that steady state spike
            if (arrLen == 0) then
                print *, 'I though this would not happen'
                stop (1)
            else if (arrLen == 1) then
                arrLen = arrLen + 1
            end if

            allocate(steadyState(arrLen))
            ! If I subtract the size of poolSomaSpikes by arrLen and add 1, I get the
            ! the first index in which spike time>90.0
            steadyState = interneuronPools(1)%poolSomaSpikes(size(interneuronPools(1)%&
                            poolSomaSpikes(:,1)) - arrLen + 1:,1)
            auxSteadyState = 0.0
            do l = 1, arrLen - 1
                auxSteadyState = auxSteadyState + real(1000, wp)/(steadyState(l+1)-steadyState(l))
            end do
            firingRate_pps(4) = auxSteadyState/(arrLen - 1)
        else 
            do l = 1, size(firingRate_pps)
                firingRate_pps(l) = 0
            end do
        end if
        call interneuronPools(1)%reset()
        if (allocated(steadyState)) deallocate(steadyState)

        ! Getting values
        st(k) = firingRate_pps(1)
        nd(k) = firingRate_pps(2)
        rd(k) = firingRate_pps(3)
        th(k) = firingRate_pps(4)
    end do


    call gp%title('1')
    call gp%xlabel('current')
    call gp%ylabel('firing rate')
    call gp%plot(currents, st, 'with line lw 2 lc rgb "#0008B0"') 
    
    call gp%title('2')
    call gp%xlabel('current')
    call gp%ylabel('firing rate')
    call gp%plot(currents, nd, 'with line lw 2 lc rgb "#0008B0"') 
    
    call gp%title('3')
    call gp%xlabel('current')
    call gp%ylabel('firing rate')
    call gp%plot(currents, rd, 'with line lw 2 lc rgb "#0008B0"') 
    
    call gp%title('4')
    call gp%xlabel('current')
    call gp%ylabel('firing rate')
    call gp%plot(currents, th, 'with line lw 2 lc rgb "#0008B0"') 
    
end program FxI