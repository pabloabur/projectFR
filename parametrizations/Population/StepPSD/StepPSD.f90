! '''
!     Neuromuscular simulator in Fortran.
!     Copyright (C) 2018  Renato Naville Watanabe
!                         Pablo Alejandro    
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     any later version.

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.

!     Contact: renato.watanabe@ufabc.edu.br
! '''

program StepPSD
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
    integer :: i, j, k
    real(wp), dimension(:), allocatable :: t, MNDend_mV, MNSoma_mV, RCv_mV, &
        synapticInput, commomDriveG, V
    real(wp), allocatable :: Vs(:,:)
    real(wp) :: sigmaSquared,  sigmaSquaredi,  chiSquared
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp), dimension(:), allocatable :: FR
    !real(wp) :: FR
    integer :: GammaOrder 
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
    real(wp), dimension(:), allocatable :: force
    ! Input parameters
    real(wp) :: dt
    real(wp) :: tf
    logical, parameter :: probDecay = .false.
    real(wp), parameter :: FFConducStrength = 0.0275_wp, & 
        declineFactorMN = real(1, wp)/6, declineFactorRC = real(3.5, wp)/3
    character(len=3), parameter :: nS = '300', nFR = '0', &
        nFF = '0', nCM = '400', nMN = '300' ! nS+nFR+nFF
    character(len=3) :: nRC
    integer :: sampleSize = 101
    integer :: startMNIndex = 100
    integer, parameter :: pulseWidth_ms = 20
    integer, parameter :: frequency_Hz = 10
    real(wp), parameter :: pulsePeriod_ms = real(1, wp)/frequency_Hz*1e3
    integer :: pulseWidth_i, pulsePeriod_i, pulseCycles

    call init_random_seed()

    conf = Configuration(filename)
    conf%simDuration_ms = 6000!30000!100000

    call get_command_argument(1, param)
    if (param.eq.'c') then
        nRC = '600'
    else if (param.eq.'o') then
        nRC = '0'
    else
        print *, 'Wrong parametrization option'
        stop (1)
    endif


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
    value2 = '6'
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'position:RC_ext-'
    value1 = '0'
    value2 = '6'
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    ! ! RC spontaneous activity 
    paramtag = 'gmax:Noise>RC_ext-@soma|excitatory'
    value1 = '0.028'
    value2 = ''
    call conf%changeconfigurationparameter(paramtag, value1, value2)
    paramtag = 'NoiseFunction_RC_ext'
    value1 = '7'
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
    paramTag = 'Con:Noise>MG-@dendrite|excitatory'
    value1 = '100'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'NoiseGammaOrder_MG'
    value1 = '7'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'NoiseTarget_MG'
    value1 = 'FR'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'NoiseFunction_MG'
    value1 = '625'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'gmax:Noise>MG-@dendrite|excitatory'
    value1 = '0.04'
    !value1 = '2.1' ! Compensated value (for closed loop)
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    !!!!!!!!!!!!!!!! Descending commands parameters
    paramTag = 'gmax:CMExt->MG-S@dendrite|excitatory'
    value1 = '0.6'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'gmax:CMExt->MG-FR@dendrite|excitatory'
    value1 = '0.6'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'gmax:CMExt->MG-FF@dendrite|excitatory'
    value1 = '0.6'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    GammaOrder = 7

    ! Removing influence of stimulus (required)
    paramTag = 'stimIntensity_PTN'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    print *, 'Building neural elements'
    allocate(neuralTractPools(1))
    pool = 'CMExt'
    neuralTractPools(1) = NeuralTract(conf, pool)

    allocate(motorUnitPools(2))
    pool = 'MG'
    motorUnitPools(1) = MotorUnitPool(conf, pool)    
    ! MNs from a SOL pool used for computing the 
    ! synaptic conductance caused by commom drive
    ! Parameters were set directly on conf file
    pool = 'SOL'
    motorUnitPools(2) = MotorUnitPool(conf, pool)    

    if (param.eq.'c') then
        allocate(interneuronPools(1))
        pool = 'RC'
        group = 'ext'    
        interneuronPools(1) = InterneuronPool(conf, pool, group)
    else if (param.eq.'o') then
        allocate(interneuronPools(0))
    endif

    allocate(afferentPools(0))

    synapticNoisePools = synapseFactory(conf, neuralTractPools, &
                                        motorUnitPools, &
                                        interneuronPools, &
                                        afferentPools, probDecay)
    
    tf = conf%simDuration_ms
    dt = conf%timeStep_ms
    timeLength = nint(tf/dt)
    
    allocate(t(timeLength))
    allocate(FR(timeLength))
    allocate(MNDend_mV(timeLength))
    allocate(MNSoma_mV(timeLength))
    allocate(V(timeLength))
    allocate(RCv_mV(timeLength))
    allocate(synapticInput(timeLength))
    allocate(commomDriveG(timeLength))
    allocate(force(timeLength))
    allocate(Vs(timeLength,sampleSize)) ! According to slice takenSize

    t = [(dt*(i-1), i=1, timeLength)]

    FR(:) = 0_wp
    pulseCycles = int(tf/pulsePeriod_ms)
    pulsePeriod_i = int(pulsePeriod_ms/dt)
    pulseWidth_i = int(pulseWidth_ms/dt)
    if (pulseWidth_i.ge.pulsePeriod_i) then
        print *, "Pulse width larger or equal pulse period"
        stop (1)
    end if
    do i = 1, pulseCycles
        FR((i-1)*pulsePeriod_i+pulsePeriod_i-pulseWidth_i:pulsePeriod_i*i) = 140_wp
    end do
    
    print *, 'Running simulation'
    call cpu_time(tic)

    do i = 1, size(t)
        do j = 1, size(neuralTractPools)
            call neuralTractPools(j)%atualizePool(t(i), FR(i), GammaOrder)
        end do
        do j = 1, size(synapticNoisePools)
            call synapticNoisePools(j)%atualizePool(t(i))
        end do
        if (param.eq.'c') then
            do j = 1, size(interneuronPools)
                call interneuronPools(j)%atualizeInterneuronPool(t(i))
            !    RCv_mV(i) = interneuronPools(j)%v_mV(1)        
            end do
        endif
        do j = 1, size(motorUnitPools)
            if (j == 1) then
                call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
                !MNDend_mV(i) = motorUnitPools(j)%v_mV(2*(1)-1)
                !MNSoma_mV(i) = motorUnitPools(j)%v_mV(2*(1))
                !synapticInput(i) = motorUnitPools(j)%iIonic(1)
                ! Slice with desired MNs
                do k=startMNIndex, startMNIndex+sampleSize-1
                    Vs(i,k-99) = motorUnitPools(j)%v_mV(2*(k))
                end do
            else
                call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
                MNDend_mV(i) = motorUnitPools(j)%v_mV(2*(1)-1)
                synapticInput(i) = motorUnitPools(j)%iIonic(1)
                commomDriveG(i) = synapticInput(i)/(MNDend_mV(i) - 70)
            end if
        end do
    end do

    call motorUnitPools(1)%listSpikes()
    call motorUnitPools(1)%getMotorUnitPoolEMG()
    !call neuralTractPools(1)%listSpikes()
    !call interneuronPools(1)%listSpikes()

    do i = 1, timeLength
        force(i) = motorUnitPools(1)%NoHillMuscle%force(i)
    end do

    if (param.eq.'c') then
        filename = "inputV_I.dat"
    else if (param.eq.'o') then
        filename = "datNoRC/inputV_I.dat"
    endif
    open(1, file=filename, status = 'replace')
    do i = 1, timeLength
        write(1, '(F15.6, 1X, F15.6)') MNDend_mV(i), synapticInput(i)
    end do
    close(1)

    if (param.eq.'c') then
        filename = "g_emg.dat"
    else if (param.eq.'o') then
        filename = "datNoRC/g_emg.dat"
    endif
    open(1, file=filename, status = 'replace')
    do i = 1, timeLength
        write(1, '(F15.6, 1X, F15.6)') commomDriveG(i), motorUnitPools(1)%emg(i)
    end do
    close(1)

    if (param.eq.'c') then
        filename = "MNspk.dat"
    else if (param.eq.'o') then
        filename = "datNoRC/MNspk.dat"
    endif
    open(1, file=filename, status = 'replace')
    do i = 1, size(motorUnitPools(1)%poolSomaSpikes, 1)
        write(1, '(F15.6, 1X, F15.1)') motorUnitPools(1)%poolSomaSpikes(i,1), &
            motorUnitPools(1)%poolSomaSpikes(i,2)
    end do
    close(1)

    if (param.eq.'c') then
        filename = "force.dat"
    else if (param.eq.'o') then
        filename = "datNoRC/force.dat"
    endif
    open(1, file=filename, status = 'replace')
    do i = 1, timeLength
           write(1, '(F15.2)') force(i)
    end do
    close(1)

    call cpu_time(toc)
    print '(F15.6, A)', toc - tic, ' seconds'
    
    !!!!!!! Measuring synchrony
    ! Synchrony measure computed according to Golomb, Hansel and Mato (2001)
    ! Global part
    V = sum(Vs, dim=2)/sampleSize
    sigmaSquared = sum(V(2000:)**2)/size(t) - (sum(V(2000:))/size(t))**2

    filename = "V.dat"
    open(1, file=filename, status = 'replace')
    do i = 1, timeLength
           write(1, '(F15.2)') V(i)
    end do
    close(1)

    ! Individual part
    sigmaSquaredi = 0
    do i=1, sampleSize
        sigmaSquaredi = sigmaSquaredi + sum(Vs(2000:,i)**2)/size(t)-&
            (sum(Vs(2000:,i))/size(t))**2
    end do
    sigmaSquaredi = sigmaSquaredi/sampleSize

    ! Chi squared
    chiSquared = sigmaSquared/sigmaSquaredi
    print *, '----------- Chi --------------'
    print *, sqrt(chiSquared)
    print *, '------------------------------'

    !call gp%title('MN1')
    !call gp%xlabel('t (ms))')
    !call gp%ylabel('Volts (mV)')
    !call gp%plot(t, FR, 'with line lw 2 lc rgb "#0008B0"') 

    !call gp%title('Voltage at soma of MN 1')
    !call gp%xlabel('t (ms))')
    !call gp%ylabel('Volts (mV)')
    !call gp%plot(t, MNSoma_mV, 'with line lw 2 lc rgb "#0008B0"') 

    !call gp%title('MN spike instants at the soma')
    !call gp%xlabel('t (s))')
    !call gp%ylabel('Motoneuron index')
    !call gp%plot(motorUnitPools(1)%poolSomaSpikes(:,1), &
    !motorUnitPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

end program StepPSD
