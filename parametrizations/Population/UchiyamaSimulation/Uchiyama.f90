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

program ForcedPSD
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
    integer :: timeLength
    integer :: i, j
    real(wp), dimension(:), allocatable :: t
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp) :: FR(41) = [(i, i = 0, 400, 10)]
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
    ! Input parameters
    real(wp) :: dt
    real(wp) :: tf
    logical, parameter :: probDecay = .false.
    real(wp), parameter :: FFConducStrength = 0.0275_wp, & 
        declineFactorMN = real(1, wp)/6, declineFactorRC = real(3.5, wp)/3
    character(len=3), parameter :: nS = '75', nFR = '75', &
        nFF = '150', nRC = '600', nCM = '400', nMN = '300' ! nS+nFR+nFF
    integer :: k

    call init_random_seed()

    conf = Configuration(filename)
    conf%simDuration_ms = 1000

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
    else
        print *, 'Wrong parametrization option'
        stop (1)
    endif

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
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'gmax:Noise>MG-@dendrite|excitatory'
    value1 = '10'
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
    GammaOrder = 10

    ! Removing influence of stimulus (required)
    paramTag = 'stimIntensity_PTN'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    print *, 'Building neural elements'
    allocate(neuralTractPools(1))
    pool = 'CMExt'
    neuralTractPools(1) = NeuralTract(conf, pool)

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

    t = [(dt*(i-1), i=1, timeLength)]
    
    print *, 'Running simulation'
    call cpu_time(tic)

    do k=1, size(FR)
        do i = 1, size(t)
            do j = 1, size(neuralTractPools)
                call neuralTractPools(j)%atualizePool(t(i), FR(k), GammaOrder)
            end do
            do j = 1, size(motorUnitPools)
                call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
            end do
            do j = 1, size(synapticNoisePools)
                call synapticNoisePools(j)%atualizePool(t(i))
            end do
            do j = 1, size(interneuronPools)
                call interneuronPools(j)%atualizeInterneuronPool(t(i))
            end do
        end do

        call motorUnitPools(1)%listSpikes()
        !call neuralTractPools(1)%listSpikes() Allocation would exceed memory limit
        call interneuronPools(1)%listSpikes()

        write(filename, '("spks", I2, ".dat")') k
        open(1, file=filename, status = 'replace')
        do i = 1, size(motorUnitPools(1)%poolSomaSpikes, 1)
            write(1, '(F15.6, 1X, F15.1)') motorUnitPools(1)%poolSomaSpikes(i,1), &
                motorUnitPools(1)%poolSomaSpikes(i,2)
        end do
        write(1, '(A2)') 'RC'
        close(1)

        open(1, file=filename, position = 'append', status = 'old')
        do i = 1, size(interneuronPools(1)%poolSomaSpikes, 1)
            write(1, '(F15.6, 1X, F15.1)') interneuronPools(1)%poolSomaSpikes(i,1), &
                interneuronPools(1)%poolSomaSpikes(i,2)
        end do
        close(1)

        !call gp%title('MN spike instants at the soma')
        !call gp%xlabel('t (s))')
        !call gp%ylabel('Motoneuron index')
        !call gp%plot(motorUnitPools(1)%poolSomaSpikes(:,1), &
        !motorUnitPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

        !call gp%title('RC spike instants at the soma')
        !call gp%xlabel('t (s))')
        !call gp%ylabel('Interneuron index')
        !call gp%plot(interneuronPools(1)%poolSomaSpikes(:,1), &
        !interneuronPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

        call neuralTractPools(1)%reset()
        call motorUnitPools(1)%reset()
        call interneuronPools(1)%reset()
    end do

    call cpu_time(toc)
    print '(F15.6, A)', toc - tic, ' seconds'

end program ForcedPSD