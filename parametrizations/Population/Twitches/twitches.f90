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

program twitches
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
    real(wp), dimension(:), allocatable :: t, MNSoma_mV
    real(wp) :: tic, toc
    type(gpf) :: gp
    integer :: GammaOrder 
    character(len = 80) :: pool, group
    character(len = 100) :: filename = '../../conf.rmto'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools     
    type(AfferentPool), dimension(:), allocatable:: afferentPools     
    character(len=80) :: paramTag
    character(len=80) :: value1, value2
    real(wp), dimension(:), allocatable :: force
    ! Input parameters
    real(wp) :: dt
    real(wp) :: tf
    logical, parameter :: probDecay = .false.
    character(len=3), parameter :: nS = '75', nFR = '75', &
        nFF = '150', nCM = '400', nMN = '300' ! nS+nFR+nFF
    character(len=4), parameter :: nRC = '600'

    call init_random_seed()

    conf = Configuration(filename)
    conf%simDuration_ms = 1000

    !Changing configuration file
    paramTag = 'Number_CMExt'
    value1 = nCM
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'MUnumber_SOL-S'
    value1 = nS
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'MUnumber_SOL-FR'
    value1 = nFR
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'MUnumber_SOL-FF'
    value1 = nFF
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'Number_RC_ext'
    value1 = nRC
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'Number_SOL'
    value1 = nMN
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'Number_RC_ext'
    value1 = nRC
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Connectivity
        paramTag = 'Con:RC_ext->SOL-S@dendrite|inhibitory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:RC_ext->SOL-FR@dendrite|inhibitory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:RC_ext->SOL-FF@dendrite|inhibitory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:SOL-S>RC_ext-@soma|excitatory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:SOL-FR>RC_ext-@soma|excitatory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:SOL-FF>RC_ext-@soma|excitatory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Conductances
        paramTag = 'gmax:RC_ext->SOL-S@dendrite|inhibitory'
        value1 = '0.119'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:RC_ext->SOL-FR@dendrite|inhibitory'
        value1 = '0.119'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:RC_ext->SOL-FF@dendrite|inhibitory'
        value1 = '0.094'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:SOL-S>RC_ext-@soma|excitatory'
        value1 = '0.00681'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:SOL-FR>RC_ext-@soma|excitatory'
        value1 = '0.00833'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:SOL-FF>RC_ext-@soma|excitatory'
        value1 = '0.01500'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Decay factors
        paramTag = 'dec:SOL-S>RC_ext-@soma|excitatory'
        value1 = '0.166666'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'dec:SOL-FR>RC_ext-@soma|excitatory'
        value1 = '0.166666'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'dec:SOL-FF>RC_ext-@soma|excitatory'
        value1 = '0.166666'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'dec:RC_ext->SOL-S@dendrite|inhibitory'
        value1 = '1.166666'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'dec:RC_ext->SOL-FR@dendrite|inhibitory'
        value1 = '1.166666'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'dec:RC_ext->SOL-FF@dendrite|inhibitory'
        value1 = '1.166666'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Twitches
        paramTag = 'twitchTimePeak:SOL-S'
        value1 = '110'
        value2 = '58'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'twitchTimePeak:SOL-FR'
        value1 = '55'
        value2 = '30'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'twitchTimePeak:SOL-FF'
        value1 = '47'
        value2 = '20'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'twitchPeak:SOL-S'
        value1 = '0.002'
        value2 = '0.053'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'twitchPeak:SOL-FR'
        value1 = '0.002'
        value2 = '0.074'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'twitchPeak:SOL-FF'
        value1 = '0.012'
        value2 = '0.380'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'twTetSOCDS:SOL-S'
        value1 = '8.5'
        value2 = '7.4'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'twTetSOCDS:SOL-FR'
        value1 = '18.5'
        value2 = '7.7'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'twTetSOCDS:SOL-FF'
        value1 = '2.5'
        value2 = '2.7'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'bSatSOCDS:SOL-S'
        value1 = '0.24'
        value2 = '0.27'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'bSatSOCDS:SOL-FR'
        value1 = '0.1'
        value2 = '0.25'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'bSatSOCDS:SOL-FF'
        value1 = '0.84'
        value2 = '0.77'
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Axon conductions
        paramTag = 'axonDelayCondVel:SOL-S'
        value1 = '64'
        value2 = '102'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'axonDelayCondVel:SOL-FR'
        value1 = '82'
        value2 = '114'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'axonDelayCondVel:SOL-FF'
        value1 = '74'
        value2 = '122'
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Some parameters were very high. Probably need it for another simulation
        paramTag = 'threshold:SOL-S'
        value1 = '12.35'
        value2 = '16.45'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:CMExt->SOL-S@dendrite|excitatory'
        value1 = '30'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)

    !!!!!!!!!!!!!!!! Turn off independent noise on MNs
    paramTag = 'NoiseTarget_SOL'
    value1 = 'FR'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'NoiseFunction_SOL'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    GammaOrder = 1

    !! Stimulus
    !paramTag = 'stimStart_PTN'
    !value1 = '0'
    !value2 = ''
    !call conf%changeConfigurationParameter(paramTag, value1, value2)
    !paramTag = 'stimStop_PTN'
    !value1 = '4000'
    !value2 = ''
    !call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'stimIntensity_PTN'
    value1 = '0' ! 0 to turn off, 80 for supramaximal
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    !paramTag = 'stimPulseDuration_PTN'
    !value1 = '0.2'
    !value2 = ''
    !call conf%changeConfigurationParameter(paramTag, value1, value2)
    !paramTag = 'stimModulationStart_PTN'
    !value1 = '0'
    !value2 = ''
    !call conf%changeConfigurationParameter(paramTag, value1, value2)
    !paramTag = 'stimModulationStop_PTN'
    !value1 = '0'
    !value2 = ''
    !call conf%changeConfigurationParameter(paramTag, value1, value2)

    print *, 'Building neural elements'
    allocate(neuralTractPools(0))

    allocate(motorUnitPools(1))
    pool = 'SOL'
    motorUnitPools(1) = MotorUnitPool(conf, pool)    

    allocate(interneuronPools(0))
    !pool = 'RC'
    !group = 'ext'    
    !interneuronPools(1) = InterneuronPool(conf, pool, group)

    allocate(afferentPools(0))

    synapticNoisePools = synapseFactory(conf, neuralTractPools, &
                                        motorUnitPools, &
                                        interneuronPools, &
                                        afferentPools, probDecay)
    
    tf = conf%simDuration_ms
    dt = conf%timeStep_ms
    timeLength = nint(tf/dt)
    
    allocate(t(timeLength))
    allocate(MNSoma_mV(timeLength))
    allocate(force(timeLength))

    t = [(dt*(i-1), i=1, timeLength)]
    
    print *, 'Running simulation'
    call cpu_time(tic)

    do i = 1, size(t)
        if (t(i)>2.and.t(i)<1000) then ! 3 and 1000
            motorUnitPools(1)%iInjected(2*(1)) = 85_wp
        else
            motorUnitPools(1)%iInjected(2*(1)) = 0_wp
        endif
        do j = 1, size(motorUnitPools)
            call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
            MNSoma_mV(i) = motorUnitPools(j)%v_mV(2*(1))
        end do
    end do

    call motorUnitPools(1)%listSpikes()

    do i = 1, timeLength
        force(i) = motorUnitPools(1)%NoHillMuscle%force(i)
    end do

    call gp%title('MN spike instants at the soma')
    call gp%xlabel('t (s))')
    call gp%ylabel('Motoneuron index')
    call gp%plot(motorUnitPools(1)%poolSomaSpikes(:,1), &
    motorUnitPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    !call gp%title('RC spike instants at the soma')
    !call gp%xlabel('t (s))')
    !call gp%ylabel('Interneuron index')
    !call gp%plot(interneuronPools(1)%poolSomaSpikes(:,1), &
    !interneuronPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    call gp%title('force')
    call gp%xlabel('t (ms))')
    call gp%ylabel('force (N)')
    call gp%plot(t, force, 'with line lw 2 lc rgb "#0008B0"') 

    !call gp%title('membrane')
    !call gp%xlabel('t (ms))')
    !call gp%ylabel('volts (mV)')
    !call gp%plot(t, MNSoma_mV, 'with line lw 2 lc rgb "#0008B0"') 

    call cpu_time(toc)
    print '(F15.6, A)', toc - tic, ' seconds'

end program twitches
