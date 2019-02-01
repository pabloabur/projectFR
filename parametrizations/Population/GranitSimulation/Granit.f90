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

program Granit
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
    integer :: i, j, k, l
    real(wp), dimension(:), allocatable :: t
    real(wp) :: tic, toc
    !real(wp) :: FR(15) = [(i, i = 5, 537, 38)]
    integer, dimension(7), parameter :: fs = [10, 20, 30, 40, 50, 60, 70]
    integer :: GammaOrder 
    character(len = 80) :: pool, group
    character(len = 100) :: filename = '../../conf.rmto'
    character(len = 100) :: path = '/home/pablo/osf/Master-Thesis-Data/population/'
    character(len = 100) :: folderName = 'granit/trial1/'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools     
    type(AfferentPool), dimension(:), allocatable:: afferentPools     
    character(len=1), dimension(2) :: param
    character(len=80) :: paramTag
    character(len=80) :: value1, value2
    real(wp), dimension(:), allocatable :: force
    ! Input parameters
    real(wp) :: dt
    real(wp) :: tf
    logical, parameter :: probDecay = .false.
    character(len=3), parameter :: nS = '500', nFR = '0', &
        nFF = '0', nCM = '400', nMN = '500'! nS+nFR+nFF
    character(len=4), parameter :: nRC = '1000' 

    call init_random_seed()

    conf = Configuration(filename)
    conf%simDuration_ms = 4000

    param = ['c', 'o']

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
        value1 = '0.076'
        value2 = '0.8'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'twitchPeak:SOL-FR'
        value1 = '0.06'
        value2 = '0.74'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'twitchPeak:MG-FF'
        value1 = '3.35'
        value2 = '14.52'
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

        ! Columnar length
        paramTag = 'position:SOL-'
        value1 = '0'
        value2 = '10'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'position:RC_ext-'
        value1 = '0'
        value2 = '10'
        call conf%changeConfigurationParameter(paramTag, value1, value2)

    ! Dynamics of MN-RC synapse
    paramtag = 'dyn:SOL-S>RC_ext-@soma|excitatory'
    value1 = 'None'
    value2 = ''
    call conf%changeconfigurationparameter(paramtag, value1, value2)
    paramTag = 'dyn:SOL-FR>RC_ext-@soma|excitatory'
    value1 = 'None'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'dyn:SOL-FF>RC_ext-@soma|excitatory'
    value1 = 'None'
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

    ! Stimulus
    paramTag = 'stimStart_PTN'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'stimStop_PTN'
    value1 = '4000'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'stimIntensity_PTN'
    value1 = '80' ! 0 to turn off
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'stimPulseDuration_PTN'
    value1 = '0.2'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'stimModulationStart_PTN'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'stimModulationStop_PTN'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    print *, 'Building neural elements'
    allocate(neuralTractPools(0))
    !pool = 'CMExt'
    !neuralTractPools(1) = NeuralTract(conf, pool)

    allocate(motorUnitPools(1))
    pool = 'SOL'
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
    allocate(force(timeLength))

    t = [(dt*(i-1), i=1, timeLength)]
    
    print *, 'Running simulation'
    call cpu_time(tic)

    do l = 1, size(param)
        do k=1, size(fs)
            paramTag = 'stimFrequency_PTN'
            write(value1, '(I15)') fs(k)
            value2 = ''
            call conf%changeConfigurationParameter(paramTag, value1, value2)

            ! Apply stimulus to the nerve
            do i = 1, size(motorUnitPools(1)%unit)
                call motorUnitPools(1)%unit(i)%createStimulus()
            end do
        
            do i = 1, size(t)
                do j = 1, size(motorUnitPools)
                    call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
                end do
                do j = 1, size(synapticNoisePools)
                    call synapticNoisePools(j)%atualizePool(t(i))
                end do
                if (param(l).eq.'c') then
                    do j = 1, size(interneuronPools)
                        call interneuronPools(j)%atualizeInterneuronPool(t(i))
                    end do
                endif
            end do

            call motorUnitPools(1)%listSpikes()
            if (param(l).eq.'c') then
                call interneuronPools(1)%listSpikes()
            end if

            do i = 1, timeLength
                force(i) = motorUnitPools(1)%NoHillMuscle%force(i)
            end do

            if (param(l).eq.'c') then
                write(filename, '("force", I2, "c.dat")') k
                filename = trim(path) // trim(folderName) // filename
            else if (param(l).eq.'o') then
                write(filename, '("force", I2, "o.dat")') k
                filename = trim(path) // trim(folderName) // filename
            endif
            open(1, file=filename, status = 'replace')
            do i = 1, timeLength
                write(1, '(F9.3, 1X, F15.2)') t(i), force(i)
            end do
            close(1)

            if (param(l).eq.'c') then
                write(filename, '("spks", I2, "c.dat")') k
                filename = trim(path) // trim(folderName) // filename
            else if (param(l).eq.'o') then
                write(filename, '("spks", I2, "o.dat")') k
                filename = trim(path) // trim(folderName) // filename
            endif
            open(1, file=filename, status = 'replace')
            do i = 1, size(motorUnitPools(1)%poolSomaSpikes, 1)
                write(1, '(F15.6, 1X, F15.1)') motorUnitPools(1)%poolSomaSpikes(i,1), &
                    motorUnitPools(1)%poolSomaSpikes(i,2)
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

            if (param(l).eq.'c') then
                call interneuronPools(1)%reset()
                call synapticNoisePools(2)%reset()
            endif
            call motorUnitPools(1)%reset()
            call synapticNoisePools(1)%reset()
        end do
    end do

    call cpu_time(toc)
    print '(F15.6, A)', toc - tic, ' seconds'

end program Granit
