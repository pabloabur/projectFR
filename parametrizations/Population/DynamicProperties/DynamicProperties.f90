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

program DynamicProperties
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
    real(wp) :: dt
    real(wp) :: tf
    integer :: timeLength
    integer :: i, j, k
    real(wp), dimension(:), allocatable :: t, MNv_mV, RCv_mV
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp) :: FR
    integer :: GammaOrder 
    character(len = 80) :: pool, group
    character(len = 100) :: filename = '../../conf.rmto'
    character(len = 100) :: path = '/home/pablo/osf/Master-Thesis-Data/population/'
    character(len = 100) :: folderName = 'dynamic/false_decay/trial10/'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools     
    type(AfferentPool), dimension(:), allocatable:: afferentPools     
    character(len=80) :: paramTag
    character(len=80) :: value1, value2
    ! Input parameters
    logical, parameter :: probDecay = .false.
    real(wp), parameter :: FFConducStrength = 0.033_wp, & 
        declineFactorMN = real(1, wp)/6, declineFactorRC = real(3.5, wp)/3
    character(len=3), parameter :: stimAmp = '80', nS = '75', nFR = '75', &
        nFF = '150', nRC = '600'
    integer, parameter :: frequency1 = 5
    ! values for frequency2 are already substracted to yield proper modulation
    integer, dimension(4), parameter :: frequency2 = [5, 15, 28, 45]

    call init_random_seed()

    conf = Configuration(filename)
    conf%simDuration_ms = 2000

    !Changing configuration file
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

        ! Threshold
        paramTag = 'threshold:RC_ext-'
        value1 = '22.9608'
        value2 = '22.9608'
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
        value1 = '0.128'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:RC_ext->MG-FR@dendrite|inhibitory'
        value1 = '0.119'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:RC_ext->MG-FF@dendrite|inhibitory'
        value1 = '0.094'
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
        value1 = '8500'
        value2 = '8500'
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
    
    ! Stimulus
    paramTag = 'stimStart_PTN'
    value1 = '0.1'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'stimStop_PTN'
    value1 = '2000'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'stimIntensity_PTN'
    value1 = stimAmp
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'stimFrequency_PTN'
    write(value1, '(I1)') frequency1
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'stimPulseDuration_PTN'
    value1 = '0.2'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'stimModulationStart_PTN'
    value1 = '599'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'stimModulationStop_PTN'
    value1 = '2000'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

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

    ! Removing RC spontaneous activity 
    paramtag = 'gmax:Noise>RC_ext-@soma|excitatory'
    value1 = '0.03015'
    value2 = ''
    call conf%changeconfigurationparameter(paramtag, value1, value2)
    paramtag = 'NoiseFunction_RC_ext'
    value1 = '0'! it was '7', but in the experiment there was none
    value2 = ''
    call conf%changeconfigurationparameter(paramtag, value1, value2)

    print *, 'Building neural elements'
    allocate(neuralTractPools(0))
    pool = 'MG'
    allocate(motorUnitPools(1))
    motorUnitPools(1) = MotorUnitPool(conf, pool)    
    pool = 'RC'
    group = 'ext'    
    allocate(interneuronPools(1))
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
    
    print *, 'Running simulation'
    call cpu_time(tic)
    do k=1, size(frequency2)
        paramTag = 'stimModulation_PTN'
        write(value1, '(I15)') frequency2(k)
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Apply stimulus to the nerve
        do i = 1, size(motorUnitPools(1)%unit)
            call motorUnitPools(1)%unit(i)%createStimulus()
        end do
        
        do i = 1, size(t)        
            !do j = 1, size(synapticNoisePools)
            !    call synapticNoisePools(j)%atualizePool(t(i))
            !end do
            do j = 1, size(interneuronPools)
                call interneuronPools(j)%atualizeInterneuronPool(t(i))
                RCv_mV(i) = interneuronPools(j)%v_mV(1)        
            end do
            do j = 1, size(motorUnitPools)
                call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
                MNv_mV(i) = motorUnitPools(j)%v_mV(2*(1))      
            end do
        end do

        call motorUnitPools(1)%listSpikes()
        call interneuronPools(1)%listSpikes()

        write(filename, '("output", I1, ".dat")') k
        filename = trim(path) // trim(folderName) // filename
        open(1, file=filename, status = 'replace')
        do i = 1, size(interneuronPools(1)%poolSomaSpikes, 1)
            write(1, '(F15.6, 1X, F15.1, 1X, F15.2)') interneuronPools(1)%poolSomaSpikes(i,1), &
                interneuronPools(1)%poolSomaSpikes(i,2)
        end do
        close(1)
        write(filename, '("stimulus", I1, ".dat")') k
        filename = trim(path) // trim(folderName) // filename
        open(1, file=filename, status = 'replace')
        do i = 1, size(motorUnitPools(1)%unit(1)%nerveStimulus_mA)
            write(1, '(F15.2)') motorUnitPools(1)%unit(1)%nerveStimulus_mA(i)
        end do
        close(1)

        call motorUnitPools(1)%reset()
        call interneuronPools(1)%reset()
        !call synapticNoisePools(1)%reset()

        !call gp%title('Membrane potential of the soma of the RC #1')
        !call gp%xlabel('t (ms))')
        !call gp%ylabel('Membrane potential (mV)')
        !call gp%plot(t, RCv_mV, 'with line lw 2 lc rgb "#0008B0"') 

        !call gp%title('Membrane potential of the soma of the MN #1')
        !call gp%xlabel('t (ms))')
        !call gp%ylabel('Membrane potential (mV)')
        !call gp%plot(t, MNv_mV, 'with line lw 2 lc rgb "#0008B0"')

        !call gp%title('PTN stimulus')
        !call gp%xlabel('t (ms))')
        !call gp%ylabel('Stimulus (mA)')
        !call gp%plot(t, motorUnitPools(1)%unit(1)%nerveStimulus_mA, 'with line lw 2 lc rgb "#0008B0"')
    end do
    
    call cpu_time(toc)
    print '(F15.6, A)', toc - tic, ' seconds'
    
end program DynamicProperties
