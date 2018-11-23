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

program StaticProperties
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
    real(wp) :: dt
    real(wp) :: tf
    integer :: timeLength
    integer :: i, j, k
    real(wp), dimension(:), allocatable :: t, MNv_mV, RCv_mV
    real(wp) :: tic, toc
    type(gpf) :: gp
    character(len = 80) :: pool, group
    character(len = 100) :: filename = '../../conf.rmto'
    character(len = 100) :: path = '/home/pablo/osf/Master-Thesis-Data/population/'
    character(len = 100) :: folderName = 'Static/false_decay/trial1/'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools     
    type(AfferentPool), dimension(:), allocatable:: afferentPools     
    character(len=5) :: param
    character(len=80) :: paramTag
    character(len=80) :: value1, value2
    ! Input parameters
    logical, parameter :: probDecay = .false.
    real(wp), parameter :: declineFactorMN = real(1, wp)/6, declineFactorRC = real(3.5, wp)/3
    real(wp), parameter :: FFConducStrength = 0.033_wp ! Simple weight decay
    !real(wp), parameter :: FFConducStrength = 0.025_wp ! Decay of probability of connection (Pconn) with distance
    !real(wp), parameter :: FFConducStrength = 0.04_wp ! Decay of weight and Pconn with distance
    !real(wp), parameter :: FFConducStrength = 0.055_wp ! Sparse connections and weight decay
    !real(wp), parameter :: FFConducStrength = 0.06_wp ! Sparse connections, weight and Pcon decay
    character(len=3), parameter :: stimAmp = '80', nS = '75', nFR = '75', &
        nFF = '150', nRC = '600'
    integer, dimension(7), parameter :: fs = [10, 20, 30, 40, 50, 60, 70]

    call init_random_seed()

    conf = Configuration(filename)

    conf%simDuration_ms = 1000

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
        value1 = '22.9608'
        value2 = '22.9608'
        call conf%changeConfigurationParameter(paramTag, value1, value2)

        ! Connectivity
        paramTag = 'Con:RC_ext->MG-S@dendrite|inhibitory'
        value1 = '100'!70
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:RC_ext->MG-FR@dendrite|inhibitory'
        value1 = '100'!70
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:RC_ext->MG-FF@dendrite|inhibitory'
        value1 = '100'!'70'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:MG-S>RC_ext-@soma|excitatory'
        value1 = '100'!'41'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:MG-FR>RC_ext-@soma|excitatory'
        value1 = '100'!'41'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:MG-FF>RC_ext-@soma|excitatory'
        value1 = '100'!'41'
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

        ! RC spontaneous activity
        paramtag = 'gmax:Noise>RC_ext-@soma|excitatory'
        value1 = '0.03015'
        value2 = ''
        call conf%changeconfigurationparameter(paramtag, value1, value2)
        paramtag = 'NoiseFunction_RC_ext'
        value1 = '0'! it was '7', but in the experiment there was none
        value2 = ''
        call conf%changeconfigurationparameter(paramtag, value1, value2)
    else 
        print *, 'Wrong parametrization option'
        stop (1)
    endif
    
    ! Stimulus
    paramTag = 'stimStart_PTN'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'stimStop_PTN'
    value1 = '1000'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'stimIntensity_PTN'
    value1 = stimAmp
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
            !motorUnitPools(1)%iInjected(2*(21)) = 70.0
            !motorUnitPools(1)%iInjected(2*(22)) = 70.0
            !motorUnitPools(1)%iInjected(2*(23)) = 70.0
            !motorUnitPools(1)%iInjected(2*(24)) = 70.0
            !motorUnitPools(1)%iInjected(2*(25)) = 70.0
            !motorUnitPools(1)%iInjected(:) = 70.0
            do j = 1, size(interneuronPools)
                call interneuronPools(j)%atualizeInterneuronPool(t(i))
                !RCv_mV(i) = interneuronPools(j)%v_mV(301)
            end do
            do j = 1, size(motorUnitPools)
                call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
                !MNv_mV(i) = motorUnitPools(j)%v_mV(2*(1))
            end do
        end do

        call motorUnitPools(1)%listSpikes()
        call interneuronPools(1)%listSpikes()

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

        !call gp%title('PTN stimulus')
        !call gp%xlabel('t (ms))')
        !call gp%ylabel('Stimulus (mA)')
        !call gp%plot(t, motorUnitPools(1)%unit(1)%nerveStimulus_mA, 'with line lw 2 lc rgb "#0008B0"')  

        !call gp%title('Membrane potential of the soma of the MN #1')
        !call gp%xlabel('t (ms))')
        !call gp%ylabel('Descending command index')
        !call gp%plot(t, MNv_mV, 'with line lw 2 lc rgb "#0008B0"')  

        !call gp%title('Membrane potential of the soma of the RC')
        !call gp%xlabel('t (ms))')
        !call gp%ylabel('Potential (mV)')
        !call gp%plot(t, RCv_mV, 'with line lw 2 lc rgb "#0008B0"')  

        write(filename, '("output", I2, ".dat")') fs(k)
        filename = trim(path) // trim(folderName) // filename
        open(1, file=filename, status = 'replace')
        do i = 1, size(interneuronPools(1)%poolSomaSpikes, 1)
            write(1, '(F15.6, 1X, F15.1)') interneuronPools(1)%poolSomaSpikes(i,1), &
                interneuronPools(1)%poolSomaSpikes(i,2)
        end do
        close(1)

        call motorUnitPools(1)%reset()
        call interneuronPools(1)%reset()
    end do
    
    call cpu_time(toc)
    print '(F15.6, A)', toc - tic, ' seconds'
    
end program StaticProperties
