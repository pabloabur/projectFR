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
program FrequencyAnalysis
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
    real(wp), dimension(:), allocatable :: t, dendPotSOL, dendPotMG, &
        synCurrentSOL, synCurrentMG, synConductSOL, synConductMG
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp), dimension(:), allocatable :: FR
    integer :: GammaOrder 
    character(len = 80) :: pool, group
    character(len = 100) :: filename = '../../../conf.rmto'
    character(len = 100) :: path = '/home/pablo/osf/Master-Thesis-Data/population/'
    character(len = 100) :: folderName = 'psd/cancel/trial1/'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools     
    type(AfferentPool), dimension(:), allocatable:: afferentPools     
    character(len=2) :: inputTrial
    character(len=2) :: inputParam
    character(len=80) :: paramTag
    character(len=80) :: value1, value2
    real(wp), dimension(:), allocatable :: force
    real(wp) :: dt
    real(wp) :: tf
    logical, parameter :: probDecay = .false.
    ! Noise parameters
    real(wp) :: noiseFR
    character(len=3), parameter :: noiseCon = '100'
    character(len=1), parameter :: noiseOrder = '1'
    character(len=2), parameter :: noiseTarget = 'FR'
    character(len=2), parameter :: noiseg = '12' !  Step simulation
    real(wp) :: Esyn
    character(len=3), parameter :: nS = '75', nFR = '75', &
        nFF = '150', nCM = '400', nMN = '300', nRC = '600' ! nS+nFR+nFF
    real(wp), dimension(:,:), allocatable :: Vs
    integer :: poolSize = 300 ! Quantify of used MNs

    call init_random_seed()

    !*************************************
    !******* Getting input
    !*************************************
    call get_command_argument(2, inputTrial)
    call get_command_argument(1, inputParam)
    if (inputParam.eq.'dc') then
        print *, 'descending command only'
    else if (inputParam.eq.'rc') then
        print *, 'renshaw cell effect only'
    else
        print *, 'Wrong parametrization option. Available options are:'
        print *, '- dc'
        print *, '- rc'
        stop (1)
    endif

    !*************************************
    !******* Changing basic parameters
    !*************************************
    conf = Configuration(filename)
    conf%simDuration_ms = 11000

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


    !!!!!!!!!!!!!!!! Independent noise
    ! N.B. noiseFunc is not being used anymore, it is passed as argument now
    paramTag = 'Con:Noise>MG-@dendrite|excitatory'
    value1 = noiseCon
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'NoiseGammaOrder_MG'
    value1 = noiseOrder
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'NoiseTarget_MG'
    value1 = noiseTarget
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'gmax:Noise>MG-@dendrite|excitatory'
    value1 = noiseg
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    GammaOrder = 7
    ! Removing influence of stimulus (required)
    paramTag = 'stimIntensity_PTN'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'gmax:RC_ext->SOL-S@dendrite|inhibitory'
    value1 = '0.1190'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'MUnumber_SOL-S'
    value1 = '300'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'Number_SOL'
    value1 = '300'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'Con:SOL-S>RC_ext-@soma|excitatory'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'threshold:SOL-S'
    value1 = '500'
    value2 = '500'
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'position:SOL-'
    value1 = '0'
    value2 = '6'
    call conf%changeConfigurationParameter(paramTag, value1, value2)


    if (inputParam.eq.'dc') then
        paramTag = 'Con:CMExt->SOL-S@dendrite|excitatory'
        value1 = '30'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:RC_ext->SOL-S@dendrite|inhibitory'
        value1 = '0'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:Noise>SOL-@dendrite|excitatory'
        value1 = noiseCon
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'NoiseGammaOrder_SOL'
        value1 = noiseOrder
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'NoiseTarget_SOL'
        value1 = noiseTarget
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'gmax:Noise>SOL-@dendrite|excitatory'
        value1 = noiseg
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        Esyn = 70._wp
    else if (inputParam.eq.'rc') then
        paramTag = 'Con:CMExt->SOL-S@dendrite|excitatory'
        value1 = '0'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:RC_ext->SOL-S@dendrite|inhibitory'
        value1 = '100'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramTag = 'Con:Noise>SOL-@dendrite|excitatory'
        value1 = '0'
        value2 = ''
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        Esyn = -16._wp
    endif

    print *, '*************************************'
    print *, '************ Building neural elements'
    print *, '*************************************'
    allocate(neuralTractPools(1))
    pool = 'CMExt'
    neuralTractPools(1) = NeuralTract(conf, pool)

    allocate(motorUnitPools(2))
    pool = 'MG'
    motorUnitPools(1) = MotorUnitPool(conf, pool)    
    pool = 'SOL'
    motorUnitPools(2) = MotorUnitPool(conf, pool)    

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
    allocate(FR(timeLength))
    allocate(force(timeLength))
    allocate(dendPotSOL(timeLength))
    allocate(dendPotMG(timeLength))
    allocate(synCurrentSOL(timeLength))
    allocate(synCurrentMG(timeLength))
    allocate(synConductSOL(timeLength))
    allocate(synConductMG(timeLength))
    allocate(Vs(poolSize, timeLength)) ! For all MNs in simulation

    t = [(dt*(i-1), i=1, timeLength)]
    
    print *, 'Running simulation'
    call cpu_time(tic)

    !*************************************
    !*************** Preparing proper input
    !*************************************
    FR(:) = 230_wp + 11.5_wp*sin(2*pi*10*t*1e-3)
    noiseFR = 90_wp

    !*************************************
    !*************** Running simulation
    !*************************************
    do i = 1, size(t)        
        ! Updating elements
        do j = 1, size(neuralTractPools)
            call neuralTractPools(j)%atualizePool(t(i), FR(i), GammaOrder)
        end do

        do j = 1, size(synapticNoisePools)
            if (synapticNoisePools(j)%pool.eq.'MG') then
                call synapticNoisePools(j)%atualizePool(t(i), noiseFR)
            else if (synapticNoisePools(j)%pool.eq.'SOL') then
                call synapticNoisePools(j)%atualizePool(t(i), noiseFR)
            else if (synapticNoisePools(j)%pool.eq.'RC_ext') then
                call synapticNoisePools(j)%atualizePool(t(i), 7.0_wp)
            else
                print *, 'Error assigning noise value to pool'
                stop (1)
            endif
        end do

        do j = 1, size(interneuronPools)
            call interneuronPools(j)%atualizeInterneuronPool(t(i))
        end do

        do j = 1, size(motorUnitPools)
            call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
            if (trim(motorUnitPools(j)%pool).eq.'MG') then
                do k = 1, poolSize
                    Vs(k,i) = motorUnitPools(j)%v_mV(2*(k))
                end do
                dendPotMG(i) = motorUnitPools(j)%v_mV(2*(1)-1)
                synCurrentMG(i) = motorUnitPools(j)%iIonic(1)
                synConductMG(i) = synCurrentMG(i)/(dendPotMG(i) - Esyn)
            else if (trim(motorUnitPools(j)%pool).eq.'SOL') then
                dendPotSOL(i) = motorUnitPools(j)%v_mV(2*(40)-1)
                synCurrentSOL(i) = motorUnitPools(j)%iIonic(2*(40)-1)
                synConductSOL(i) = synCurrentSOL(i)/(dendPotSOL(i) - Esyn)
            end if
        end do
    end do

    !*************************************
    !*************** Saving data
    !*************************************
    folderName = trim(folderName) // 'trial' // trim(inputTrial) // '/'

    ! Saving input conductance and EMG
    filename = trim(path) // trim(folderName) // "gsyn" // trim(inputParam) // ".dat"
    open(1, file=filename, status = 'replace')
    do i = 1, timeLength
        write(1, '(F15.2, 1X, F15.6, 1X, F15.6)') t(i), synConductMG(i), synConductSOL(i)
    end do
    close(1)

    !*************************************
    !*************** Plotting
    !*************************************
    !if (inputParam.eq.'s') then
    !    call gp%title('RC spike instants at the soma')
    !    call gp%xlabel('t (s))')
    !    call gp%ylabel('Interneuron index')
    !    call gp%plot(interneuronPools(1)%poolSomaSpikes(:,1), &
    !    interneuronPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')
    !end if

    !call gp%title('Descending command')
    !call gp%xlabel('t (ms))')
    !call gp%ylabel('descending command rate (pps)')
    !call gp%plot(t, FR, 'with line lw 2 lc rgb "#0008B0"')

    !call gp%title('Independent noise')
    !call gp%xlabel('t (ms))')
    !call gp%ylabel('command rate (pps)')
    !call gp%plot(t, noiseFR, 'with line lw 2 lc rgb "#0008B0"')

    !call gp%title('MN spike instants at the soma')
    !call gp%xlabel('t (s))')
    !call gp%ylabel('Motoneuron index')
    !call gp%plot(motorUnitPools(1)%poolSomaSpikes(:,1), &
    !motorUnitPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    !call gp%title('Soma voltage of MN 1')
    !call gp%xlabel('t (ms))')
    !call gp%ylabel('Volts (mV)')
    !call gp%plot(t, Vs(1,:), 'with line lw 2 lc rgb "#0008B0"')

    !call gp%title('force')
    !call gp%xlabel('t (ms))')
    !call gp%ylabel('force (N)')
    !call gp%plot(t, motorUnitPools(1)%NoHillMuscle%force, 'with line lw 2 lc rgb "#0008B0"')

    !*************************************
    !*************** Reset elements
    !*************************************
    !if (inputParam.eq.'s') then
    !    call interneuronPools(1)%reset()
    !endif
    !do j = 1, size(motorUnitPools)
    !    call motorUnitPools(j)%reset()
    !end do
    !call neuralTractPools(1)%reset()
    !do j = 1, size(synapticNoisePools)
    !    call synapticNoisePools(j)%reset()
    !end do

    call cpu_time(toc)
    print '(F15.6, A)', toc - tic, ' seconds'
    
end program FrequencyAnalysis
