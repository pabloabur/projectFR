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
    !use MKL_VSL
    
    implicit none 
    !integer, parameter :: wp = kind(1.0d0)
    type(Configuration) :: conf
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    integer :: timeLength
    integer :: i, j, k
    real(wp), dimension(:), allocatable :: t, dendPotSOL, dendPotMG
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp), dimension(:), allocatable :: FR
    integer :: GammaOrder 
    character(len = 80) :: pool, group
    character(len = 100) :: filename = '../../../conf.rmto'
    character(len = 100) :: path = '/home/pablo/osf/Master-Thesis-Data/population/'
    character(len = 100) :: folderName = 'psd/dc/trial4/'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools     
    type(AfferentPool), dimension(:), allocatable:: afferentPools     
    character(len=2) :: inputTrial
    character(len=4) :: inputParam
    character(len=2) :: inputMVC
    character(len=20) :: gmaxS, gmaxFR, gmaxFF
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
    ! Quantity parameters
    character(len=3), parameter :: nS = '75', nFR = '75', &
        nFF = '150', nCM = '400', nMN = '300', nRC = '600' ! nS+nFR+nFF
    real(wp), dimension(:,:), allocatable :: Vs
    integer :: poolSize = 300 ! Quantify of used MNs

    call init_random_seed()

    !*************************************
    !******* Getting input
    !*************************************
    call get_command_argument(3, inputTrial)
    call get_command_argument(1, inputParam)
    if (inputParam.eq.'o') then
        print *, 'no renshaw cell'
    else if (inputParam.eq.'h') then
        print *, 'renshaw cell with half conductance'
    else if (inputParam.eq.'s') then
        print *, 'renshaw cell with standard conductance'
    else if (inputParam.eq.'d') then
        print *, 'renshaw cell with double conductance'
    else
        print *, 'Wrong parametrization option. Available options are:'
        print *, '- o'
        print *, '- h'
        print *, '- s'
        print *, '- d'
        stop (1)
    endif

    call get_command_argument(2, inputMVC)
    if (inputMVC.eq.'05') then
        print *, '5% MVC'
    else if (inputMVC.eq.'70') then
        print *, '70% MVC'
    else
        print *, 'Wrong MVC option. Available options are:'
        print *, '- 05'
        print *, '- 70'
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

    if (inputMVC.eq.'05') then
        GammaOrder = 7
    else if (inputMVC.eq.'70') then
        GammaOrder = 1
    endif
    ! Removing influence of stimulus (required)
    paramTag = 'stimIntensity_PTN'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    ! Conductances
    if (inputParam.eq.'d') then
        gmaxS = '0.2380'
        gmaxFR = '0.2380'
        gmaxFF = '0.1880'
    else if (inputParam.eq.'s') then
        gmaxS = '0.1190'
        gmaxFR = '0.1190'
        gmaxFF = '0.0940'
    else if (inputParam.eq.'h') then
        gmaxS = '0.0595'
        gmaxFR = '0.0595'
        gmaxFF = '0.0470'
    else if (inputParam.eq.'o') then
        continue
    else
        print *, 'Wrong modulation option.'
        stop (1)
    endif
    paramTag = 'gmax:RC_ext->MG-S@dendrite|inhibitory'
    value1 = gmaxS
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'gmax:RC_ext->MG-FR@dendrite|inhibitory'
    value1 = gmaxFR
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'gmax:RC_ext->MG-FF@dendrite|inhibitory'
    value1 = gmaxFF
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    !*************************************
    !******* Parameters for SOL pools used for computing 
    !******* synaptic conductance caused by commom drive
    !*************************************
    paramTag = 'MUnumber_SOL-S'
    value1 = '300'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'Number_SOL'
    value1 = '300'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'Con:CMExt->SOL-S@dendrite|excitatory'
    value1 = '30'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'Con:RC_ext->SOL-S@dendrite|inhibitory'
    value1 = '0'
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
    !!!!!!!!!!!!!!!! Independent noise
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

    if (inputParam.ne.'o') then
        allocate(interneuronPools(1))
        pool = 'RC'
        group = 'ext'    
        interneuronPools(1) = InterneuronPool(conf, pool, group)
    else
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
    allocate(force(timeLength))
    allocate(dendPotSOL(timeLength))
    allocate(dendPotMG(timeLength))
    allocate(Vs(poolSize, timeLength)) ! For all MNs in simulation

    t = [(dt*(i-1), i=1, timeLength)]
    
    print *, 'Running simulation'
    call cpu_time(tic)

    !*************************************
    !*************** Preparing proper input
    !*************************************
    if (inputParam.eq.'s') then 
        if (inputMVC.eq.'05') then
            FR(:) = 230_wp + 2.3_wp*sin(2*pi*10*t*1e-3)
        else if (inputMVC.eq.'70') then
            FR(:) = 930_wp + 2.3_wp*sin(2*pi*10*t*1e-3)
        end if
    else if (inputParam.eq.'o') then
        if (inputMVC.eq.'05') then
            FR(:) = 141_wp + 2.3_wp*sin(2*pi*10*t*1e-3)
        else if (inputMVC.eq.'70') then
            FR(:) = 665_wp + 2.3_wp*sin(2*pi*10*t*1e-3)
        end if
    else if (inputParam.eq.'d') then
        if (inputMVC.eq.'05') then
            FR(:) = 333_wp + 2.3_wp*sin(2*pi*10*t*1e-3)
        else if (inputMVC.eq.'70') then
            FR(:) = 1180_wp + 2.3_wp*sin(2*pi*10*t*1e-3)
        end if
    else if (inputParam.eq.'h') then
        if (inputMVC.eq.'05') then
            FR(:) = 190_wp + 2.3_wp*sin(2*pi*10*t*1e-3)
        else if (inputMVC.eq.'70') then
            FR(:) = 795_wp + 2.3_wp*sin(2*pi*10*t*1e-3)
        end if
    endif

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

        if (inputParam.ne.'o') then
            do j = 1, size(interneuronPools)
                call interneuronPools(j)%atualizeInterneuronPool(t(i))
            end do
        endif

        do j = 1, size(motorUnitPools)
            call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
            if (trim(motorUnitPools(j)%pool).eq.'MG') then
                do k = 1, poolSize
                    Vs(k,i) = motorUnitPools(j)%v_mV(2*(k))
                end do
                dendPotMG(i) = motorUnitPools(j)%v_mV(2*(1)-1)
            else if (trim(motorUnitPools(j)%pool).eq.'SOL') then
                dendPotSOL(i) = motorUnitPools(j)%v_mV(2*(1)-1)
            end if
        end do
    end do

    !call neuralTractPools(1)%listSpikes()
    call motorUnitPools(1)%listSpikes()
    call motorUnitPools(1)%getMotorUnitPoolEMG()
    if (inputParam.ne.'o') then
        call interneuronPools(1)%listSpikes()
    endif

    do i = 1, timeLength
        force(i) = motorUnitPools(1)%NoHillMuscle%force(i)
    end do


    !*************************************
    !*************** Saving data
    !*************************************
    ! Saving spikes
    folderName = trim(folderName) // 'trial' // trim(inputTrial) // '/'

    filename = trim(path) // trim(folderName) // "spike" // trim(inputMVC) // trim(inputParam) // ".dat"
    open(1, file=filename, status = 'replace')
    do i = 1, size(motorUnitPools(1)%poolSomaSpikes, 1)
        write(1, '(F15.6, 1X, F15.1)') motorUnitPools(1)%poolSomaSpikes(i,1), &
            motorUnitPools(1)%poolSomaSpikes(i,2)
    end do
    close(1)

    ! Saving force
    filename = trim(path) // trim(folderName) // "force" // trim(inputMVC) // trim(inputParam) // ".dat"
    open(1, file=filename, status = 'replace')
    do i = 1, timeLength
           write(1, '(F15.2, 1X, F15.6)') t(i), force(i)
    end do
    close(1)

    ! Saving input conductance and EMG
    filename = trim(path) // trim(folderName) // "g_emg" // trim(inputMVC) // trim(inputParam) // ".dat"
    open(1, file=filename, status = 'replace')
    do i = 1, timeLength
        write(1, '(F15.6, 1X, F15.6)') dendPotSOL(i), motorUnitPools(1)%emg(i)
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

    !call gp%title('Commom drive')
    !call gp%xlabel('t (ms))')
    !call gp%ylabel('rate (pps)')
    !call gp%plot(t, commomDriveG, 'with line lw 2 lc rgb "#0008B0"')

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
