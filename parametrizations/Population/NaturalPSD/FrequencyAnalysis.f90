! '''
!   README: This program performs two distinct simulations: noise 
!   and descending command strategy. The former basically tries to 
!   reproduce results from Williams and Baker (2009), while the latter
!   is a reinterpretation of it. Commented parts are used to choose
!   between simulations. Note that noise strategy only uses s and o
!   cases, both for 5% MVC.
! '''
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
    real(wp), dimension(:), allocatable :: t, &
        synapticInput, MNDend_mV, commomDriveG
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp), dimension(:), allocatable :: FR
    integer :: GammaOrder 
    character(len = 80) :: pool, group
    character(len = 100) :: filename = '../../conf.rmto'
    character(len = 100) :: path = '/home/pablo/osf/Master-Thesis-Data/population/'
    character(len = 100) :: folderName = 'psd/natural/trial4/'
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
    integer, dimension(:), allocatable :: firingMNsIdx, auxFiring
    ! Quantity parameters
    character(len=3), parameter :: nS = '75', nFR = '75', &
        nFF = '150', nCM = '400', nMN = '300', nRC = '600' ! nS+nFR+nFF
    real(wp), dimension(:), allocatable :: V
    real(wp), dimension(:,:), allocatable :: Vs, Vsf
    real(wp) :: delta, sumDeltai, sync
    integer :: sampleSize ! Quantify of used MNs
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
    else if (inputMVC.eq.'30') then
        print *, '30% MVC'
    else if (inputMVC.eq.'70') then
        print *, '70% MVC'
    else
        print *, 'Wrong MVC option. Available options are:'
        print *, '- 05'
        print *, '- 30'
        print *, '- 70'
        stop (1)
    endif
    !*************************************
    !******* Changing basic parameters
    !*************************************
    conf = Configuration(filename)
    conf%simDuration_ms = 9000

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

    ! I gave up on coherence, so I will not be using SOL pool anymore
    ! N.B. I decided to comment this because there was error
    ! when an SOL update was called and I did not know how to fix it
    !*************************************
    !******* Parameters for SOL pools used for computing 
    !******* synaptic conductance caused by commom drive
    !*************************************
    paramTag = 'MUnumber_SOL-S'
    value1 = '1'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'Number_SOL'
    value1 = '1'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'Con:CMExt->SOL-S@dendrite|excitatory'
    value1 = '100'
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
    !!!!!!!!!!!!!!!! Independent noise
    ! N.B. noiseFunc is not being used anymore, it is passed as argument now
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
    allocate(synapticInput(timeLength))
    allocate(force(timeLength))
    allocate(commomDriveG(timeLength))
    allocate(MNDend_mV(timeLength))
    allocate(V(timeLength))
    allocate(Vs(poolSize, timeLength)) ! For all MNs in simulation

    t = [(dt*(i-1), i=1, timeLength)]
    
    print *, 'Running simulation'
    call cpu_time(tic)

    !*************************************
    !*************** Preparing proper input
    !*************************************
    if (inputParam.eq.'s') then 
        ! DCI+IN (noise strategy)
        !FR(:) = 25_wp + 25_wp*sin(2*pi*10*t*1e-3)!square wave was 52.5_wp + 52.5_wp*sin(2*pi*10*t*1e-3)
        ! DCI+IN (descending command strategy)
        if (inputMVC.eq.'05') then
            FR(:) = 210_wp + (20_wp + 20_wp*sin(2*pi*10*t*1e-3))
        else if (inputMVC.eq.'70') then
            FR(:) = 910_wp + (75_wp + 75_wp*sin(2*pi*10*t*1e-3))
        end if
    else if (inputParam.eq.'o') then
        ! DCI+IN (noise strategy)
        !FR(:) = 25_wp + 25_wp*sin(2*pi*10*t*1e-3)!square wave was 52.5_wp + 52.5_wp*sin(2*pi*10*t*1e-3)
        ! DCI+IN (descending command strategy)
        if (inputMVC.eq.'05') then
            FR(:) = 121_wp + (20_wp + 20_wp*sin(2*pi*10*t*1e-3))
        else if (inputMVC.eq.'70') then
            FR(:) = 645_wp + (75_wp + 75_wp*sin(2*pi*10*t*1e-3))
        end if
    else if (inputParam.eq.'d') then
        if (inputMVC.eq.'05') then
            FR(:) = 323_wp + (20_wp + 20_wp*sin(2*pi*10*t*1e-3))
        else if (inputMVC.eq.'70') then
            FR(:) = 1160_wp + (75_wp + 75_wp*sin(2*pi*10*t*1e-3))
        end if
    else if (inputParam.eq.'h') then
        if (inputMVC.eq.'05') then
            FR(:) = 170_wp + (20_wp + 20_wp*sin(2*pi*10*t*1e-3))
        else if (inputMVC.eq.'70') then
            FR(:) = 775_wp + (75_wp + 75_wp*sin(2*pi*10*t*1e-3))
        end if
    endif

    ! DCI+IN (noise strategy)
    !if (inputParam.eq.'s') then
    !    noiseFR = 1040_wp
    !else if (inputParam.eq.'o') then
    !    noiseFR = 575_wp
    !endif
    ! DCI+IN (descending command strategy), chosen arbitrarily to be small, after FR definition
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
            if (j == 1) then ! MG pool, in this case
                do k = 1, poolSize
                    Vs(k,i) = motorUnitPools(j)%v_mV(2*(k))
                end do
            else ! j=2, Sol pool
                MNDend_mV(i) = motorUnitPools(j)%v_mV(2*(1)-1)
                synapticInput(i) = motorUnitPools(j)%iIonic(1)
                commomDriveG(i) = synapticInput(i)/(MNDend_mV(i) - 70)
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

    !****************************
    !******** Measuring synchrony
    !****************************
    ! Synchrony measure computed according to Golomb, Hansel and Mato (2001)
    ! Global part
    ! Gathering only firing MNs
    allocate(firingMNsIdx(1))
    allocate(auxFiring(1))
    k = 1
    firingMNsIdx(1) = int(motorUnitPools(1)%poolSomaSpikes(1,2))
    do i = 2, size(motorUnitPools(1)%poolSomaSpikes(:,2))
        if (any(firingMNsIdx == int(motorUnitPools(1)%poolSomaSpikes(i,2)))) cycle
        ! Dynamically allocate new values
        k = k + 1
        do j = 1, size(auxFiring)
            auxFiring(j) = firingMNsIdx(j)
        end do
        deallocate(firingMNsIdx)
        allocate(firingMNsIdx(k))
        do j = 1, size(auxFiring)
            firingMNsIdx(j) = auxFiring(j)
        end do
        firingMNsIdx(k) = int(motorUnitPools(1)%poolSomaSpikes(i,2))
        deallocate(auxFiring)
        allocate(auxFiring(k))
    end do
    sampleSize = size(firingMNsIdx)

    ! Each line of Vs is the membrane potential of each MN
    allocate(Vsf(sampleSize, timeLength)) ! For all MNs in simulation
    do i = 1, sampleSize
        Vsf(i, :) = Vs(firingMNsIdx(i), :)
    end do
    V = sum(Vsf, dim=1)/sampleSize
    ! Eliminating first 100 ms = n*0.05 ms => n = 2000
    delta = sum(V(2000:)**2)/size(t(2000:)) - (sum(V(2000:))/size(t(2000:)))**2

    ! Individual part
    sumDeltai = 0
    do i=1, sampleSize
        sumDeltai = sumDeltai + sum(Vsf(i, 2000:)**2)/size(t(2000:))-&
            (sum(Vsf(i, 2000:))/size(t(2000:)))**2
    end do
    sumDeltai = sumDeltai/sampleSize

    ! Synchrony measure
    sync = delta/sumDeltai
    print *, '------ Sync coefficient ------'
    print *, sync
    print *, '------------------------------'

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

    ! Saving sync coefficient
    filename = trim(path) // trim(folderName) // "sync" // trim(inputMVC) // trim(inputParam) // ".dat"
    open(1, file=filename, status = 'replace')
    write(1, '(F15.8)') sync
    close(1)


    ! Saving input conductance and EMG
    filename = trim(path) // trim(folderName) // "g_emg" // trim(inputMVC) // trim(inputParam) // ".dat"
    open(1, file=filename, status = 'replace')
    do i = 1, timeLength
        write(1, '(F15.6, 1X, F15.6)') commomDriveG(i), motorUnitPools(1)%emg(i)
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
    if (inputParam.eq.'s') then
        call interneuronPools(1)%reset()
    endif
    do j = 1, size(motorUnitPools)
        call motorUnitPools(j)%reset()
    end do
    call neuralTractPools(1)%reset()
    do j = 1, size(synapticNoisePools)
        call synapticNoisePools(j)%reset()
    end do

    call cpu_time(toc)
    print '(F15.6, A)', toc - tic, ' seconds'
    
end program FrequencyAnalysis
