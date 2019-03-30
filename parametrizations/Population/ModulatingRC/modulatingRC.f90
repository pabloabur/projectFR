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

program modulatingRC
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
    integer :: i, j, k
    real(wp), dimension(:), allocatable :: t, force, soma_mV
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp), dimension(:), allocatable :: FR
    integer :: GammaOrder 
    character(len = 80) :: pool, group
    character(len = 100) :: filename = '../../conf.rmto'
    character(len = 100) :: path = '/home/pablo/osf/Master-Thesis-Data/population/'
    character(len = 100) :: folderName = 'modulation/'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools     
    type(AfferentPool), dimension(:), allocatable:: afferentPools     
    character(len=7) :: inputTrial
    character(len=5) :: inputParam
    character(len=1) :: inputMod
    character(len=20) :: gmaxS, gmaxFR, gmaxFF
    character(len=80) :: paramTag
    character(len=80) :: value1, value2
    ! Input parameters
    real(wp), dimension(:), allocatable :: V
    real(wp), dimension(:,:), allocatable :: Vs, Vsf
    real(wp) :: delta, sumDeltai, sync
    integer :: sampleSize ! Quantity of MNs firing
    integer, dimension(:), allocatable :: firingMNsIdx, auxFiring
    integer :: poolSize = 300 ! Quantify of used MNs
    real(wp) :: dt
    real(wp) :: tf
    logical, parameter :: probDecay = .false.
    character(len=3), parameter :: nS = '75', nFR = '75', &
        nFF = '150', nCM = '400', nMN = '300', nRC = '600' ! nS+nFR+nFF

    call init_random_seed()

    !*************************************
    !******* Getting input
    !*************************************
    call get_command_argument(3, inputTrial)
    call get_command_argument(2, inputMod)
    if (inputMod.eq.'d') then
        print *, 'double modulation'
    else if (inputMod.eq.'s') then
        print *, 'standard modulation'
    else if (inputMod.eq.'h') then
        print *, 'half modulation'
    else if (inputMod.eq.'o') then
        print *, 'open loop'
    else
        print *, 'Wrong modulation option. Available options are:'
        print *, '- d'
        print *, '- s'
        print *, '- h'
        print *, '- o'
        stop (1)
    endif
    
    call get_command_argument(1, inputParam)
    if (inputParam.eq.'high') then
        print *, 'high input'
    else if (inputParam.eq.'low') then
        print *, 'low input'
    else if (inputParam.eq.'max') then
        print *, 'max input'
    else
        print *, 'Wrong parametrization option. Available options are:'
        print *, '- max'
        print *, '- high'
        print *, '- low'
        stop (1)
    endif

    !*************************************
    !******* Changing basic parameters
    !*************************************
    conf = Configuration(filename)
    if (inputParam.eq.'high') then
        conf%simDuration_ms = 2000
    else if (inputParam.eq.'low') then
        conf%simDuration_ms = 2000
    else if (inputParam.eq.'max') then
        conf%simDuration_ms = 500
    else
        print *, 'Wrong parametrization option. Available options are:'
        print *, '- max'
        print *, '- high'
        print *, '- low'
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

    if (inputParam.eq.'high') then
        GammaOrder = 1
    else if (inputParam.eq.'max') then
        GammaOrder = 1
    else if (inputParam.eq.'low') then
        GammaOrder = 7
    else
        print *, 'Wrong parametrization option. Available options are:'
        print *, '- high'
        print *, '- low'
        stop (1)
    endif

    ! Removing influence of stimulus (required)
    paramTag = 'stimIntensity_PTN'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    ! Conductances
    if (inputMod.eq.'d') then
        gmaxS = '0.2380'
        gmaxFR = '0.2380'
        gmaxFF = '0.1880'
    else if (inputMod.eq.'s') then
        gmaxS = '0.1190'
        gmaxFR = '0.1190'
        gmaxFF = '0.0940'
    else if (inputMod.eq.'h') then
        gmaxS = '0.0595'
        gmaxFR = '0.0595'
        gmaxFF = '0.0470'
    else if (inputMod.eq.'o') then
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

    print *, 'Building neural elements'
    allocate(neuralTractPools(1))
    pool = 'CMExt'
    neuralTractPools(1) = NeuralTract(conf, pool)

    allocate(motorUnitPools(1))
    pool = 'MG'
    motorUnitPools(1) = MotorUnitPool(conf, pool)    

    if (inputMod.ne.'o') then
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
    allocate(force(timeLength))
    allocate(FR(timeLength))
    allocate(soma_mV(timeLength))
    allocate(Vs(poolSize, timeLength)) ! For all MNs in simulation
    allocate(V(timeLength))

    t = [(dt*(i-1), i=1, timeLength)]
    
    if (inputParam.eq.'low') then
        !!!!!!!!! 5% MVC
        if (inputMod.eq.'d') then
            FR(:) = 348_wp ! strong
        else if (inputMod.eq.'s') then
            FR(:) = 255_wp ! medium
        else if (inputMod.eq.'h') then
            FR(:) = 210_wp ! weak
        else if (inputMod.eq.'o') then
            FR(:) = 160_wp
        endif
    else if (inputParam.eq.'max') then
        ! uncompensated case
        !FR(:) = 1500_wp
        ! compensated case
        if (inputMod.eq.'d') then
            FR(:) = 2190_wp ! strong
        else if (inputMod.eq.'s') then
            FR(:) = 1650_wp ! medium
        else if (inputMod.eq.'h') then
            FR(:) = 1450_wp ! weak
        else if (inputMod.eq.'o') then
            FR(:) = 1245_wp
        endif
    else if (inputParam.eq.'high') then
        !!!!!!!!! 70% MVC
        if (inputMod.eq.'d') then
            FR(:) = 1187_wp ! strong
        else if (inputMod.eq.'s') then
            FR(:) = 950_wp ! medium (and uncompensated case)
        else if (inputMod.eq.'h') then
            FR(:) = 825_wp ! weak
        else if (inputMod.eq.'o') then
            FR(:) = 705_wp
        endif
    else
        print *, 'Wrong parametrization option'
        stop (1)
    endif

    call cpu_time(tic)

    !*************************************
    !*************** Running simulation
    !*************************************
    do i = 1, size(t)
        do j = 1, size(neuralTractPools)
            call neuralTractPools(j)%atualizePool(t(i), FR(i), GammaOrder)
        end do
        do j = 1, size(motorUnitPools)
            call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
            if (j==1) then
                soma_mV(i) = motorUnitPools(j)%v_mV(2*(1))
            endif
            do k = 1, poolSize
                Vs(k,i) = motorUnitPools(j)%v_mV(2*(k))
            end do
        end do
        do j = 1, size(synapticNoisePools)
            if (synapticNoisePools(j)%pool.eq.'MG') then
                call synapticNoisePools(j)%atualizePool(t(i), 0.0_wp)
            else if (synapticNoisePools(j)%pool.eq.'RC_ext') then
                call synapticNoisePools(j)%atualizePool(t(i), 7.0_wp)
            else
                print *, 'Error assigning noise value to pool'
                stop (1)
            endif
        end do
        if (inputMod.ne.'o') then
            do j = 1, size(interneuronPools)
                call interneuronPools(j)%atualizeInterneuronPool(t(i))
            end do
        endif
    end do

    call motorUnitPools(1)%listSpikes()
    if (inputMod.ne.'o') then
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
    folderName = trim(folderName) // trim(inputParam) // '/trial' // trim(inputTrial) // '/'

    filename = trim(path) // trim(folderName) // "force" // trim(inputMod) // ".dat"
    open(1, file=filename, status = 'replace')
    do i = 1, timeLength
           write(1, '(F9.3, 1X, F15.6)') t(i), force(i)
    end do
    close(1)

    filename = trim(path) // trim(folderName) // "sync" // trim(inputMod) // ".dat"
    open(1, file=filename, status = 'replace')
    write(1, '(F15.8)') sync
    close(1)

    ! Saving spikes
    filename = trim(path) // trim(folderName) // "spike" // trim(inputMod) // ".dat"
    open(1, file=filename, status = 'replace')
    do i = 1, size(motorUnitPools(1)%poolSomaSpikes, 1)
        write(1, '(F15.6, 1X, F15.1)') motorUnitPools(1)%poolSomaSpikes(i,1), &
            motorUnitPools(1)%poolSomaSpikes(i,2)
    end do
    close(1)

    !*************************************
    !*************** Plotting
    !*************************************
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

    !call gp%title('force')
    !call gp%xlabel('t (s))')
    !call gp%ylabel('force')
    !call gp%plot(t, force, 'with points pt 5 lc rgb "#0008B0"')

    !call gp%title('descending command')
    !call gp%xlabel('t (s))')
    !call gp%ylabel('fr')
    !call gp%plot(t, FR, 'with points pt 5 lc rgb "#0008B0"')

    !call gp%title('membrane')
    !call gp%xlabel('t (ms))')
    !call gp%ylabel('volts (mV)')
    !call gp%plot(t, soma_mV, 'with line lw 2 lc rgb "#0008B0"') 

    !*************************************
    !*************** Resets
    !*************************************
    !call interneuronPools(1)%reset()
    !do j = 1, size(synapticNoisePools)
    !    call synapticNoisePools(j)%reset()
    !end do
    !call neuralTractPools(1)%reset()
    !call motorUnitPools(1)%reset()
    !call synapticNoisePools(1)%reset()

    call cpu_time(toc)
    print '(F15.6, A)', toc - tic, ' seconds'

end program modulatingRC
