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
program Farina
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
    real(wp), dimension(:), allocatable :: dendPotMG, t
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp), dimension(:), allocatable :: FR
    integer :: GammaOrder 
    character(len = 80) :: pool, group
    character(len = 100) :: filename = '../../conf.rmto'
    character(len = 100) :: path = '/home/pablo/git/CSTAnalysis/Data/'
    character(len = 100) :: folderName
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools     
    type(AfferentPool), dimension(:), allocatable:: afferentPools     
    character(len=2) :: inputTrial
    character(len=4) :: inputParam
    character(len=80) :: paramTag
    character(len=80) :: value1, value2
    real(wp), dimension(:), allocatable :: force
    real(wp), dimension(:,:), allocatable :: Vs
    real(wp) :: dt
    real(wp) :: tf
    real(wp) :: phase
    logical, parameter :: probDecay = .false.
    integer :: poolSize = 300 ! Quantify of used MNs
    ! Noise parameters
    real(wp) :: noiseFR
    character(len=1), parameter :: noiseOrder = '1'
    character(len=2), parameter :: noiseg = '12' !  Step simulation
    ! Quantity parameters

    call init_random_seed()

    !*************************************
    !******* Getting input
    !*************************************
    call get_command_argument(2, inputTrial)
    call get_command_argument(1, inputParam)
    if (inputParam.eq.'o') then
        print *, 'no renshaw cell'
    else if (inputParam.eq.'s') then
        print *, 'renshaw cell with standard conductance'
    else
        print *, 'Wrong parametrization option. Available options are:'
        print *, '- o'
        print *, '- s'
        stop (1)
    endif

    !*************************************
    !******* Changing basic parameters
    !*************************************
    conf = Configuration(filename)
    conf%simDuration_ms = 10000

    !!!!!!!!!!!!!!!! Independent noise
    paramTag = 'NoiseGammaOrder_MG'
    value1 = noiseOrder
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'gmax:Noise>MG-@dendrite|excitatory'
    value1 = noiseg
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    !*************************************
    !******* Open loop case
    !*************************************
    ! Opening the loop
    paramTag = 'Con:MG-S>RC_ext-@soma|excitatory'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'Con:MG-FR>RC_ext-@soma|excitatory'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'Con:MG-FF>RC_ext-@soma|excitatory'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    print *, '*************************************'
    print *, '************ Building neural elements'
    print *, '*************************************'
    allocate(neuralTractPools(1))
    pool = 'CMExt'
    neuralTractPools(1) = NeuralTract(conf, pool)

    allocate(motorUnitPools(1))
    pool = 'MG'
    motorUnitPools(1) = MotorUnitPool(conf, pool)    

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
    allocate(dendPotMG(timeLength))
    allocate(force(timeLength))
    allocate(Vs(poolSize, timeLength)) ! For all MNs in simulation

    t = [(dt*(i-1), i=1, timeLength)]
    
    print *, 'Running simulation'
    call cpu_time(tic)

    !*************************************
    !*************** Preparing proper input
    !*************************************
    noiseFR = 0_wp
    GammaOrder = 7
    if (inputParam.ne.'o') then
        FR(:) = 438_wp
    else
        FR(:) = 300_wp
    endif
    do i = 1, 40
        call random_number(phase)
        FR(:) = FR(:) + 2_wp*sin(2*pi*0.5_wp*i*t*1e-3 + phase*360_wp)
    end do

    !*************************************
    !*************** Running simulation
    !*************************************
    do i = 1, size(t)        
        ! Inject current to RC
        if (inputParam.ne.'o') then
            interneuronPools(1)%iInjected(:) = 3.5_wp
        endif

        ! Updating elements
        do j = 1, size(neuralTractPools)
            call neuralTractPools(j)%atualizePool(t(i), FR(i), GammaOrder)
        end do

        do j = 1, size(synapticNoisePools)
            if (synapticNoisePools(j)%pool.eq.'MG') then
                call synapticNoisePools(j)%atualizePool(t(i), noiseFR)
            else if (synapticNoisePools(j)%pool.eq.'RC_ext') then
                call synapticNoisePools(j)%atualizePool(t(i), 0.0_wp)
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
            do k = 1, poolSize
                Vs(k,i) = motorUnitPools(j)%v_mV(2*(k))
            end do
            dendPotMG(i) = motorUnitPools(j)%v_mV(2*(40)-1)
        end do
    end do

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
    folderName = 'trial' // trim(inputTrial) // '/'

    filename = trim(path) // trim(folderName) // "spike" // trim(inputParam) // ".dat"
    open(1, file=filename, status = 'replace')
    do i = 1, size(motorUnitPools(1)%poolSomaSpikes, 1)
        write(1, '(F15.6, 1X, F15.1)') motorUnitPools(1)%poolSomaSpikes(i,1), &
            motorUnitPools(1)%poolSomaSpikes(i,2)
    end do
    close(1)
    if (inputParam.ne.'o') then
        filename = trim(path) // trim(folderName) // "inspike" // trim(inputParam) // ".dat"
        open(1, file=filename, status = 'replace')
        do i = 1, size(interneuronPools(1)%poolSomaSpikes, 1)
            write(1, '(F15.6, 1X, F15.1)') interneuronPools(1)%poolSomaSpikes(i,1), &
                interneuronPools(1)%poolSomaSpikes(i,2)
        end do
        close(1)
    endif

    !! Saving force
    filename = trim(path) // trim(folderName) // "inout" // trim(inputParam) // ".dat"
    open(1, file=filename, status = 'replace')
    do i = 1, timeLength
           write(1, '(F15.2, 1X, F15.6, 1X, F15.6, 1X, F15.6)') t(i), &
               dendPotMG(i), force(i), FR(i)
    end do
    close(1)

    !! Saving input conductance and EMG
    !filename = trim(path) // trim(folderName) // "g_emg" // trim(inputParam) // ".dat"
    !open(1, file=filename, status = 'replace')
    !do i = 1, timeLength
    !    write(1, '(F15.6)') motorUnitPools(1)%emg(i)
    !end do
    !close(1)

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
    !call gp%xlabel('t (ms)')
    !call gp%ylabel('Motoneuron index')
    !call gp%plot(motorUnitPools(1)%poolSomaSpikes(:,1), &
    !motorUnitPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

    !call gp%title('Soma voltage')
    !call gp%xlabel('t (ms))')
    !call gp%ylabel('Volts (mV)')
    !call gp%plot(t, Vs(40,:), 'with line lw 2 lc rgb "#0008B0"')

    !call gp%title('Dendrite voltage')
    !call gp%xlabel('t (ms))')
    !call gp%ylabel('Volts (mV)')
    !call gp%plot(t, dendPotMG, 'with line lw 2 lc rgb "#0008B0"')

    !call gp%title('force')
    !call gp%xlabel('t (ms))')
    !call gp%ylabel('force (N)')
    !call gp%plot(t, motorUnitPools(1)%NoHillMuscle%force, 'with line lw 2 lc rgb "#0008B0"')

    call cpu_time(toc)
    print '(F15.6, A)', toc - tic, ' seconds'
    
end program Farina
