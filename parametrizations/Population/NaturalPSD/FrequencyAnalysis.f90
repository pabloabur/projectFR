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
    integer :: i, j, k, l
    real(wp), dimension(:), allocatable :: t, &
        synapticInput, MNDend_mV, commomDriveG, V
    real(wp), allocatable :: Vs(:,:), inFR(:,:)
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp), dimension(:), allocatable :: FR
    integer :: GammaOrder 
    character(len = 80) :: pool, group
    character(len = 100) :: filename = '../../conf.rmto'
    character(len = 100) :: path = '/home/pablo/osf/Master-Thesis-Data/population/'
    character(len = 100) :: folderName = 'psd/natural/trial3/'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools     
    type(AfferentPool), dimension(:), allocatable:: afferentPools     
    character(len=1), dimension(2) :: param
    character(len=4) :: inputParam
    character(len=80) :: paramTag
    character(len=80) :: value1, value2
    real(wp), dimension(:), allocatable :: force
    real(wp) :: dt
    real(wp) :: tf
    logical, parameter :: probDecay = .false.
    ! Noise parameters
    character(len=3), parameter :: noiseCon = '100'
    character(len=1), parameter :: noiseOrder = '1'
    character(len=2), parameter :: noiseTarget = 'FR'
    !character(len=3), parameter :: noiseFunc = '275' ! Natural PSD
    character(len=3), parameter :: noiseFunc = '625' ! Step simulation
    !character(len=2), parameter :: noiseg = '10' ! '2.1' for compensation
    character(len=4), parameter :: noiseg = '0.04' !  Step simulation
    ! Quantity parameters
    character(len=3), parameter :: nS = '75', nFR = '75', &
        nFF = '150', nCM = '400', nMN = '300', nRC = '600' ! nS+nFR+nFF
    !! Free inputs
    !! Step inputs
    integer :: sampleSize = 105 ! Quantify of used MNs
    integer :: poolSize = 300 ! Quantify of MNs
    real(wp) :: gc, gld, gls, Vth
    real(wp), dimension(:), allocatable :: Irh
    ! wgn parameters
    !real(kind=4) :: Irh(1000)
    !integer(kind=4) :: errcode
    !TYPE (VSL_STREAM_STATE) :: stream
    !real(kind=4) a,sigma
    !integer :: n

    call init_random_seed()

    ! wgn parameters
    !a = 0.0
    !sigma = 1.0
    !n = 1000
    !errcode = vsrnggaussian(VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER, stream, n, Irh, a, sigma)
    !print *, Irh

    !*************************************
    !******* Getting input
    !*************************************
    call get_command_argument(1, inputParam)
    if (inputParam.eq.'free') then
        print *, 'natural input'
    else if (inputParam.eq.'step') then
        print *, 'step input'
    else if (inputParam.eq.'injc') then
        print *, 'injected current input'
    else
        print *, 'Wrong parametrization option'
        stop (1)
    endif

    !*************************************
    !******* Changing basic parameters
    !*************************************
    conf = Configuration(filename)
    conf%simDuration_ms = 9000

    param = ['c', 'o']

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


    ! This is not being used anymore
    !!!!!!!!!!!!!!!! Independent noise
    !paramTag = 'Con:Noise>MG-@dendrite|excitatory'
    !value1 = noiseCon
    !value2 = ''
    !call conf%changeConfigurationParameter(paramTag, value1, value2)
    !paramTag = 'NoiseGammaOrder_MG'
    !value1 = noiseOrder
    !value2 = ''
    !call conf%changeConfigurationParameter(paramTag, value1, value2)
    !paramTag = 'NoiseTarget_MG'
    !value1 = noiseTarget
    !value2 = ''
    !call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'NoiseFunction_MG'
    value1 = '0'!noiseFunc
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    !paramTag = 'gmax:Noise>MG-@dendrite|excitatory'
    !value1 = noiseg
    !value2 = ''
    !call conf%changeConfigurationParameter(paramTag, value1, value2)

    GammaOrder = 1
    ! Removing influence of stimulus (required)
    paramTag = 'stimIntensity_PTN'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

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
    !paramTag = 'Con:Noise>SOL-@dendrite|excitatory'
    !value1 = noiseCon
    !value2 = ''
    !call conf%changeConfigurationParameter(paramTag, value1, value2)
    !paramTag = 'NoiseGammaOrder_SOL'
    !value1 = noiseOrder
    !value2 = ''
    !call conf%changeConfigurationParameter(paramTag, value1, value2)
    !paramTag = 'NoiseTarget_SOL'
    !value1 = noiseTarget
    !value2 = ''
    !call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'NoiseFunction_SOL'
    value1 = '0'!noiseFunc
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    !paramTag = 'gmax:Noise>SOL-@dendrite|excitatory'
    !value1 = noiseg
    !value2 = ''
    !call conf%changeConfigurationParameter(paramTag, value1, value2)

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
    allocate(synapticInput(timeLength))
    allocate(force(timeLength))
    allocate(commomDriveG(timeLength))
    allocate(MNDend_mV(timeLength))
    allocate(V(timeLength))
    allocate(Vs(sampleSize, timeLength)) ! For all MNs in simulation
    allocate(inFR(poolSize, timeLength)) ! For all MNs in simulation
    allocate(Irh(sampleSize))

    t = [(dt*(i-1), i=1, timeLength)]
    
    print *, 'Running simulation'
    do k = 1, size(param)
        call cpu_time(tic)

        !*************************************
        !*************** Preparing proper input
        !*************************************
        if (inputParam.eq.'free') then
            if (param(k).eq.'c') then
                FR(:) = 200_wp
            else if (param(k).eq.'o') then
                FR(:) = 160_wp
            endif
        else if (inputParam.eq.'step') then
            ! Loading commom synaptic input
            open(1, file='inputs/csi.txt', status = 'old')
                do i = 1, timeLength
                    read(1, '(F9.6)') FR(i) ! Only works if I use format, but requires
                    FR(i) = 1000000*FR(i)   ! this gambiarra
                end do
            close(1)
            ! Loading independent noise, depending on inhibition
            if (param(k).eq.'c') then
                open(1, file='inputs/inc.txt', status = 'old')
                    do i = 1, poolSize
                        read(1, *) inFR(i, :)
                    end do
                close(1)
            else if (param(k).eq.'o') then
                open(1, file='inputs/ino.txt', status = 'old')
                    do i = 1, poolSize
                        read(1, *) inFR(i, :)
                    end do
                close(1)
            endif
        else if (inputParam.eq.'injc') then
            if (param(k).eq.'c') then
                FR(:) = 20_wp ! This is a guess
            else if (param(k).eq.'o') then
                FR(:) = 10_wp ! This is a guess
            endif
            ! Investigating current needed to make MNs fire at a 
            ! certain rate. I just used this to find aproximately
            ! what that current would be, so the following piece of 
            ! code is just for debuging. Loop in the simulation loop
            ! should then do something similar to a linear regression 
            ! on the values chosen
            Irh(:) = 0_wp ! Initialization
            do i = 1, 1!sampleSize
                gls = motorUnitPools(1)%unit(i)%compartments(1)%gLeak_muS
                gld = motorUnitPools(1)%unit(i)%compartments(2)%gLeak_muS
                ! cytR is 70 for all MNs
                gc = 200.0_wp/(&
                    70_wp*motorUnitPools(1)%unit(i)%compartments(2)%length_mum/&
                        (pi*(motorUnitPools(1)%unit(i)%compartments(2)%diameter_mum/2)**2)+&
                    70_wp*motorUnitPools(1)%unit(i)%compartments(1)%length_mum/&
                        (pi*(motorUnitPools(1)%unit(i)%compartments(1)%diameter_mum/2)**2))
                Vth = motorUnitPools(1)%unit(i)%threshold_mV
                Irh(:) = 5.8_wp!Vth*(gls+gld*gc/(gld+gc))
            end do
        end if

        !*************************************
        !*************** Running simulation
        !*************************************
        do i = 1, size(t)        
            ! Preparations for each case
            if (inputParam.eq.'injc') then
                do j = 1, sampleSize
                    motorUnitPools(1)%iInjected(2*(j)) = 2_wp/15_wp*(j-16_wp)+7.8_wp
                end do
            else if (inputParam.eq.'step') then
                do j = 1, poolSize
                    motorUnitPools(1)%iInjected(2*(j)-1) = inFR(j,i)
                end do
                ! 200 is arbitrary. It was picked from IN that could not activate MNs
                motorUnitPools(2)%iInjected(2*(1)-1) = inFR(250,i)
            end if

            ! Updating elements
            do j = 1, size(neuralTractPools)
                call neuralTractPools(j)%atualizePool(t(i), FR(i), GammaOrder)
            end do
            do j = 1, size(synapticNoisePools)
                call synapticNoisePools(j)%atualizePool(t(i))
            end do
            if (param(k).eq.'c') then
                do j = 1, size(interneuronPools)
                    call interneuronPools(j)%atualizeInterneuronPool(t(i))
                end do
            endif
            do j = 1, size(motorUnitPools)
                if (j == 1) then ! MG pool, in this case
                    call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
                    do l = 1, sampleSize
                        Vs(l,i) = motorUnitPools(j)%v_mV(2*(l))
                    end do
                else ! j=2, Sol pool
                    call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
                    MNDend_mV(i) = motorUnitPools(j)%v_mV(2*(1)-1)
                    synapticInput(i) = motorUnitPools(j)%iIonic(1)
                    commomDriveG(i) = synapticInput(i)/(MNDend_mV(i) - 70)
                end if
            end do
        end do

        call neuralTractPools(1)%listSpikes()
        call motorUnitPools(1)%listSpikes()
        call motorUnitPools(1)%getMotorUnitPoolEMG()
        if (param(k).eq.'c') then
            call interneuronPools(1)%listSpikes()
        endif

        do i = 1, timeLength
            force(i) = motorUnitPools(1)%NoHillMuscle%force(i)
        end do

        !*************************************
        !*************** Saving data
        !*************************************
        ! Saving spikes
        if (param(k).eq.'c') then
            filename = trim(path) // trim(folderName) // "spksc.dat"
        else if (param(k).eq.'o') then
            filename = trim(path) // trim(folderName) // "spkso.dat"
        endif
        open(1, file=filename, status = 'replace')
        do i = 1, size(motorUnitPools(1)%poolSomaSpikes, 1)
            write(1, '(F15.6, 1X, F15.1)') motorUnitPools(1)%poolSomaSpikes(i,1), &
                motorUnitPools(1)%poolSomaSpikes(i,2)
        end do
        close(1)

        ! Saving force
        if (param(k).eq.'c') then
            filename = trim(path) // trim(folderName) // "forcec.dat"
        else if (param(k).eq.'o') then
            filename = trim(path) // trim(folderName) // "forceo.dat"
        endif
        open(1, file=filename, status = 'replace')
        do i = 1, timeLength
               write(1, '(F15.2, 1X, F15.6)') t(i), force(i)
        end do
        close(1)

        ! Saving input conductance and EMG
        if (param(k).eq.'c') then
            filename = trim(path) // trim(folderName) // "g_emgc.dat"
        else if (param(k).eq.'o') then
            filename = trim(path) // trim(folderName) // "g_emgo.dat"
        endif
        open(1, file=filename, status = 'replace')
        do i = 1, timeLength
            write(1, '(F15.6, 1X, F15.6)') commomDriveG(i), motorUnitPools(1)%emg(i)
        end do
        close(1)

        !*************************************
        !*************** Plotting
        !*************************************
        !if (param(k).eq.'c') then
        !    call gp%title('RC spike instants at the soma')
        !    call gp%xlabel('t (s))')
        !    call gp%ylabel('Interneuron index')
        !    call gp%plot(interneuronPools(1)%poolSomaSpikes(:,1), &
        !    interneuronPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')
        !end if

        call gp%title('descending command')
        call gp%xlabel('t (ms))')
        call gp%ylabel('descending command rate (pps)')
        call gp%plot(t, FR, 'with line lw 2 lc rgb "#0008B0"')

        call gp%title('independent noise')
        call gp%xlabel('t (ms))')
        call gp%ylabel('rate (pps)')
        call gp%plot(t, inFR(239, :), 'with line lw 2 lc rgb "#0008B0"')

        !call gp%title('MN spike instants at the soma')
        !call gp%xlabel('t (s))')
        !call gp%ylabel('Motoneuron index')
        !call gp%plot(motorUnitPools(1)%poolSomaSpikes(:,1), &
        !motorUnitPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

        !call gp%title('force')
        !call gp%xlabel('t (ms))')
        !call gp%ylabel('force (N)')
        !call gp%plot(t, motorUnitPools(1)%NoHillMuscle%force, 'with line lw 2 lc rgb "#0008B0"')

        !*************************************
        !*************** Reset elements
        !*************************************
        if (param(k).eq.'c') then
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
    end do
    
    !!!!!!!! Measuring synchrony
    !! Synchrony measure computed according to Golomb, Hansel and Mato (2001)
    !! Global part
    !V = sum(Vs, dim=2)/sampleSize
    !sigmaSquared = sum(V(2000:)**2)/size(t) - (sum(V(2000:))/size(t))**2

    !filename = "V.dat"
    !open(1, file=filename, status = 'replace')
    !do i = 1, timeLength
    !       write(1, '(F15.2)') V(i)
    !end do
    !close(1)

    !! Individual part
    !sigmaSquaredi = 0
    !do i=1, sampleSize
    !    sigmaSquaredi = sigmaSquaredi + sum(Vs(2000:,i)**2)/size(t)-&
    !        (sum(Vs(2000:,i))/size(t))**2
    !end do
    !sigmaSquaredi = sigmaSquaredi/sampleSize

    !! Chi squared
    !chiSquared = sigmaSquared/sigmaSquaredi
    !print *, '----------- Chi --------------'
    !print *, sqrt(chiSquared)
    !print *, '------------------------------'

end program FrequencyAnalysis
