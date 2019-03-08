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

program OnionSkin
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
    real(wp), dimension(:), allocatable :: t, FR
    integer, dimension(:), allocatable :: gammaOrder
    real(wp) :: tic, toc
    type(gpf) :: gp
    character(len = 80) :: pool, group
    character(len = 100) :: filename = '../../conf.rmto'
    character(len = 100) :: path = '/home/pablo/osf/Master-Thesis-Data/population/'
    character(len = 100) :: folderName = 'onion/false_decay/trial5/'
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
    character(len=3), parameter :: nS = '75', nFR = '75', &
        nFF = '150', nCM = '400', nMN = '300', nRC = '600' ! nS+nFR+nFF
    ! 340 gives 100% MVC
    ! 20% MVC/s, reaching 80% MVC
    integer, dimension(2) :: mvc = [310, 220]
    
    call init_random_seed()

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

    !!!!!!!!!!!!!!!! Independent noise
    paramTag = 'NoiseTarget_MG'
    value1 = 'FR'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'NoiseFunction_MG'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    ! Removing influence of stimulus (required)
    paramTag = 'stimIntensity_PTN'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    print *, 'Building neural elements'
    allocate(neuralTractPools(1))
    pool = 'CMExt'
    neuralTractPools(1) = NeuralTract(conf, pool)

    allocate(motorUnitPools(1))
    pool = 'MG'
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
    allocate(FR(timeLength))
    allocate(force(timeLength))
    allocate(gammaOrder(timeLength))

    t = [(dt*(i-1), i=1, timeLength)]

    ! Simply make all values compatible with poisson without changing rest of the code
    gammaOrder(:) = 1
    print *, 'Running simulation'
    call cpu_time(tic)
    ! Defining gamma order as a vector depending on the force of this simulation
    !do i=1, size(t)
    !    if (i*dt>0.and.i*dt<1550) then
    !        gammaOrder(i) = 7
    !    else if (i*dt>1550.and.i*dt<1910) then
    !        gammaOrder(i) = 5
    !    else if (i*dt>1910.and.i*dt<2850) then
    !        gammaOrder(i) = 4
    !    else if (i*dt>2850.and.i*dt<4050) then
    !        gammaOrder(i) = 3
    !    else if (i*dt>4050.and.i*dt<5050) then
    !        gammaOrder(i) = 2
    !    else if (i*dt>5050.and.i*dt<6164) then
    !        gammaOrder(i) = 3
    !    else if (i*dt>6164.and.i*dt<7088) then
    !        gammaOrder(i) = 4
    !    else if (i*dt>7088.and.i*dt<7475) then
    !        gammaOrder(i) = 5
    !    else if (i*dt>7475.and.i*dt<7893) then
    !        gammaOrder(i) = 7
    !    end if
    !end do

    do k = 1, size(param)
        ! Creating descending command vector
        do i=1, timeLength
            if (i*dt < 4000) then
                FR(i) = mvc(k)/(4000/dt)*i
            else if (i*dt < 5000) then
                FR(i) = mvc(k)
            else
                FR(i) = -mvc(k)/(4000/dt)*(i - 9000/dt)
            end if
        end do
        do i = 1, size(t)
            ! Updating elements
            do j = 1, size(neuralTractPools)
                call neuralTractPools(j)%atualizePool(t(i), FR(i), gammaOrder(1))
            end do
            do j = 1, size(motorUnitPools)
                call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
            end do
            if (param(k).eq.'c') then
                do j = 1, size(interneuronPools)
                    call interneuronPools(j)%atualizeInterneuronPool(t(i))
                end do
            endif
        end do

        call motorUnitPools(1)%listSpikes()
        if (param(k).eq.'c') then
            call interneuronPools(1)%listSpikes()
        end if

        call gp%title('MN spike instants at the soma')
        call gp%xlabel('t (s))')
        call gp%ylabel('Motoneuron index')
        call gp%plot(motorUnitPools(1)%poolSomaSpikes(:,1), &
        motorUnitPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

        call gp%title('force')
        call gp%xlabel('t (ms))')
        call gp%ylabel('force (N)')
        call gp%plot(t, motorUnitPools(1)%NoHillMuscle%force, 'with line lw 2 lc rgb "#0008B0"')

        call gp%title('descending command')
        call gp%xlabel('t (ms))')
        call gp%ylabel('descending command rate (pps)')
        call gp%plot(t, FR, 'with line lw 2 lc rgb "#0008B0"')

        if (param(k).eq.'c') then
            call gp%title('RC spike instants at the soma')
            call gp%xlabel('t (s))')
            call gp%ylabel('Interneuron index')
            call gp%plot(interneuronPools(1)%poolSomaSpikes(:,1), &
            interneuronPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')
        end if

        ! Saving spike times
        if (param(k).eq.'c') then
            filename = trim(path) // trim(folderName) // "outputc.dat"
        else if (param(k).eq.'o') then
            filename = trim(path) // trim(folderName) // "outputo.dat"
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
            write(1, '(F15.2, 1X, F15.6)') t(i), motorUnitPools(1)%NoHillMuscle%force(i)
        end do
        close(1)

        if (param(k).eq.'c') then
            call interneuronPools(1)%reset()
        endif
        call neuralTractPools(1)%reset()
        call motorUnitPools(1)%reset()
    end do
    call cpu_time(toc)
    print '(F15.6, A)', toc - tic, ' seconds'

end program OnionSkin
