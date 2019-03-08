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

program RecruitmentOrder
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
    real(wp), dimension(:), allocatable :: t
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp), dimension(:), allocatable :: FR
    character(len = 80) :: pool, group
    integer, dimension(6) :: GammaOrder 
    character(len = 100) :: filename = '../../conf.rmto'
    character(len = 100) :: path = '/home/pablo/osf/Master-Thesis-Data/population/'
    character(len = 100) :: folderName = 'recruitment/false_decay/trial2/'
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
    
    call init_random_seed()

    conf = Configuration(filename)
    conf%simDuration_ms = 1000

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

    !!!!!!!!!!!!!!!! Descending commands parameters
    GammaOrder = [7, 5, 4, 3, 2, 1]

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

    t = [(dt*(i-1), i=1, timeLength)]

    do i = 1, size(t)
        FR(i) = 1500_wp/1000.0_wp*t(i)
    end do
        
    print *, 'Running simulation'
    call cpu_time(tic)
    do k = 1, size(param)
        do i = 1, size(t)
            !! Injected current in the soma of MNs
            !if (i>1) then ! If there is a value in the first
            !              ! iteration, program throws error
            !    do j = 1, size(motorUnitPools(1)%unit)
            !        motorUnitPools(1)%iInjected(2*(j)) = 40.0_wp/1000.0_wp*t(i)
            !    end do
            !end if

            ! Updating elements
            !do j = 1, size(synapticNoisePools)
            !    call synapticNoisePools(j)%atualizePool(t(i))
            !end do
            do j = 1, size(neuralTractPools)
                call neuralTractPools(j)%atualizePool(t(i), FR(i), GammaOrder(6))
            end do
            do j = 1, size(motorUnitPools)
                call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
            end do
            if (param(k).eq.'c') then
                do j = 1, size(interneuronPools)
                    call interneuronPools(j)%atualizeInterneuronPool(t(i))
                !    RCv_mV(i) = interneuronPools(j)%v_mV(1)        
                end do
            endif
        end do

        call motorUnitPools(1)%listSpikes()
        if (param(k).eq.'c') then
            call interneuronPools(1)%listSpikes()
        end if

        do i = 1, size(t)
            force(i) = motorUnitPools(1)%NoHillMuscle%force(i)
        end do

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

        if (param(k).eq.'c') then
            filename = trim(path) // trim(folderName) // "forcec.dat"
        else if (param(k).eq.'o') then
            filename = trim(path) // trim(folderName) // "forceo.dat"
        end if
        open(1, file=filename, status = 'replace')
        do i = 1, size(t)
            write(1, '(F15.6, 1X, F15.6)') t(i), force(i)
        end do
        close(1)

        if (param(k).eq.'c') then
            filename = trim(path) // trim(folderName) // "MNc.dat"
        else if (param(k).eq.'o') then
            filename = trim(path) // trim(folderName) // "MNo.dat"
        endif
        open(1, file=filename, status = 'replace')
        do i = 1, size(motorUnitPools(1)%poolSomaSpikes, 1)
            write(1, '(F15.6, 1X, F15.1)') motorUnitPools(1)%poolSomaSpikes(i,1), &
                motorUnitPools(1)%poolSomaSpikes(i,2)
        end do
        close(1)
        
        if (param(k).eq.'c') then
            filename = trim(path) // trim(folderName) // "INc.dat"
            open(1, file=filename, status = 'replace')
            do i = 1, size(interneuronPools(1)%poolSomaSpikes, 1)
                write(1, '(F15.6, 1X, F15.1)') interneuronPools(1)%poolSomaSpikes(i,1), &
                    interneuronPools(1)%poolSomaSpikes(i,2)
            end do
            close(1)
        endif

        if (param(k).eq.'c') then
            call interneuronPools(1)%reset()
        endif
        call neuralTractPools(1)%reset()
        call motorUnitPools(1)%reset()
    end do

    call cpu_time(toc)
    print '(F15.6, A)', toc - tic, ' seconds'

end program RecruitmentOrder
