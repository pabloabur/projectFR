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
    real(wp), dimension(:), allocatable :: t
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp), dimension(10) :: FR
    character(len = 80) :: pool, group
    integer, dimension(10) :: GammaOrder 
    character(len = 100) :: filename = '../../conf.rmto'
    character(len = 100) :: path = '/home/pablo/osf/Master-Thesis-Data/population/'
    character(len = 100) :: folderName = 'onion/false_decay/trial1/'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools     
    type(AfferentPool), dimension(:), allocatable:: afferentPools     
    character(len=5) :: param
    character(len=80) :: paramTag
    character(len=80) :: value1, value2
    real(wp), dimension(:), allocatable :: force
    ! Input parameters
    real(wp) :: dt
    real(wp) :: tf
    logical, parameter :: probDecay = .false.
    real(wp), parameter :: FFConducStrength = 0.033_wp, & 
        declineFactorMN = real(1, wp)/6, declineFactorRC = real(3.5, wp)/3
    character(len=3), parameter :: nS = '75', nFR = '75', &
        nFF = '150', nCM = '400', nMN = '300' ! nS+nFR+nFF
    character(len=3) :: nRC
    
    call init_random_seed()

    conf = Configuration(filename)
    conf%simDuration_ms = 2000

    call get_command_argument(1, param)
    if (param.eq.'c') then
        nRC = '600'
    else if (param.eq.'o') then
        nRC = '0'
    else
        print *, 'Wrong parametrization option'
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

        ! Columnar length
        paramTag = 'position:MG-'
        value1 = '0'
        value2 = '6'
        call conf%changeConfigurationParameter(paramTag, value1, value2)
        paramtag = 'NoiseFunction_RC_ext'
        value1 = '0'!'7' ! Spontaneuous activity change mean MN membrane potential
        value2 = ''
        call conf%changeconfigurationparameter(paramtag, value1, value2)

    !!!!!!!!!!!!!!!! Independent noise
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'NoiseTarget_MG'
    value1 = 'FR'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'NoiseFunction_MG'
    value1 = '0'
    value2 = ''

    !!!!!!!!!!!!!!!! Descending commands parameters
    GammaOrder = [7, 5, 4, 4, 4, 3, 2, 2, 1, 1]
    !FR = [150, 300, 450, 600, 750, 900, 1050, 1200, 1350, 1500]
    !FR = [200, 210, 220, 230, 240, 900, 1050, 1200, 1350, 1500]
    !FR = [20, 40, 60, 80, 100, 120, 140, 160, 180, 200]
    FR = [34, 68, 102, 136, 170, 204, 238, 272, 306, 340]

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

    if (param.eq.'c') then
        allocate(interneuronPools(1))
        pool = 'RC'
        group = 'ext'    
        interneuronPools(1) = InterneuronPool(conf, pool, group)
    else if (param.eq.'o') then
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

    t = [(dt*(i-1), i=1, timeLength)]

    print *, 'Running simulation'
    call cpu_time(tic)
    do k=1, size(FR)
        do i = 1, size(t)
            ! Updating elements
            !do j = 1, size(synapticNoisePools)
            !    call synapticNoisePools(j)%atualizePool(t(i))
            !end do
            do j = 1, size(neuralTractPools)
                call neuralTractPools(j)%atualizePool(t(i), FR(k), GammaOrder(1))
            end do
            do j = 1, size(motorUnitPools)
                call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
            end do
            if (param.eq.'c') then
                do j = 1, size(interneuronPools)
                    call interneuronPools(j)%atualizeInterneuronPool(t(i))
                !    RCv_mV(i) = interneuronPools(j)%v_mV(1)        
                end do
            endif
        end do

        call motorUnitPools(1)%listSpikes()
        if (param.eq.'c') then
            call interneuronPools(1)%listSpikes()
        end if

        call gp%title('MN spike instants at the soma')
        call gp%xlabel('t (s))')
        call gp%ylabel('Motoneuron index')
        call gp%plot(motorUnitPools(1)%poolSomaSpikes(:,1), &
        motorUnitPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

        !call gp%title('RC spike instants at the soma')
        !call gp%xlabel('t (s))')
        !call gp%ylabel('Interneuron index')
        !call gp%plot(interneuronPools(1)%poolSomaSpikes(:,1), &
        !interneuronPools(1)%poolSomaSpikes(:,2), 'with points pt 5 lc rgb "#0008B0"')

        call gp%title('force')
        call gp%xlabel('t (ms))')
        call gp%ylabel('force (N)')
        call gp%plot(t, motorUnitPools(1)%NoHillMuscle%force, 'with line lw 2 lc rgb "#0008B0"')

        if (param.eq.'c') then
            write(filename, '("output", I1, "c.dat")') k-1
            filename = trim(path) // trim(folderName) // filename
        else if (param.eq.'o') then
            write(filename, '("output", I1, "o.dat")') k-1
            filename = trim(path) // trim(folderName) // filename
        endif
        open(1, file=filename, status = 'replace')
        do i = 1, size(motorUnitPools(1)%poolSomaSpikes, 1)
            write(1, '(F15.6, 1X, F15.1)') motorUnitPools(1)%poolSomaSpikes(i,1), &
                motorUnitPools(1)%poolSomaSpikes(i,2)
        end do
        close(1)

        if (param.eq.'c') then
            call interneuronPools(1)%reset()
        endif
        call neuralTractPools(1)%reset()
        call motorUnitPools(1)%reset()
    end do

    call cpu_time(toc)
    print '(F15.6, A)', toc - tic, ' seconds'

end program OnionSkin
