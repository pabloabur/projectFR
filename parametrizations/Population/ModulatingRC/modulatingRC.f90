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
    integer :: i, j
    real(wp), dimension(:), allocatable :: t, force, soma_mV
    real(wp) :: tic, toc
    type(gpf) :: gp
    real(wp), dimension(:), allocatable :: FR
    integer :: GammaOrder 
    character(len = 80) :: pool, group
    character(len = 100) :: filename = '../../conf.rmto'
    character(len = 100) :: path = '/home/pablo/osf/Master-Thesis-Data/population/'
    character(len = 100) :: folderName = 'modulation/trial1/'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools     
    type(AfferentPool), dimension(:), allocatable:: afferentPools     
    character(len=7) :: params
    character(len=5) :: inputParam
    character(len=1) :: inputMod
    character(len=6) :: gmaxS, gmaxFR, gmaxFF
    character(len=80) :: paramTag
    character(len=80) :: value1, value2
    ! Input parameters
    real(wp) :: dt
    real(wp) :: tf
    logical, parameter :: probDecay = .false.
    character(len=3), parameter :: nS = '75', nFR = '75', &
        nFF = '150', nCM = '400', nMN = '300', nRC = '600' ! nS+nFR+nFF

    call init_random_seed()

    !*************************************
    !******* Getting input
    !*************************************
    call get_command_argument(2, inputMod)
    if (inputMod.eq.'d') then
        print *, 'double modulation'
    else if (inputMod.eq.'s') then
        print *, 'standard modulation'
    else if (inputMod.eq.'h') then
        print *, 'half modulation'
    else
        print *, 'Wrong modulation option. Available options are:'
        print *, '- d'
        print *, '- s'
        print *, '- h'
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
    else
        print *, 'Wrong modulation option. Available options are:'
        print *, '- high'
        print *, '- low'
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
    allocate(force(timeLength))
    allocate(FR(timeLength))
    allocate(soma_mV(timeLength))

    t = [(dt*(i-1), i=1, timeLength)]
    
    if (inputParam.eq.'low') then
        !!!!!!!!! 5% MVC
        if (inputMod.eq.'d') then
            FR(:) = 320_wp ! strong
        else if (inputMod.eq.'s') then
            FR(:) = 250_wp ! medium (and uncompensated case)
        else if (inputMod.eq.'h') then
            FR(:) = 210_wp ! weak
        endif
    else if (inputParam.eq.'max') then
        ! uncompensated case
        !FR(:) = 1500_wp
        ! compensated case
        if (inputMod.eq.'d') then
            FR(:) = 2150_wp ! strong
        else if (inputMod.eq.'s') then
            FR(:) = 1650_wp ! medium
        else if (inputMod.eq.'h') then
            FR(:) = 1450_wp ! weak
        endif
    else if (inputParam.eq.'high') then
        !!!!!!!!! 70% MVC
        if (inputMod.eq.'d') then
            FR(:) = 1150_wp ! strong
        else if (inputMod.eq.'s') then
            FR(:) = 950_wp ! medium (and uncompensated case)
        else if (inputMod.eq.'h') then
            FR(:) = 825_wp ! weak
        endif
    else
        print *, 'Wrong parametrization option. Available options are:'
        print *, '- max'
        print *, '- high'
        print *, '- low'
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
        do j = 1, size(interneuronPools)
            call interneuronPools(j)%atualizeInterneuronPool(t(i))
        end do
    end do

    call motorUnitPools(1)%listSpikes()
    call interneuronPools(1)%listSpikes()

    do i = 1, timeLength
        force(i) = motorUnitPools(1)%NoHillMuscle%force(i)
    end do

    !*************************************
    !*************** Saving data
    !*************************************
    if (inputParam.eq.'high') then
        if (inputMod.eq.'d') then
            params = "hid.dat"
        else if (inputMod.eq.'s') then
            params = "his.dat"
        else if (inputMod.eq.'h') then
            params = "hih.dat"
        endif
    else if (inputParam.eq.'max') then
        if (inputMod.eq.'d') then
            params = "mad.dat"
        else if (inputMod.eq.'s') then
            params = "mas.dat"
        else if (inputMod.eq.'h') then
            params = "mah.dat"
        endif
    else if (inputParam.eq.'low') then
        if (inputMod.eq.'d') then
            params = "lod.dat"
        else if (inputMod.eq.'s') then
            params = "los.dat"
        else if (inputMod.eq.'h') then
            params = "loh.dat"
        endif
    endif
    filename = trim(path) // trim(folderName) // "force" // params
    open(1, file=filename, status = 'replace')
    do i = 1, timeLength
           write(1, '(F9.3, 1X, F15.6)') t(i), force(i)
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
