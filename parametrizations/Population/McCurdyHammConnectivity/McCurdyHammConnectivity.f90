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

program McCurdyHammConnectivity
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
    real(wp), dimension(:,:), allocatable :: RC_mV, MN_mV
    real(wp), dimension(:), allocatable :: t, positions
    real(wp) :: tic, toc
    type(gpf) :: gp
    integer :: GammaOrder 
    character(len = 80) :: pool, group
    character(len = 100) :: filename = '../../conf.rmto'
    character(len = 100) :: path = '/home/pablo/osf/Master-Thesis-Data/population/'
    character(len = 100) :: folderName = 'McHamm/false_decay/trial4/'
    type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools     
    type(AfferentPool), dimension(:), allocatable:: afferentPools     
    character(len=5) :: param
    character(len=80) :: paramTag
    character(len=80) :: value1, value2
    ! Input parameters
    real(wp) :: dt
    real(wp) :: tf
    logical, parameter :: probDecay = .false.
    character(len=3), parameter :: nS = '75', nFR = '75', &
        nFF = '150', nRC = '600', nCM = '0', nMN = '300' ! nS+nFR+nFF
    integer, parameter :: MNNumber = 300
    integer, parameter :: RCNumber = 600
    integer, parameter :: MNs = 14 ! This variable to represent
                                   ! MN selected for stimulation
    integer :: availableMNs(MNNumber), stimulatedMNs(MNs)
    real(wp) :: rgn
    
    ! Determining MNs to be stimulated
    i = 1
    do while (i<15)
        ! Generate random number from 1 to 300
        call random_number(rgn)
        stimulatedMNs(i) = int(1 + (300 - 1)*rgn)
        if (any(availableMNs == availableMNs(i))) then
            continue
        end if
        i = i + 1
    end do

    call init_random_seed()

    conf = Configuration(filename)
    conf%simDuration_ms = 50

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

        ! ! RC spontaneous activity 
        paramtag = 'NoiseFunction_RC_ext'
        value1 = '0' ! Spontaneuous activity change mean MN membrane potential
        value2 = ''
        call conf%changeconfigurationparameter(paramtag, value1, value2)

    ! Dynamics of MN-RC synapse
    paramtag = 'dyn:MG-S>RC_ext-@soma|excitatory'
    value1 = 'None'
    value2 = ''
    call conf%changeconfigurationparameter(paramtag, value1, value2)
    paramTag = 'dyn:MG-FR>RC_ext-@soma|excitatory'
    value1 = 'None'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'dyn:MG-FF>RC_ext-@soma|excitatory'
    value1 = 'None'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    !!!!!!!!!!!!!!!! Independent noise
    paramTag = 'NoiseTarget_MG'
    value1 = 'FR'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)
    paramTag = 'NoiseFunction_MG'
    value1 = '0' ! Turned off
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    ! Removing influence of stimulus (required)
    paramTag = 'stimIntensity_PTN'
    value1 = '0'
    value2 = ''
    call conf%changeConfigurationParameter(paramTag, value1, value2)

    print *, 'Building neural elements'
    allocate(neuralTractPools(0))

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
    allocate(MN_mV(timeLength, MNNumber))
    allocate(RC_mV(timeLength, RCNumber))
    allocate(positions(MNNumber))

    t = [(dt*(i-1), i=1, timeLength)]
        
    print *, 'Running simulation'
    call cpu_time(tic)
    do k=1, size(stimulatedMNs)
        do i = 1, size(t)
            ! Stimulus to soma of stimulatedMNs(k)
            if (t(i)>10.and.t(i)<11) then
                motorUnitPools(1)%iInjected(2*(stimulatedMNs(k))) = 50
            else
                motorUnitPools(1)%iInjected(2*(stimulatedMNs(k))) = 0
            end if

            ! Updating elements
            do j = 1, size(motorUnitPools)
                call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 32.0_wp, 32.0_wp)
            end do
            do j = 1, size(interneuronPools)
                call interneuronPools(j)%atualizeInterneuronPool(t(i))
            end do

            ! Getting values
            do j = 1, MNNumber
                MN_mV(i, j) = motorUnitPools(1)%v_mV(2*(j))
            end do
            do j = 1, RCNumber
                RC_mV(i, j) = interneuronPools(1)%v_mV(j)
            end do
        end do
        if (k==1) then
            do j = 1, MNNumber
                positions(j) = motorUnitPools(1)%unit(j)%position_mm
            end do
        end if

        call motorUnitPools(1)%listSpikes()
        call interneuronPools(1)%listSpikes()

        write(filename, '("MNV", I3, ".dat")') stimulatedMNs(k)
        filename = trim(path) // trim(folderName) // filename
        open(1, file=filename, status = 'replace')
        do i = 1, timeLength
            do j = 1, MNNumber
                write(1, '(F15.6)', advance='no') MN_mV(i, j)
            end do
            write(1, *) '' ! Line break
        end do
        close(1)

        filename = trim(path) // trim(folderName) // "positions.dat"
        open(1, file=filename, status = 'replace')
        do i = 1, MNNumber
               write(1, '(F15.2)') positions(i)
        end do
        close(1)

        call motorUnitPools(1)%reset()
        call interneuronPools(1)%reset()
    end do
    
    call cpu_time(toc)
    print '(F15.6, A)', toc - tic, ' seconds'

end program McCurdyHammConnectivity
