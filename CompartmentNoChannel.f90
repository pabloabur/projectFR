! # '''
! #     Neuromuscular simulator in Fortran.
! #     Copyright (C) 2018  Renato Naville Watanabe

! #     This program is free software: you can redistribute it and/or modify
! #     it under the terms of the GNU General Public License as published by
! #     the Free Software Foundation, either version 3 of the License, or
! #     any later version.

! #     This program is distributed in the hope that it will be useful,
! #     but WITHOUT ANY WARRANTY; without even the implied warranty of
! #     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! #     GNU General Public License for more details.

! #     You should have received a copy of the GNU General Public License
! #     along with this program.  If not, see <http://www.gnu.org/licenses/>.

! #     Contact: renato.watanabe@ufabc.edu.br
! # '''

module CompartmentNoChannelClass
    ! '''
    ! Class that implements a neural compartment. For now it is implemented
    ! *dendrite* and *soma*.
    ! '''
    use ConfigurationClass
    use ChannelConductanceClass    
    use SynapseClass
    implicit none
    private
    integer, parameter :: wp = kind(1.0d0)
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    public :: CompartmentNoChannel

    type CompartmentNoChannel
        character(len = 9) :: compKind
        type(ChannelConductance), dimension(:), allocatable :: Channels ! ## List of ChannelConductance objects in the Compartment.
        character(len = 2)::  neuronKind
        integer :: index, numberChannels
        type(Configuration), pointer :: conf
        real(wp) :: length_mum, diameter_mum, capacitance_nF
        real(wp) :: EqPot_mV, IPump_nA, gLeak_muS
        type(Synapse), dimension(2) :: SynapsesIn
        
        contains

            procedure :: reset
            procedure :: computeCurrent

    end type CompartmentNoChannel

    interface CompartmentNoChannel
        module procedure init_CompartmentNoChannel
    end interface CompartmentNoChannel

    contains

    type(CompartmentNoChannel) function init_CompartmentNoChannel(compKind, conf, pool, index, neuronKind)
        ! '''
        ! Constructor

        ! - Inputs:
        !     + **kind**: The kind of compartment. For now, it can be *soma*, *dendrite*, 
        !     *node* or *internode*.

        !     + **conf**: Configuration object with the simulation parameters.

        !     + **pool**: string with Motor unit pool to which the motor unit belongs.

        !     + **index**: integer corresponding to the motor unit order in the pool, according to 
        !     the Henneman's principle (size principle).

        !     + **neuronKind**: string with the type of the motor unit. It can be *S* (slow), *FR* (fast and resistant), 
        !     and *FF* (fast and fatigable).
        ! '''
        character(*), intent(in) :: compKind
        class(Configuration), intent(in), target :: conf
        character(*), intent(in) :: pool
        integer, intent(in) :: index
        character(*), intent(in) ::  neuronKind
        character(len=80) :: paramTag, paramChar
        real(wp) :: area_cm2, specifRes_Ohmcm2, membraneCapacitance
        character(len = 6) :: channelKind
        character(len = 10) :: synapseKind
        
        
        init_CompartmentNoChannel%conf => conf
        ! ## String with the type of the motor unit. It can be *S* (slow), *FR* (fast and resistant), 
        ! ## and *FF* (fast and fatigable).
        init_CompartmentNoChannel%neuronKind = neuronKind
        ! ## List of summed synapses (see Lytton, 1996) that the Compartment do with other neural components.
        
        
        
        ! ## List of summed synapses (see Lytton, 1996) that the Compartment receive from other neural components.
        synapseKind = 'excitatory'
        init_CompartmentNoChannel%SynapsesIn(1) = Synapse(conf, pool, index, compKind, synapseKind, neuronKind)
        synapseKind = 'inhibitory'
        init_CompartmentNoChannel%SynapsesIn(2) = Synapse(conf, pool, index, compKind, synapseKind, neuronKind)
        
        ! ## The kind of compartment. For now, it can be *soma* or *dendrite*, *node* or *internode*..
        init_CompartmentNoChannel%compKind = compKind
        
        ! ## Integer corresponding to the motor unit order in the pool, according to 
        ! ## the Henneman's principle (size principle).
        init_CompartmentNoChannel%index = index
        
        ! ## Length of the compartment, in \f$\mu\f$m.
        paramTag = 'l@' // trim(init_CompartmentNoChannel%compKind)
        paramChar =  init_CompartmentNoChannel%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_CompartmentNoChannel%length_mum
        ! ## Diameter of the compartment, in \f$\mu\f$m.
        paramTag = 'd@' // trim(init_CompartmentNoChannel%compKind)
        paramChar =  init_CompartmentNoChannel%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_CompartmentNoChannel%diameter_mum
        area_cm2 = init_CompartmentNoChannel%length_mum * pi * init_CompartmentNoChannel%diameter_mum * 1e-8
        paramTag = 'res@' // trim(init_CompartmentNoChannel%compKind)
        paramChar = init_CompartmentNoChannel%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)specifRes_Ohmcm2
        ! ## Capacitance of the compartment, in nF.
        paramTag = 'membCapac'
        paramChar = init_CompartmentNoChannel%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)membraneCapacitance
        init_CompartmentNoChannel%capacitance_nF = membraneCapacitance * area_cm2 * 1e3
        
        ! ## Equilibrium potential, in mV.
        paramTag = 'EqPot@' // trim(init_CompartmentNoChannel%compKind)
        paramChar = init_CompartmentNoChannel%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_CompartmentNoChannel%EqPot_mV

        ! ## Pump current in the compartment, in nA.
        paramTag = 'IPump@' // trim(init_CompartmentNoChannel%compKind)
        paramChar = init_CompartmentNoChannel%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_CompartmentNoChannel%IPump_nA

        ! ## Leak conductance of the compartment, in \f$\mu\f$S.
        init_CompartmentNoChannel%gLeak_muS = (1e6 * area_cm2) / specifRes_Ohmcm2

        init_CompartmentNoChannel%numberChannels = 0    
        
        

    end function

    real(wp) function computeCurrent(self, t, V_mV) result(current)
        ! '''
        ! Computes the active currents of the compartment. Active currents are the currents from the ionic channels
        ! and from the synapses.

        ! - Inputs:
        !     + **t**: current instant, in ms.

        !     + **V_mV**: membrane potential, in mV.
        ! '''
        class(CompartmentNoChannel), intent(inout):: self
        real(wp), intent(in) :: t, V_mV
        integer :: i
        
        
        current = 0.0
        
        
        if (self%SynapsesIn(1)%numberOfIncomingSynapses > 0) current = current + self%SynapsesIn(1)%computeCurrent(t, V_mV)
        
        if (self%SynapsesIn(2)%numberOfIncomingSynapses > 0) current = current + self%SynapsesIn(2)%computeCurrent(t, V_mV)
        
        
    end function

    subroutine reset(self)
        ! '''

        ! '''
        class(CompartmentNoChannel), intent(inout):: self
        integer :: i

        do i = 1, size(self%SynapsesIn)
            call self%SynapsesIn(i)%reset()
        end do
        do i = 1, self%numberChannels
            call self%Channels(i)%reset()
        end do
    end subroutine


end module CompartmentNoChannelClass



        

        
    
    