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

module NeuralTractUnitClass
    ! '''
    ! Class that implements a neural tract unit. 
    ! It consists of a point process generator.
    ! '''
    use PointProcessGeneratorClass
    implicit none
    private
    integer, parameter :: wp = kind( 1.0d0 )
    public :: NeuralTractUnit

    type NeuralTractUnit
        integer :: index
        character(len=6) :: pool
        character(len = 6) :: neuronKind
        type(PointProcessGenerator) :: spikesGenerator
        

        contains
            procedure :: atualizeNeuralTractUnit
            !procedure :: transmitSpikes !TODO:
            procedure :: reset
    end type NeuralTractUnit

    interface NeuralTractUnit
        module procedure init_NeuralTractUnit
    end interface

    contains

        type(NeuralTractUnit) function init_NeuralTractUnit(pool, index)
            ! '''
            ! Constructor

            ! - Inputs:
            !     + **pool**: string with the name of the Neural tract.

            !     + **index**: integer corresponding to the neural tract unit identification.

            ! '''     
            character(len = 6), intent(in) :: pool
            integer, intent(in) :: index
            

            init_NeuralTractUnit%pool = pool
            init_NeuralTractUnit%neuronKind = ' ' 
            !  Integer corresponding to the neural tract unit identification.
            init_NeuralTractUnit%index = index  
            ! A PointProcessGenerator object, corresponding the generator of
            ! spikes of the neural tract unit.   
            init_NeuralTractUnit%spikesGenerator = PointProcessGenerator(index)  
            
            ! TODO: 
            ! # Build synapses       
            ! ## 
            ! self.SynapsesOut = []
            ! self.transmitSpikesThroughSynapses = []
            ! self.indicesOfSynapsesOnTarget = []
            ! TODO
        
        end function init_NeuralTractUnit

        subroutine atualizeNeuralTractUnit(self, t, FR, GammaOrder)
            ! '''

            ! - Inputs:
            !     + **t**: current instant, in ms.

            !     + **FR**:
            ! '''
            class(NeuralTractUnit), intent(inout) :: self
            real(wp), intent(in) :: t, FR
            integer, intent(in) :: GammaOrder

            call self%spikesGenerator%atualizeGenerator(t, FR, GammaOrder)
            ! TODO:
            ! if (self%terminalSpikeTrain and -1e-3 < (t - self.terminalSpikeTrain[-1][0]) < 1e-3) then
            !     self.transmitSpikes(t)
            ! end if
            ! TODO
        end subroutine

        ! TODO:        
        ! def transmitSpikes(self, t):
        !     '''
        !     - Inputs:
        !         + **t**: current instant, in ms.
        !     '''
        !     for i in xrange(len(self.indicesOfSynapsesOnTarget)):
        !         self.transmitSpikesThroughSynapses[i].receiveSpike(t, self.indicesOfSynapsesOnTarget[i])
      
        subroutine reset(self)
            class(NeuralTractUnit), intent(inout) :: self

            call self%spikesGenerator%reset()        
        end subroutine

end module NeuralTractUnitClass





      