! '''
!     Neuromuscular simulator in Python.
!     Copyright (C) 2018  Renato Naville Watanabe

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

module SynapsesFactoryModule
    ! '''
    ! Module to build all the synapses in the system.
    ! '''
    use ConfigurationClass
    use MotorUnitPoolClass
    use NeuralTractClass
    use DynamicalArrays
    use SynapsePointerClass
    implicit none
    integer, parameter :: wp = kind( 1.0d0 )
    private :: wp
    contains

        subroutine synapseFactory(conf, neuralTractPools, motorUnitPools)
            ! '''
            ! Constructor

            ! - Inputs:
            !     + **conf**: Configuration object with the simulation parameters.

            !     + **pools**: list of all the pools in the system.

            
            ! '''
            class(Configuration), intent(inout) :: conf
            class(NeuralTract), intent(inout) :: neuralTractPools(:)
            class(MotorUnitPool), intent(inout), target:: motorUnitPools(:)
            integer :: numberOfSynapses
            integer :: poolOut, unitOut, synapseIn, poolIn, unitIn, compartmentIn, synapseComp
            character(len = 80) :: neuralSource, paramTag, paramChar
            real(wp) :: tau, var, gmax, conn, delay, declineFactor
            character(len=80) :: dyn
            integer :: pos, i, newIndex, j
            real(wp) :: randomNumber
            real(wp) :: neuronsDistance, weight
            type(SynapsePointer), dimension(:), allocatable:: tempTransmitSpikes

            ! ## Total number of synapses in the system.
            numberOfSynapses = 0
            
            !NeuralTract to MotorUnitPool
            do poolOut = 1, size(neuralTractPools)
                do unitOut = 1, size(neuralTractPools(poolOut)%unit)
                    neuralSource = trim(neuralTractPools(poolOut)%pool) // '-' // neuralTractPools(poolOut)%unit(unitOut)%neuronKind
                    neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut = conf%determineSynapses(neuralSource)
                    ! #print pools[poolOut].pool
                    ! #print pools[poolOut].unit[unitOut].SynapsesOut
                    do synapseIn = 1, size(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item)
                        paramTag = 'Con:' // trim(neuralTractPools(poolOut)%pool) // '-' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string) &
                                    // '-' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string) &
                                    // '@' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string) &
                                    // '|' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                        paramChar = conf%parameterSet(paramTag,neuralTractPools(poolOut)%pool, 0)
                        read(paramChar, *)conn
                        conn = conn / 100.0

                        paramTag = 'gmax:' // trim(neuralTractPools(poolOut)%pool) // '-' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string) &
                                    // '-' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string) &
                                    // '@' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string) &
                                    // '|' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                        paramChar = conf%parameterSet(paramTag,neuralTractPools(poolOut)%pool, 0)
                        read(paramChar,*)gmax

                        paramTag = 'delay:' // trim(neuralTractPools(poolOut)%pool) // '-' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string) &
                                    // '-' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string) &
                                    // '@' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string) &
                                    // '|' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                        
                        paramChar = conf%parameterSet(paramTag, neuralTractPools(poolOut)%pool, 0)
                        read(paramChar, *)delay
                        
                        paramTag = 'dec:' // trim(neuralTractPools(poolOut)%pool) // '-' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string) &
                                    // '-' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string) &
                                    // '@' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string) &
                                    // '|' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                        paramChar = conf%parameterSet(paramTag, neuralTractPools(poolOut)%pool, 0)
                        if (trim(paramChar)=='inf') then
                            declineFactor = 1e6
                        else 
                            read(paramChar, *)declineFactor
                        end if
                        
                        paramTag = 'dyn:' // trim(neuralTractPools(poolOut)%pool) // '-' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string) &
                                    // '-' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string) &
                                    // '@' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string) &
                                    // '|' &
                                    // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                        dyn = conf%parameterSet(paramTag, neuralTractPools(poolOut)%pool, 0)
                        
                        if (trim(dyn).ne.'None') then
                            paramTag = 'var:' // trim(neuralTractPools(poolOut)%pool) // '-' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            
                            paramChar = conf%parameterSet(paramTag,neuralTractPools(poolOut)%pool, 0)
                            read(paramChar, *)var

                            paramTag = 'tau:' // trim(neuralTractPools(poolOut)%pool) // '-' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            paramChar = conf%parameterSet(paramTag,neuralTractPools(poolOut)%pool, 0)
                            read(paramChar, *)tau
                        else
                            var = 0.0
                            tau = 1e6
                        end if
                        
                        do poolIn = 1, size(motorUnitPools)
                            if (trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)==&
                                trim(motorUnitPools(poolIn)%pool)) then
                                do unitIn = 1, size(motorUnitPools(poolIn)%unit)
                                    do compartmentIn = 1, size(motorUnitPools(poolIn)%unit(unitIn)%Compartments)
                                        if ((trim(neuralTractPools(poolOut)%unit(unitOut)%&
                                                  SynapsesOut%item(synapseIn)%item(1)%string)==&
                                             trim(motorUnitPools(poolIn)%pool)) .and.&
                                            (trim(neuralTractPools(poolOut)%unit(unitOut)%&
                                                  SynapsesOut%item(synapseIn)%item(2)%string)==&
                                             trim(motorUnitPools(poolIn)%unit(unitIn)%neuronKind)) .and.&
                                            (trim(neuralTractPools(poolOut)%unit(unitOut)%&
                                                  SynapsesOut%item(synapseIn)%item(3)%string)==&
                                             trim(motorUnitPools(poolIn)%unit(unitIn)%Compartments(compartmentIn)%compKind))) then
                                            call random_number(randomNumber)
                                            if (randomNumber <= conn) then
                                                do synapseComp = 1, size(motorUnitPools(poolIn)%unit(unitIn)%&
                                                                        Compartments(compartmentIn)%SynapsesIn)
                                                    if (motorUnitPools(poolIn)%unit(unitIn)%Compartments(compartmentIn)%&
                                                        SynapsesIn(synapseComp)%synapseKind ==&
                                                        neuralTractPools(poolOut)%unit(unitOut)%&
                                                        SynapsesOut%item(synapseIn)%item(4)%string) then
                                                        
                                                        if (declineFactor<1e5) then
                                                            neuronsDistance = abs(motorUnitPools(poolIn)%unit(unitIn)%position_mm &
                                                                           - neuralTractPools(poolOut)%unit(unitOut)%position_mm)
                                                            weight = declineFactor / (declineFactor + neuronsDistance**2)
                                                            gmax = gmax * weight
                                                        end if
                                                        call motorUnitPools(poolIn)%unit(unitIn)%&
                                                        Compartments(compartmentIn)%SynapsesIn(synapseComp)%&
                                                        addConductance(gmax, delay, dyn, var, tau)
                                                        
                                                        if (allocated(neuralTractPools(poolOut)%&
                                                            unit(unitOut)%transmitSpikesThroughSynapses)) then
                                                            
                                                            
                                                            allocate(tempTransmitSpikes(size(neuralTractPools(poolOut)%&
                                                                unit(unitOut)%transmitSpikesThroughSynapses)))
                                                            
                                                            do j = 1, size(tempTransmitSpikes)
                                                                call tempTransmitSpikes(j)%assignSynapse(neuralTractPools(poolOut)%&
                                                                    unit(unitOut)%transmitSpikesThroughSynapses(j)%synapse)
                                                            end do
                                                            
                                                            
                                                            deallocate(neuralTractPools(poolOut)%&
                                                                    unit(unitOut)%transmitSpikesThroughSynapses)
                                                            
                                                            allocate(neuralTractPools(poolOut)%&
                                                                    unit(unitOut)%&
                                                                    transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1))
                                                            
                                                            do j = 1, size(tempTransmitSpikes)
                                                                call neuralTractPools(poolOut)%unit(unitOut)%&
                                                                transmitSpikesThroughSynapses(j)%&
                                                                assignSynapse(tempTransmitSpikes(j)%synapse)
                                                            end do

                                                            neuralTractPools(poolOut)%unit(unitOut)%&
                                                                    transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1)=&
                                                                    SynapsePointer()        

                                                            call neuralTractPools(poolOut)%unit(unitOut)%&
                                                            transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1)%&
                                                            assignSynapse(motorUnitPools(poolIn)%unit(unitIn)%&
                                                            Compartments(compartmentIn)%SynapsesIn(synapseComp))  
                                                            
                                                            deallocate(tempTransmitSpikes)
                                                        else 
                                                            allocate(neuralTractPools(poolOut)%&
                                                                    unit(unitOut)%&
                                                                    transmitSpikesThroughSynapses(1))
                                                            neuralTractPools(poolOut)%unit(unitOut)%&
                                                                    transmitSpikesThroughSynapses(1)=&
                                                                    SynapsePointer()        

                                                            call neuralTractPools(poolOut)%unit(unitOut)%&
                                                                    transmitSpikesThroughSynapses(1)%&
                                                                    assignSynapse(motorUnitPools(poolIn)%unit(unitIn)%&
                                                                    Compartments(compartmentIn)%&
                                                                    SynapsesIn(synapseComp))
                                                        end if            
                                                        
                                                        newIndex = size(motorUnitPools(poolIn)%&
                                                            unit(unitIn)%Compartments(compartmentIn)%&
                                                            SynapsesIn(synapseComp)%gmax_muS)
                                                        
                                                        call integerAddToList(neuralTractPools(poolOut)%&
                                                        unit(unitOut)%indicesOfSynapsesOnTarget,newIndex)
                                                        
                                                        numberOfSynapses = numberOfSynapses + 1
                                                    end if
                                                end do
                                            end if
                                        end if
                                    end do
                                end do
                            end if
                        end do
                    end do
                end do
            end do

            print '(A, I10, A)', 'All the ', numberOfSynapses,  ' synapses were built'


        ! ## Total number of synaptic noises in the system.
        ! self.numberOfSynapticNoise = 0

        ! NoiseSynapsesOut = conf.determineSynapses('Noise')
        ! for synapseIn in xrange(len(NoiseSynapsesOut)):
        !     pools[len(pools)] = SynapticNoise(conf, NoiseSynapsesOut[synapseIn][0])
        !     poolOut = len(pools) - 1
        !     gmax = float(conf.parameterSet('gmax:Noise>' + NoiseSynapsesOut[synapseIn][0] 
        !                                    + '-' + NoiseSynapsesOut[synapseIn][1]
        !                                    + '@' + NoiseSynapsesOut[synapseIn][2] + '|'
        !                                    + NoiseSynapsesOut[synapseIn][3],
        !                                    '', 0))
        !     delay = float(conf.parameterSet('delay:Noise>' + NoiseSynapsesOut[synapseIn][0]
        !                                     + '-' + NoiseSynapsesOut[synapseIn][1]
        !                                     + '@' + NoiseSynapsesOut[synapseIn][2] + '|'
        !                                     + NoiseSynapsesOut[synapseIn][3],
        !                                     '', 0))
        !     declineFactor = float(conf.parameterSet('dec:Noise>' + NoiseSynapsesOut[synapseIn][0]
        !                                             + '-' + NoiseSynapsesOut[synapseIn][1]
        !                                             + '@' + NoiseSynapsesOut[synapseIn][2] + '|'
        !                                             + NoiseSynapsesOut[synapseIn][3],
        !                                             '', 0))
        !     dyn = conf.parameterSet('dyn:Noise>' + NoiseSynapsesOut[synapseIn][0] 
        !                             + '-' + NoiseSynapsesOut[synapseIn][1]
        !                             + '@' + NoiseSynapsesOut[synapseIn][2] + '|'
        !                             + NoiseSynapsesOut[synapseIn][3],
        !                             '', 0)
        !     if dyn != 'None':
        !         var = float(conf.parameterSet('var:Noise>' + NoiseSynapsesOut[synapseIn][0]
        !                                       + '-' + NoiseSynapsesOut[synapseIn][1]
        !                                       + '@' + NoiseSynapsesOut[synapseIn][2] + '|' 
        !                                       + NoiseSynapsesOut[synapseIn][3],
        !                                       '', 0))
        !         tau = float(conf.parameterSet('tau:Noise>' + NoiseSynapsesOut[synapseIn][0]
        !                                       + '-' + NoiseSynapsesOut[synapseIn][1]
        !                                       + '@' + NoiseSynapsesOut[synapseIn][2]
        !                                       + '|' + NoiseSynapsesOut[synapseIn][3],
        !                                       '', 0))
        !     else:
        !         var = 0
        !         tau = 10000
        !     for unitOut in xrange(len(pools[poolOut].unit)):
        !         for poolIn in xrange(len(pools)):
        !             if NoiseSynapsesOut[synapseIn][0] == pools[poolIn].pool and pools[poolIn].kind != 'SN':
        !                 for unitIn in xrange(len(pools[poolIn].unit)):
        !                     for compartmentIn in xrange(len(pools[poolIn].unit[unitIn].compartment)):
        !                         if NoiseSynapsesOut[synapseIn][1] == pools[poolIn].unit[unitIn].kind and NoiseSynapsesOut[synapseIn][2] == pools[poolIn].unit[unitIn].compartment[compartmentIn].kind and pools[poolIn].unit[unitIn].index == pools[poolOut].unit[unitOut].index:
        !                             do synapse = 1, size(motorUnitPools(poolIn).unit[unitIn].compartment[compartmentIn].SynapsesIn)
        !                                 if pools[poolIn].unit[unitIn].compartment[compartmentIn].SynapsesIn[synapse].kind == NoiseSynapsesOut[synapseIn][3]:
        !                                     if np.isfinite(declineFactor):
        !                                         neuronsDistance = np.abs(pools[poolIn].unit[unitIn].position_mm
        !                                                                     - pools[poolOut].unit[unitOut].position_mm)
        !                                         weight = declineFactor / (declineFactor + neuronsDistance**2)
        !                                     else:
        !                                         weight = 1
        !                                     pools[poolIn].unit[unitIn].compartment[compartmentIn].SynapsesIn[synapse].addConductance(gmax*weight, delay, dyn, var, tau)
        !                                     pools[poolOut].unit[unitOut].transmitSpikesThroughSynapses.append(pools[poolIn].unit[unitIn].compartment[compartmentIn].SynapsesIn[synapse])                                                            
        !                                     pools[poolOut].unit[unitOut].indicesOfSynapsesOnTarget.append(len(pools[poolIn].unit[unitIn].compartment[compartmentIn].SynapsesIn[synapse].gmax_muS) - 1)
        !                                     self.numberOfSynapticNoise += 1
                  

        

        !             print 'All the ' + str(self.numberOfSynapticNoise) +  ' synaptic noises were built'                   
        



        end subroutine

end module SynapsesFactoryModule


    


    
        