! '''
!     Neuromuscular simulator in Fortran.
!     Copyright (C) 2018  Renato Naville Watanabe
!                         
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

module QueueClass
    implicit none
    private
    public :: Queue

    type Queue
        integer, dimension(:), allocatable :: item
        integer :: endQueue
        

        contains
            procedure :: append
            procedure :: popleft
            procedure :: remove
            procedure :: extend
            procedure :: clear
    end type Queue

    interface Queue
        module procedure init_Queue
    end interface Queue

    contains

        type(Queue) function init_Queue()
            if (allocated(init_Queue%item)) deallocate(init_Queue%item)
            init_Queue%endQueue = 0
        end function

        subroutine append(self, newItem)
            class(Queue), intent(inout) :: self
            integer, intent(in) :: newItem
            integer, dimension(:), allocatable :: tempList

            if (allocated(tempList)) deallocate(tempList)
                       
            if (self%endQueue > 0) then
                allocate(tempList(self%endQueue))
                tempList = self%item
                
                deallocate(self%item)
                
                allocate(self%item(self%endQueue+1))
                
                self%item(1:self%endQueue) = tempList(:)
                self%item(self%endQueue+1) = newItem
                self%endQueue = self%endQueue + 1
            else
                allocate(self%item(1))
                self%item(1) = newItem
                self%endQueue = 1
            end if
        end subroutine

        integer function popleft(self) result(item)
            class(Queue), intent(inout) :: self
            integer, dimension(:), allocatable :: tempList

            item = self%item(1)
            if (allocated(tempList)) deallocate(tempList)
            self%endQueue = self%endQueue - 1
            if (self%endQueue > 0) then
                allocate(tempList(self%endQueue))
                tempList(:) = self%item(2:self%endQueue+1)
                deallocate(self%item)
                allocate(self%item(self%endQueue))
                self%item(:) = tempList(:)
            else
                deallocate(self%item)
            end if
        end function

        subroutine remove(self, removeItem)
            class(Queue), intent(inout) :: self
            integer, intent(in) :: removeItem
            integer, dimension(:), allocatable :: tempList1, tempList2
            integer :: index, i

            index = -1
            i = 1
            do while ((index == -1).and.(i<=self%endQueue))
                i = i + 1
                if (self%item(i-1) == removeItem) index = i-1
            end do

            if (index > 0) then
                if (allocated(templist1)) deallocate(tempList1)
                if (allocated(templist2)) deallocate(tempList2)

                allocate(tempList1(index - 1))
                allocate(tempList2(self%endQueue - index))
                tempList1(:) = self%item(1:index-1)
                tempList2(:) = self%item(index+1:self%endQueue)
                
                deallocate(self%item)
                allocate(self%item(self%endQueue-1))

                self%item(1:index-1) = tempList1(:)
                self%item(index:self%endQueue-1) = tempList2
                self%endQueue = self%endQueue - 1
            end if



        end subroutine

        subroutine extend(self, newItens)
            class(Queue), intent(inout) :: self
            integer, intent(in) :: newItens(:)
            integer, dimension(:), allocatable :: tempList
            integer :: sizeNewList

            sizeNewList = size(newItens)

            if (allocated(tempList)) deallocate(tempList)
            
            if (self%endQueue > 0) then
                allocate(tempList(self%endQueue))
                tempList = self%item
                
                deallocate(self%item)
                
                allocate(self%item(self%endQueue+sizeNewList))
                
                self%item(1:self%endQueue) = tempList
                self%item(self%endQueue+1:self%endQueue+sizeNewList) = newItens
                self%endQueue = self%endQueue + sizeNewList
            else
                allocate(self%item(sizeNewList))
                self%item = newItens
                self%endQueue = sizeNewList
            end if
        end subroutine

        subroutine clear(self)
            class(Queue), intent(inout) :: self
            
            deallocate(self%item)
            self%endQueue = 0            
        end subroutine
end module QueueClass