  module DynamicalArrays

  contains

      subroutine AddToList(list, element)

          implicit none

          integer :: i, isize
          double precision, intent(in) :: element
          double precision, dimension(:), allocatable, intent(inout) :: list
          double precision, dimension(:), allocatable :: clist


          if(allocated(list)) then
              isize = size(list)
              allocate(clist(isize+1))
              do i=1,isize          
              clist(i) = list(i)
              end do
              clist(isize+1) = element

              deallocate(list)
              
              allocate(list(isize+1))
              
              do i=1,isize + 1
                list(i) = clist(i)
              end do
              deallocate(clist)

          else
              allocate(list(1))
              list(1) = element
          end if


      end subroutine AddToList

      subroutine integerAddToList(list, element)

          implicit none

          integer :: i, isize
          integer, intent(in) :: element
          integer, dimension(:), allocatable, intent(inout) :: list
          integer, dimension(:), allocatable :: clist


          if(allocated(list)) then
              isize = size(list)
              allocate(clist(isize+1))
              do i=1,isize          
                clist(i) = list(i)
              end do

              clist(isize+1) = element

              deallocate(list)
              allocate(list(isize+1))
              
              do i=1,isize + 1
                list(i) = clist(i)
              end do
              deallocate(clist)
          else
              allocate(list(1))
              list(1) = element
          end if


      end subroutine integerAddToList
  

  end module DynamicalArrays