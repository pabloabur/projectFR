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

module MusclePointerNoChannelClass
    use MotorUnitPoolNoChannelClass
    implicit none
    private
    public :: MusclePointerNoChannel

    type MusclePointerNoChannel
        type(MotorUnitPoolNoChannel), pointer :: muscle

        contains
            procedure :: assignMuscle

    end type MusclePointerNoChannel

    interface MusclePointerNoChannel
        module procedure init_MusclePointer
    end interface MusclePointerNoChannel

    contains

        type(MusclePointerNoChannel) function init_MusclePointer()
            if (associated(init_MusclePointer%muscle)) nullify(init_MusclePointer%muscle) 
        end function

        subroutine assignMuscle(self, newMuscle)
            class(MusclePointerNoChannel), intent(inout) :: self
            class(MotorUnitPoolNoChannel), intent(in), target :: newMuscle

            self%muscle => newMuscle
        end subroutine

end module MusclePointerNoChannelClass