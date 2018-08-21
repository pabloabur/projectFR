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

module StringClass
    implicit none
    private
    public :: String

    type String
        character(len = 80) :: string

       contains     
    end type String

    interface String
        module procedure init_String
    end interface String
    
    contains

        type (String) function init_String(newString)
            character(*), intent(in) :: newString 
            init_String%string = trim(newString)
        end function

end module StringClass
