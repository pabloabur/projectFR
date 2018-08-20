module UtilitiesMod
    use mkl_spblas
    use DynamicalArrays
    implicit none
    integer, parameter :: wp = kind(1.0d0)
    private :: wp

    contains

        subroutine convertToMKLSparse(A, spA)
            real(wp), intent(in) :: A(:,:)
            type(sparse_matrix_t), intent(inout) :: spA
            integer :: spIndexing, spRows, spCols, spNumberOfElements
            integer, dimension(:), allocatable :: spRowStart, spRowEnd, spColIdx
            real(wp), dimension(:), allocatable :: spValues
            integer :: i, j, stat

            spIndexing = 1
            spRows = size(A,1)
            spCols = size(A,2)
            print '(I4)', spRows
            spNumberOfElements = 0
            allocate(spRowStart(spRows))
            allocate(spRowEnd(spRows))
            spRowStart(:) = 0
            do i = 1, spRows
                do j = 1, spCols
                    if (abs(A(i,j))>1e-10) then
                        spNumberOfElements = spNumberOfElements + 1
                        if (spRowStart(i) == 0) then
                            spRowStart(i) = spNumberOfElements
                            spRowEnd(i-1) = spNumberOfElements
                        end if                    
                        call AddToList(spValues, A(i,j))
                        call integerAddToList(spColIdx, j)
                    end if                
                end do
            end do
            spRowEnd(spRows) = spNumberOfElements + 1
           
            stat = mkl_sparse_d_create_csr(spA, &
                                           spIndexing, &
                                           spRows, &
                                           spCols, &
                                           spRowStart, &
                                           spRowEnd, &
                                           spColIdx, &
                                           spValues)
            print '(I1)', stat

        end subroutine convertToMKLSparse

end module UtilitiesMod