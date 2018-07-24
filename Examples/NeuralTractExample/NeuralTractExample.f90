program NeuralTractUnitExample
    !use NeuralTractClass
    use ConfigurationClass
    use ogpf 
    use randomSeedInitialize
    implicit none 
    !integer, parameter :: wp = kind(1.0d0)
    type(Configuration) :: conf
    real(wp), parameter :: ISI_ms = 12.0
    real(wp), parameter :: FR_Hz = 1000.0/ISI_ms
    integer, parameter :: GammaOrder = 10
    real(wp), parameter :: PI = 4 * atan(1.0_wp)    
    !use, intrinsic :: iso_fortran_env, only : wp => real64
    real(wp), parameter :: dt = 0.05
    real(wp), parameter :: tf = 10000
    integer, parameter :: timeLength = int(tf/dt)
    integer :: i
    real(wp), dimension(timeLength) :: t
    real(wp) :: tic, toc
    integer :: index = 2
    type(gpf) :: gp
    real(wp), dimension(:), allocatable :: indices
    real(wp), dimension(:), allocatable :: ISI
    real(wp), dimension(:), allocatable :: ISIhist, ISIhistLimits
    real(wp) :: histogramMax
    real(wp) :: M, SD, CV, MFR, SD_FR
    character(80) :: line1
    character(80) :: param1, param2, param3
    integer :: il, ierr, j
    integer :: stop1
    real(wp) :: param2Real
    character(len = 6) :: pool = 'SOL'
    character(len = 80) :: filename = 'confNeuralTractExample.rmto'
    character(len=30)::paramTag
    real(wp) :: parammReal
    integer :: parammInt
    character(len = 80) :: paramChar
    

    conf = Configuration(filename)

    
    t = [(dt*i, i=1, timeLength)]
   
    call init_random_seed()

    print '(F15.6)', conf%timeStep_ms
    paramTag = 'timeStep'
    paramChar = conf%parameterSet(paramTag, pool, index)
    print '(A)', paramChar
    read(paramChar, *)parammReal
    parammInt = int(parammReal)
    print '(I3)', parammInt
    ! open(1,file = 'confNeuralTractExample.rmto', status='old',iostat=ierr)

    ! do while (ierr.eq.0)
    !     read(1, '(A)',iostat=ierr) line1
    !     il=len_trim(line1)
    !     !print *, line1(1:il)
    !     j = 1
    !     do i = 1, il
    !         if (line1(i:i) == ',') then
    !             if (j.eq.1) then 
    !                 param1 = line1(1:i-1)
    !                 j = j + 1
    !                 stop1 = i
    !             else if (j.eq.2) then 
    !                 param2 = line1(stop1+1:i-1)
    !                 param3 = line1(i+1:il)
    !             end if 
    !         end if            
    !     end do
    !     print '(A, A, A)',  param1, param2, param3 
    !     print '(I3, I3, I3)', len_trim(param1), len_trim(param2), len_trim(param3)
    !     print '(A)', param2(1:len_trim(param2))

    !     !read(param2(1:len_trim(param2)), *)param2Real
    !     !print '(F15.6)', param2Real
    ! end do

    ! close(unit = 1)

    ! open(1,file = 'confNeuralTractExample.rmto', status='old',iostat=ierr)

    ! do while (ierr.eq.0)
    !     read(1, '(A)',iostat=ierr) line1
    !     il=len_trim(line1)
    !     !print *, line1(1:il)
    !     j = 1
    !     do i = 1, il
    !         if (line1(i:i) == ',') then
    !             if (j.eq.1) then 
    !                 param1 = line1(1:i-1)
    !                 j = j + 1
    !                 stop1 = i
    !             else if (j.eq.2) then 
    !                 param2 = line1(stop1+1:i-1)
    !                 param3 = line1(i+1:il)
    !             end if 
    !         end if            
    !     end do
    !     print '(A, A, A)',  param1, param2, param3 
    !     print '(I3, I3, I3)', len_trim(param1), len_trim(param2), len_trim(param3)
    !     print '(A)', param2(1:len_trim(param2))

    !     !read(param2(1:len_trim(param2)), *)param2Real
    !     !print '(F15.6)', param2Real
    ! end do
    ! pool = 'SOL'
    ! line1 = trim(pool) // '-F'
    ! if (line1.eq.('SOL'//'-F')) print "(A)", line1

    ! print '(I3)', len(trim(line1))




    ! call cpu_time(tic)
    ! do i = 1, timeLength
    !     call NT_unit%atualizeNeuralTractUnit(t(i), FR_Hz*dt/1000.0, GammaOrder)
    ! enddo
    ! call cpu_time(toc)
    ! print "(E15.6, A)", toc - tic,  " seconds"

    

    ! allocate(indices(size(NT_unit%spikesGenerator%points)))
    
    ! do i = 1, size(NT_unit%spikesGenerator%points)
    !     indices(i) = index    
    ! enddo
    
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! call gp%title('Neural tract spikes instants')
    ! call gp%xlabel('t (s))')
    ! call gp%ylabel('Descending command index')
    ! call gp%plot(NT_unit%spikesGenerator%points/1000, indices, 'with points pt 5 lc rgb "#0008B0"')
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! allocate(ISI(size(NT_unit%spikesGenerator%points)-1))
    ! do i = 1, size(ISI)
    !     ISI(i) = NT_unit%spikesGenerator%points(i+1) - NT_unit%spikesGenerator%points(i) 
    ! end do

    ! allocate(ISIhist(10))
    ! allocate(ISIhistLimits(10))
    ! histogramMax = maxval(ISI)
    
    ! do i = 1, size(ISIhist)
    !     ISIhist(i) = count(ISI.gt.((i-1)*histogramMax/10.0).and.ISI.lt.(i*histogramMax/10.0))
    !     ISIhistLimits(i) = (i)*histogramMax/10.0
    ! end do

    ! call gp%title('ISI histogram')
    ! call gp%xlabel('ISI')
    ! call gp%ylabel('count')

    ! call gp%plot(ISIhistLimits,ISIhist, 'with impulses lw 2.5', &
    !         x2=ISIhistLimits, y2=ISIhist,  ls2='with points pt 7')
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! M = sum(ISI)/size(ISI)
    ! SD = SQRT(sum((ISI-M)**2)/size(ISI))
    ! CV = SD / M
    ! MFR = 1000/M
    
    ! print "(A, F15.6, A)", 'Mean = ', M, ' ms'
    ! print "(A, F15.6, A)", 'STD = ', SD, ' ms'
    ! print "(A, F15.6)", 'CV = ', CV


end program NeuralTractUnitExample