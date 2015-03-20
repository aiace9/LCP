program rantest_intrinsic
!
! test program, call to intrinsic f90 random number generator
!  generate random numbers in [0,1[ ; then,
!  generate random integers between n_min and n_max.
!
  implicit none
  real :: rnd
  integer, dimension (:), allocatable :: x
  integer :: L,i,j,n_min,n_max
  integer :: err,ios

  logical :: flag = .false.
  logical :: debug = .false.

  ! generates L random numbers in [0,1[
  print*,' Insert the length of the sequence >'
  read(*,*)L        ! length of sequence

  !  generates integer random numbers between n_min and n_max
  print*,' insert n_min, n_max >'
  read(*,*)n_min,n_max

  
  if ( debug ) print*, 'debug1'
  
  allocate(x(L), stat=err)
  if (err /= 0) print *, "x: Allocation request denied"
  
  if ( debug ) print*, 'debug2'
  
  open(unit=2, file="data.dat", iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file data.dat"
  
  if ( debug ) print*, 'debug3' 
  
  do i = 1, L, 1
    call random_number(rnd)
    x(i) = (n_max - n_min + 1)*rnd + n_min
    if ( debug ) print*, 'debug = ', x(i)
  end do

  if ( debug ) print*, 'debug4'

  do i = 1, L, 1
    if ( debug ) print*, 'debug = ', i
    write(unit=2, fmt=*, iostat=ios) i, x(i)
    if ( ios /= 0 ) stop "Write error in file unit 1" 
  end do

  if ( debug ) print*, 'done, starting period'

  !period searcing
  do i = 1, L, 1
    do j = i+1, L, 1
      if ( x(i) == x(j) ) then
        print*, 'periodo trovato start:', i,' end:' , j
        print*, abs(i-j)
        flag = .true.
        exit
      end if
    end do
    if ( flag ) exit
  end do
  if (flag .eqv. .false.) print*, 'periodo non trovato'


  close(unit=2, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 1"
  
  if ( debug ) print*, 'done'

  deallocate(x, stat=err)
  if (err /= 0) print *, "x: Deallocation request denied"

end program rantest_intrinsic

