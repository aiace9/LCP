!
!     generation of pseudorandom numbers
!     with the uniformity test
!
program random_lcm
  implicit  none
  !
  !     declaration of variables
  integer, parameter :: dp = selected_real_kind(14,200)
  ! variabili cicliche
  integer :: i,j,z
  integer :: number, k, x

  real(dp), dimension(:), allocatable :: rvalue
  real (kind=dp), dimension(:), allocatable :: test
  
  !error variable and flag
  integer :: err,ios
  logical :: debug = .true.

  !     supply initial values of some variables:
  !     seed:   to start; a: 
  !     number: how many numbers we want to generate   

  print*,' wich range of N do you want to simulate?'
  read(*,*)number

  print*, 'wich K do you want to simulate?'
  read(*,*)k

  print*, 'at wich distatce do you want to compute the correlation?'
  read(*,*)x

  allocate(rvalue(number), stat=err)
  if (err /= 0) print *, "rvalue: Allocation request denied"
  rvalue = 0

  allocate(test(number), stat=err)
  if (err /= 0) print *, "test: Allocation request denied"
  test = 0
  

  !
  open(unit=1, file="uniformity.dat", iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file data.dat"

  open(unit=2, file="correlation.dat", iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file data.dat"
  print*, "starting generation of random sequance"
  
  !generation of the numbers

  call random_number(rvalue)

  if ( debug )  print*, 'debug1'

  !uniformity test

  do i = 1, number, 1
    do j = 1, i, 1
        test(i) = test (i) + rvalue(j)**k 
    end do
    test(i) = test (i) / real(i) - 1._dp / (k+1)
  end do
  
  do i = 1, number, 1
    write(unit=1, fmt=*, iostat=ios) i, test(i)
    if ( ios /= 0 ) stop "Write error in file unit 1"
  end do

  if ( debug ) print*, 'debug2'

  !correlation test
  test = 0

  do i = 1, number - x, 1
    do j = 1, i, 1
      test(i) = test(i) + rvalue(i) * rvalue(i + x)
    end do
    test(i) = test(i) / real(i) - 1 / 4._dp
  end do

  do i = 1, number, 1
    write(unit=2, fmt=*, iostat=ios) i, test(i)
    if ( ios /= 0 ) stop "Write error in file unit 2"
  end do

  !deallocation stuff

  if (allocated(rvalue)) deallocate(rvalue, stat=err)
  if (err /= 0) print *, "rvalue: Deallocation request denied"
  
  if (allocated(test)) deallocate(test, stat=err)
  if (err /= 0) print *, "test: Deallocation request denied"  
  
  close(unit=1, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 1"
  
  close(unit=2, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 2"

  print*, "Deallocation complete"
  stop
end program random_lcm
