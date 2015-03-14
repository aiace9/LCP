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
  integer :: number, k

  real (kind=dp), dimension(:,:), allocatable :: rsequance
  real (kind=dp), dimension(:,:), allocatable :: test
  
  !error variable and flag
  integer :: flag = 0
  integer :: err,ios

  !     supply initial values of some variables:
  !     seed:   to start; a: 
  !     number: how many numbers we want to generate   

  print*,' wich range of N do you want to simulate?'
  read(*,*)number

  print*, 'wich range of K do you want to simulate?'
  read(*,*)k

  !inizialising the serching tool for find the period
  allocate(rsequance(number,number), stat=err)
  if (err /= 0) print *, "rsequance: Allocation request denied"
  rsequance=0

  allocate(test(k,number), stat=err)
  if (err /= 0) print *, "test: Allocation request denied"
  test=0
  

  !
  open(unit=1, file="data.dat", iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file data.dat"

  print*, "starting generation of random sequance"
  
  !generation of the numbers

  do i = 1, number, 1
      do j = 1, i , 1
        call random_number(rsequance(i,j))
      end do
  end do

  !Test for b point
  ! cicle over the k
  do i = 1, k, 1
    !cicle over the sequances of N:  1, 2, 3,...., N
    do j = 1, number, 1
      !sum the element of the sequance
      print*, 'start',j
      do z = 1, j, 1
        !print*, test(i,j)
        test(i,j) = test(i,j) + (rsequance(j,z)**i)
      end do
      !division over the number of smapes
      test (i,j) = test (i,j) / real(j) - 1 / (i+1)
    end do
  end do
  print*, 'a'
  
  do i = 1, number, 1
    write(unit=1, fmt=*, iostat=ios) i, test(7,i)
    if ( ios /= 0 ) stop "Write error in file unit 1"
  end do

  !deallocation stuff
  
  if (allocated(rsequance)) deallocate(rsequance, stat=err)
  if (err /= 0) print *, "rsequance: Deallocation request denied"
  
  if (allocated(test)) deallocate(test, stat=err)
  if (err /= 0) print *, "test: Deallocation request denied"  
  
  close(unit=1, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 1"
  
  print*, "Deallocation complete"
  stop
end program random_lcm
