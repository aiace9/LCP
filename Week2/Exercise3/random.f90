!
!     generation of pseudorandom numbers - linear congruential method
!
program random_lcm
  implicit  none
  !
  !     declaration of variables
  integer, parameter :: dp = selected_real_kind(14,200)
  integer :: i,j, number, column
  integer, dimension(:), allocatable :: histogram
 
  real (kind=dp), dimension(:), allocatable :: rsequance
  real(kind=dp) :: dr  
  
  !error variable and flag
  integer :: flag = 0
  integer :: err,ios

  !     supply initial values of some variables:
  !     seed:   to start; a: 
  !     number: how many numbers we want to generate   

  print*,' How many numbers do you want to generate ?'
  read(*,*)number

  print*, 'How many columns do you want in the plot ?'
  read(*,*)column

  !inizialising the serching tool for find the period
  allocate(rsequance(number), stat=err)
  if (err /= 0) print *, "rsequance: Allocation request denied"
  rsequance=0

  allocate(histogram(column), stat=err)
  if (err /= 0) print *, "histogram: Allocation request denied"
  histogram = 0
  !
  OPEN(unit=1, file="random.dat", status="replace", action="write")
  
  open(unit=2, file="ist.dat", iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file ist.dat"
  
  open(unit=3, file="xy.dat", iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file xy.dat"
  
  dr = 1._dp/column

  print*, "------"

  !generation of the number
  do i = 1, number , 1
     call random_number(rsequance(i))
     WRITE (unit=1,fmt=*) i, rsequance(i)
  end do

  print*,'data saved in random.dat'

  !create histogram.
  do i = 1, number, 1
    j= int(rsequance(i) / dr) + 1
    if ( j <= column ) then
      histogram(j) = histogram(j) + 1
    else
      print*, "Error during column filling"
    end if
  end do

  do i = 1, column, 1
    write(unit=2, fmt=*, iostat=ios) i, histogram(i)
    if ( ios /= 0 ) stop "Write error in file unit 2"  
  end do

  print*, "data saved in ist.dat"

  do i = 1, number, 2
    write(unit=3, fmt=*, iostat=ios) rsequance(i), rsequance (i+1)
    if ( ios /= 0 ) stop "Write error in file unit 3"    
  end do

  print*, "data saved in xy.dat"

  print*, "------"

  
  
  !deallocation stuff
  
  if (allocated(rsequance)) deallocate(rsequance, stat=err)
  if (err /= 0) print *, "rsequance: Deallocation request denied"
  
  if (allocated(histogram)) deallocate(histogram, stat=err)
  if (err /= 0) print *, "histogram: Deallocation request denied"
  
  close(unit=1, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 1"
  
  close(unit=2, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 2"
  
  close(unit=3, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 3"
  
  print*, "Deallocation complete"
  stop
end program random_lcm
