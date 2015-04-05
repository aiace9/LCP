! rw2d.f90
! A simple random walk program  in 2D. 

program drunk
  
  implicit none

  !parameter
  real, parameter :: step=1.0, twopi=2.0*3.1415926 ! step size and constants
  character(len=15), save :: format1 = "(1i5,1x,2F14.7)"

  
  integer :: i, j, steps, drunks, sizer
  real :: phi, rnd
  integer, dimension(:), allocatable :: seed
  real, dimension(:,:), allocatable :: r_drunks
  !support variables to calculate the mean square displacement
  real :: dr2, xn2, yn2, xn, yn

  !error flags
  integer :: ios

  ! initialization
  drunks = 0
  dr2 = 0
  xn2 = 0
  yn2 = 0
  xn = 0
  yn = 0

  !random initialization
  call random_seed(sizer)

  
  print*, "Enter the number of drunks"
  read*, drunks
  print*, "Enter number of steps:"
  read*, steps

  !open statement
  open(unit=1, file='pos.dat', iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file pos.dat"

  open(unit=2, file='g.dat', iostat=ios, status="unknown", action="write", position = "append")
  if ( ios /= 0 ) stop "Error opening file f.dat"

  !allocation statements
  allocate(seed(sizer))
  allocate(r_drunks(drunks,2))
  r_drunks = 0
  
  print *,'Here the seed has ',sizer,' components; insert them (or print "/") >'  
  read(*,*)seed

  call random_seed(PUT=seed)

  i = 0
  
  do i = 1, drunks, 1
    
    do  j=1, steps
      call random_number(rnd)
      select case (floor(rnd * 100))
        case ( 0 : 39 )
          r_drunks(i,1)=r_drunks(i,1) + 1
        case ( 40 : 59)
          r_drunks(i,1)=r_drunks(i,1) - 1 
        case ( 60 : 79 )
          r_drunks(i,2)=r_drunks(i,2) + 1
        case ( 80 : 99) 
          r_drunks(i,2)=r_drunks(i,2) - 1
        case default
            print*, 'error'            
      end select
    end do

    write(unit=1, fmt=format1, iostat=ios) i, r_drunks(i,:)
    if ( ios /= 0 ) stop "Write error in file unit 1"
  end do

  do i = 1, drunks, 1
    
    xn2 = xn2 + r_drunks(i,1)**2 !dotproduct can be use as improvement
    yn2 = yn2 + r_drunks(i,2)**2
    xn = xn + r_drunks(i,1)
    yn = yn + r_drunks(i,2)

  end do

  dr2 = xn2/drunks + yn2/drunks + (xn/drunks) ** 2 + (yn/drunks) ** 2
  
  print*, "the first motion sequence result in"
  print*, " <dr^2> value ", dr2
  print*, 'values for the plot:', steps, dr2

  write(unit=2, fmt=*, iostat=ios) steps, dr2
  if ( ios /= 0 ) stop "Write error in file unit 2"
  

  close(unit=1, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 1"
  close(unit=2, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 2"
  
  
  deallocate(seed, r_drunks)
 
end program drunk
