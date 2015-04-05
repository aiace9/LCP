! rw2d.f90
! A simple random walk program  in 2D. 

program drunk
  
  implicit none

  !parameter
  real, parameter :: step=1.0, twopi=2.0*3.1415926 ! step size and constants
  character(len=15), save :: format1 = "(1i5,1x,2F14.7)"

  
  integer :: i, j, steps, rdrops, sizer
  real :: heigh, rnd
  integer, dimension(:), allocatable :: seed
  real, dimension(:,:), allocatable :: r_rdrops
  real, dimension(:), allocatable :: t_rdrops
  !support variables to calculate the mean square displacement
  real :: dr2, xn2, yn2, xn, yn

  !error flags
  integer :: ios

  ! initialization
  rdrops = 0
  dr2 = 0
  xn2 = 0
  yn2 = 0
  xn = 0
  yn = 0

  !random initialization
  call random_seed(sizer)

  
  print*, "Enter the number of rain drops:"
  read*, rdrops
  print*, "the heigh in meters: (3000-4000)"
  read*, heigh 

  !open statement
  open(unit=1, file='pos.dat', iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file pos.dat"

  open(unit=2, file='g.dat', iostat=ios, status="unknown", action="write")
  if ( ios /= 0 ) stop "Error opening file f.dat"

  !allocation statements
  allocate(seed(sizer))
  allocate(r_rdrops(rdrops,2))
  allocate(t_rdrops(rdrops))
  r_rdrops(:,1) = 0
  r_rdrops(:,2) = heigh / 5 ! assuming a real velocity of 5 m/s and 1 second per step
  
  print *,'Here the seed has ',sizer,' components; insert them (or print "/") >'  
  read(*,*)seed

  call random_seed(PUT=seed)

  i = 0
  
  do i = 1, rdrops, 1
    j = 0    
    do  
      call random_number(rnd)
      select case (floor(rnd * 100))
        case ( 0 : 14 )
          r_rdrops(i,1)=r_rdrops(i,1) + 1
        case ( 15 : 29)
          r_rdrops(i,1)=r_rdrops(i,1) - 1 
        case ( 30 : 39 )
          r_rdrops(i,2)=r_rdrops(i,2) + 1
        case ( 40 : 99) 
          r_rdrops(i,2)=r_rdrops(i,2) - 1
        case default
            print*, 'error'            
      end select
      j = j + 1
      if ( r_rdrops(i,2) <= 0 ) exit
    end do
    
    t_rdrops(i) = j

    !write(unit=1, fmt=format1, iostat=ios) i, r_drunks(i,:)
    !if ( ios /= 0 ) stop "Write error in file unit 1"
  end do

  do i = 1, rdrops, 1
    
    xn2 = xn2 + r_rdrops(i,1)**2 !dotproduct can be use as improvement
    xn = xn + r_rdrops(i,1)

  end do

  dr2 = xn2/rdrops  + (xn/rdrops) ** 2
  
  print*, "the medium time take from the rain drops to reach the ground is", sum(t_rdrops)/real(rdrops), "s"
  print*, " <dx^2> value ", dr2
  print*, "the vertical average velocity is: ", heigh * real(rdrops)  / sum(t_rdrops) * 3.6, "km/h"

  !write(unit=2, fmt=*, iostat=ios) steps, dr2
  !if ( ios /= 0 ) stop "Write error in file unit 2"
  

  close(unit=1, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 1"
  close(unit=2, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 2"
  
  
  deallocate(seed, r_rdrops, t_rdrops)
 
end program drunk
