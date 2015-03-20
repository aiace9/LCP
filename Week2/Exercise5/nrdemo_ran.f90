program demo

  use ran_module
  implicit none
  integer :: i,idum
  real :: x

  print*, "idum (<0) = "
  read*,idum
  x =ran_func(idum)
  print*,"Random number: ",x

  do i=1,10
     x = ran_func(idum)
     print*,"Random number: ",x
  end do

end program demo


