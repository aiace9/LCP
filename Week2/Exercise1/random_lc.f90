!
!     generation of pseudorandom numbers - linear congruential method
!
program random_lcm
  implicit  none
  !
  !     declaration of variables
  integer :: i, number, old, seed, x, a,m,c
  integer :: max = 0, min
  ! 
  !     supply initial values of some variables:
  !     seed:   to start; a: 
  !     number: how many numbers we want to generate   
  print*,'Use LCM: x_(i+1)=mod(a*x_i+c,m)'
  print*,'Insert seed (=x_0), a, m, c >'
  read (*,*) seed, a, m, c

  min = m !inizialising the minimum number at the top of the sequance.

  print*,' How many numbers do you want to generate ?'
  read(*,*)number
  !
  OPEN(unit=1, file="random.dat", status="replace", action="write")
  old = seed
  !
  do i = 1, number
     x = mod ((a*old+c), m)
     if ( x > max ) max = x
     if ( x < min ) min = x
     if ( x==5 ) print*, i
     WRITE (unit=1,fmt=*) i, x
     old = x
  end do
  close(1)
  print*, min,max
  print*,' data saved in random.dat'
  stop
end program random_lcm
