!
!     generation of pseudorandom numbers - linear congruential method
!
program random_lcm
  implicit  none
  !
  !     declaration of variables
  integer :: i,j, number, seed, x, a,m,c
  integer, dimension(:), allocatable :: sequance
  integer :: flag = 0
  
  !error variable
  integer :: err
  !     supply initial values of some variables:
  !     seed:   to start; a: 
  !     number: how many numbers we want to generate   
  print*,'Use LCM: x_(i+1)=mod(a*x_i+c,m)'
  print*,'Insert seed (=x_0), a, m, c >'
  read (*,*) seed, a, m, c  

  print*,' How many numbers do you want to generate ?'
  read(*,*)number

  !inizialising the serching tool for find the period
  allocate(sequance(number), stat=err)
  if (err /= 0) print *, "sequance: Allocation request denied"
  sequance=0
  !
  OPEN(unit=1, file="random.dat", status="replace", action="write")
  sequance(1) = seed

  do i = 2, number , 1
     x = mod ((a*sequance(i-1)+c), m)
     sequance(i) = x
     WRITE (unit=1,fmt=*) i, x
  end do

  print*, "------"

  do i = 1, number, 1
    do j = i+1, number, 1
      if ( sequance(i) == sequance(j) ) then
        print*, 'periodo trovato'
        print*, abs(i-j)
        flag = 1
        exit
      end if
    end do
    if ( flag == 1 ) exit
  end do



  !print conclusion
  print*, minval(sequance),maxval(sequance)
  print*,' data saved in random.dat'
  
  !deallocation stuff
  if (allocated(sequance)) deallocate(sequance, stat=err)
  if (err /= 0) print *, "sequance: Deallocation request denied"
  close(1)

  stop
end program random_lcm
