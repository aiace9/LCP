!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! gauss_metropolis.f90
!
! METROPOLIS generation of random numbers with a Gaussian distribution
! P(x) = exp(-x**2/(2*sigma**2))/sqrt(2*pi*sigma**2)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program gauss_metropolis
  implicit none
  integer, parameter :: dp=selected_real_kind(13)
  integer :: i,n,ibin,maxbin,m
  real(kind=dp):: sigma,rnd,delta,x0
  real(kind=dp):: x,x1,x2,x3,x4,xp,expx,expxp,w,acc
  character(len=13), save :: format1 = "(a7,2x,2f9.5)"

  integer :: ios
  integer :: n_min,n_max, n_step
  integer :: j



  !print*,' insert n, sigma, x0, delta, max number of bin in the histogram  >'
  !read*, n, sigma,x0,delta,maxbin
  !Point 2 
  !n=10000
  sigma=1
  x0=0
  delta = 5
  !maxbin=500
  !read*,delta


   n_max = 10000  
   acc = 0.0_dp
   x = x0
   x2 = 0.0_dp
  
   write(unit=*,fmt=*)"# n_max, x0, delta = ",n,x0,delta
   
   do i=1,n_max
      x2 = x2 + x**2

      !ccccccccccccccccccccccccccccccc
      expx = - x**2 /(2*sigma**2)    !
      call random_number(rnd)        !
      xp = x + delta * (rnd-0.5_dp)  !
      expxp = - xp**2 /(2*sigma**2)  !   metropolis
      w = exp (expxp-expx)           !   algorithm
      call random_number(rnd)        !
      if (w > rnd) then              !
         x = xp                      !
      !ccccccccccccccccccccccccccccccc
         acc=acc+1.0_dp                
      endif
    write(unit=*,fmt=*)"---------------"
    write(unit=*,fmt=*)"# Status"
    write(unit=*,fmt=*)"# ", abs((x2/real(i)-sigma**2)/sigma**2) *100
    if (abs((x2/real(i)-sigma**2)/sigma**2 *100) < 5.d0) then
      print*, 'achievement at ', i, 'steps'
      stop 
    endif

   enddo

   print*, 'fail'
   

  
end program gauss_metropolis
