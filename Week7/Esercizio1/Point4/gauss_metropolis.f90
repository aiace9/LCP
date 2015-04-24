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

  n_min = 100
  n_max = 10000
  n_step= 100
  
  open(unit=1, file='variance.dat', iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file variance"

 
  do j=n_min,n_max,n_step
    n = j
    acc = 0.0_dp
    x = x0
    x1 = 0.0_dp
    x2 = 0.0_dp
    x3 = 0.0_dp
    x4 = 0.0_dp
    do i=1,n
       x1 = x1 + x
       x2 = x2 + x**2
       x3 = x3 + x**3
       x4 = x4 + x**4
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
    enddo
    write(unit=*,fmt=*)"---------------"
    write(unit=*,fmt=*)"# n, x0, delta = ",n,x0,delta
    write(unit=*,fmt=*)"# Results (simulation vs. exact results):"
    write(unit=*,fmt=*)"# acc = ",acc/n
    write(unit=*,fmt=format1)"# var2 = ",x2/n-(x1/n)**2,sigma**2

    write(unit=1, fmt=*, iostat=ios) n, x2/n-(x1/n)**2 - sigma**2
    if ( ios /= 0 ) stop "Write error in file unit 1"
  
  enddo

  close(unit=1, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 1"
  
end program gauss_metropolis
