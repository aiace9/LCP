!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! metropolis_gaussian.f90    
!
! METROPOLIS sampling of several physical observables for the
! hamiltonian:         h = -1/2 \nabla^2 + 1/2 x^2),
! comparison exact expected results with numerical results
! on psi^2(x), with  psi(x) = exp(-x^2/(4\sigma^2))
! \sigma=1 => psi^2(x) = costant * standard gaussian
!  P(x) =  exp(-x**2/(2*sigma**2))/sqrt(2*pi*sigma**2)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program metropolis_gaussian
  implicit none
  integer, parameter :: dp=selected_real_kind(13)
  integer :: i,n
!  integer, dimension(4) :: seed
  real(kind=dp):: sigma,etot,ekin,epot,rnd,ekinL,epotL,etot2
  real (kind=dp) :: pigr,pi2b,var,beta,var_th
  real(kind=dp):: x,x1,x2,xp,delta,expx,expxp,p,acc
  character(len=13), save :: format1 = "(a7,2x,2f9.5)"
  open(unit=7,file='e_var_gauss.dat',position='append')

  pigr = 2*asin(1.0_dp)
  acc = 0.0_dp
  x1 = 0.0_dp
  x2 = 0.0_dp
  ekin = 0.0_dp
  epot = 0.0_dp
  etot2= 0.0_dp
  print*, "n, sigma  (remember: beta = 1 / (4*sigma**2)), x0, delta"
  read*,   n,sigma,x,delta
  beta = 1 / (4*sigma**2)
!  call random_seed(put=seed)

  do i=1,n
     ekinL = - 0.5_dp * ((x/(2*sigma**2))**2 - 1/(2*sigma**2))
     epotL = 0.5_dp * x**2 + 0.125_dp * x**4
     ekin = ekin + ekinL
     epot = epot + epotL
     etot = ekin + epot
     etot2 = etot2 + (ekinL + epotL)**2
     x1 = x1 + x
     x2 = x2 + x**2
     !ccccccccccccccccccccccccccccccc
     expx = - x**2 /(2*sigma**2)    !
     call random_number(rnd)        !
     xp = x + delta * (rnd-0.5_dp)  !
     expxp = - xp**2 /(2*sigma**2)  !   metropolis
     p = exp (expxp-expx)           !   algorithm
     call random_number(rnd)        !
     if (p > rnd) then              !
        x = xp                      !
     !ccccccccccccccccccccccccccccccc
        acc=acc+1.0_dp                
     endif
  enddo

  var_th =  1._dp/(32*beta**2)+beta**2/2-1._dp/4

  write(unit=*,fmt=*)"acceptance ratio = ",acc/n
  write(unit=*,fmt=*)"# Results (simulation vs. exact results - not valid, check notes.):"
  write(unit=*,fmt=format1)"etot = ",etot/n,(1/2.0_dp-2.0_dp*beta**2) * 1 / (4.0_dp*beta) +beta + 3 /(128.0_dp * beta**2)
  write(unit=*,fmt=format1)"ekin = ",ekin/n,1.0_dp/(8.0_dp*sigma**2)
  write(unit=*,fmt=format1)"epot = ",epot/n,0.5_dp*sigma**2 + (3.0_dp/16.0_dp) * (1.0_dp/(8.0_dp*beta**2))
  write(unit=*,fmt=format1)"evar = ",etot2/n-(etot/n)**2,var_th
  write(unit=*,fmt=format1)"<x>  = ",x1/n,0.0_dp
  write(unit=*,fmt=format1)"<x^2>= ",x2/n,sigma**2

  write(7,*)sigma,etot/n,sqrt(abs(etot2/n-(etot/n)**2)),etot2/n-(etot/n)**2

 close(7)

end program metropolis_gaussian

