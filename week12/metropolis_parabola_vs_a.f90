!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! metropolis_sampling.f90
!
! METROPOLIS sampling of several physical observables for the
! hamiltonian:         h = -1/2 \nabla^2 + 1/2 x^2),
! comparison exact expected results with numerical results
! on psi^2(x), with  psi(x) = B(a^2-x^2) for |x|<a; 0 elsewhere
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program metropolis_sampling
  implicit none
  integer, parameter :: dp=selected_real_kind(13)
  integer :: i,n, j
!  integer, dimension(2) :: seed
  real(kind=dp):: a,amin,amax,astep,etot,ekin,epot,rnd,ekinL,epotL,etot2,var
  real(kind=dp):: x,x1,x2,xp,delta,psi,psip,p,acc,var_th,etot_th,ekin_th,epot_th
  character(len=13), save :: format1 = "(a7,2x,1f9.5)"
  character(len=17), save :: format2 = "(a17,2(2x,1f9.5))"
  open(unit=1,file='etot_var_vs_a',position='append')
  write(1,*)"# a,  etot_num, sqrt(var)_num, etot_th, sqrt(var)_th"

 10   print*, "amin, amax, astep, nmcs, x0"
  read*, amin, amax, astep,  n, x


do j = 1,int((amax-amin)/astep)

 acc = 0.0_dp
  x1 = 0.0_dp
  x2 = 0.0_dp
  ekin = 0.0_dp
  epot = 0.0_dp
  etot2= 0.0_dp

a = amin+(j-1)*astep
delta = 3. * a
  if (abs(x) >= a) go to 10
!  call random_seed(put=seed)

  do i=1,n
     ekinL = 1/(a**2-x**2)
     epotL = 0.5_dp * x**2
     ekin = ekin + ekinL
     epot = epot + epotL
     etot = ekin + epot
     etot2 = etot2 + (ekinL + epotL)**2
     x1 = x1 + x
     x2 = x2 + x**2
     !ccccccccccccccccccccccccccccccc
     psi = a**2-x**2
     call random_number(rnd)        !
     xp = x + delta * (rnd-0.5_dp)  !
     if (abs(xp) >= a) then         !
        psip = 0.0_dp !             !
        else                        !
           psip = a**2-xp**2        !   metropolis
        end if                      !
     p = (psip/psi)**2              !   algorithm
     call random_number(rnd)        !
     if (p > rnd) then              !
        x = xp                      !
     !ccccccccccccccccccccccccccccccc
        acc=acc+1.0_dp                
     endif
  enddo

  ekin_th = 5./4/a**2
  epot_th = a**2/14
  etot_th = ekin_th + epot_th
  var_th = (15./8/a**5) * ( a + 2*a**9/(7*45) + 2*a**5/15) - &
           (5./4/a**2+a**2/14)**2

  write(unit=*,fmt=*)"acceptance ratio = ",acc/n
  write(unit=*,fmt=*)"# Results (simulation vs. exact results):"
  write(unit=*,fmt=format2)"etot (num./th.)= ",etot/n, etot_th
  write(unit=*,fmt=format2)"ekin (num./th.)= ",ekin/n, ekin_th
  write(unit=*,fmt=format2)"epot (num./th.)= ",epot/n, epot_th
  write(unit=*,fmt=format2)"evar (num./th.)= ",etot2/n-(etot/n)**2,var_th
  write(unit=*,fmt=format1)"<x>  = ",x1/n
  write(unit=*,fmt=format1)"<x^2>= ",x2/n

  write(1,*)a,etot/n,sqrt(etot2/n-(etot/n)**2),etot_th,sqrt(var_th)

end do
  close(1)
end program metropolis_sampling
