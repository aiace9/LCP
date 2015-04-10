!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     int.f90: 
!     integra f(x)=exp(x) nell'intervallo [vmin,vamx]=[0,1] 
!     usando la regola trapezoidale o di Simpson
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module intmod
  public :: f, trapez, simpson
contains

  ! la funzione che si vuol integrare
  !
  function f(x)
    implicit none
    real :: f
    real, intent(in) :: x
    f = exp(x)
    return
  end function f


  ! regola trapezoidale
  !
  function trapez(i, min, max)
    implicit none
    real :: trapez
    integer, intent(in) :: i
    real, intent(in) :: min, max
    integer :: n
    real :: x, interval
    trapez = 0.
    interval= ((max-min) / (i-1))
    ! somma i punti medi
    do n = 2, i-1
       x = interval * (n-1)
       trapez = trapez + f(x) * interval
    end do
    ! aggiunge i punti estremi
    trapez = trapez + 0.5 * (f(min)+f(max)) * interval
    return
  end function trapez


  ! regola di Simpson
  !
  function simpson(i, min, max)
    implicit none
    real :: simpson
    integer, intent(in) :: i
    real, intent(in) :: min, max
    integer :: n
    real :: x, interval
    simpson = 0.
    interval = ((max-min) / (i-1))
    ! ciclo sui punti PARI
    do n = 2, i-1, 2
       x = interval * (n-1)
       simpson = simpson + 4*f(x)
    end do
    ! ciclo sui punti DISPARI
    do n = 3, i-1, 2
       x = interval * (n-1)
       simpson = simpson + 2*f(x)
    end do
    ! aggiunge i punti estremi
    simpson = simpson + f(min)+f(max)
    simpson = simpson * interval/3
    return
  end function simpson

end module intmod

program int
  use intmod
  !
  ! dichiarazioni di variabili:
  !     limite di accuratezza,
  !     minimo e massimo in x
  !
  implicit none
  real :: r1, r2, theo, vmin, vmax
  integer :: i, n
  ! valore teorico esatto
  vmin = 0.0
  vmax = 1.0
  theo = exp(vmax)-exp(vmin)
  print*,' valore teorico esatto =',theo
  open(unit=7,file='int-tra-sim.dat',status='unknown')
  !
  write(7,*)"# N, intervallino, esatto, Trapezi-esatto, Simpson-esatto"
  do i = 2,8
     n = 2**i
     r1 = trapez(n+1, vmin, vmax)
     r1 = (r1-theo)
     r2 = simpson(n+1, vmin, vmax)
     r2 = (r2-theo)
     write(7,'(i4,4(2x,f10.6))') n, 1./n, theo, r1, r2
  end do
  close(7)

  print*,' dati salvati in int-tra-sim.dat (|diff dal val. esatto|)'
  stop

end program int
