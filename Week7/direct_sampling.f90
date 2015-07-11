!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! direct_sampling.f90
!
! DIRECT sampling of several physical observables for the
! hamiltonian:         h = -1/2 \nabla^2 + 1/2 x^2),
! comparison exact expected results with numerical results
! on psi^2(x), with  psi(x) = exp(-x^2/(4 sigma^2)
! (sigma=1  => psi^2(x) = costant * standard gaussian
!  P(x) =  exp(-x**2/(2*sigma**2))/sqrt(2*pi*sigma**2)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Module gaussian
  implicit none
  public :: gasdev

contains

  SUBROUTINE gasdev(x)

    REAL, INTENT(OUT) :: x
    REAL :: rsq, v1, v2
    REAL, SAVE :: g
    LOGICAL, SAVE :: gaus_stored=.false.

    if (gaus_stored) then
       x=g
       gaus_stored=.false.
    else
       do
          call random_number(v1)
          call random_number(v2)
          v1 = 2.0*v1 - 1.0
          v2 = 2.0*v2 - 1.0
          rsq = v1**2 + v2**2
          if (rsq > 0.0 .and. rsq < 1.0) exit
       end do
       rsq  = sqrt(- 2.0*log(rsq)/rsq)
       x = v1*rsq
       g = v2*rsq
       gaus_stored = .true.
    end if

  END SUBROUTINE gasdev

end module gaussian

program direct_sampling
  use gaussian
  implicit none
  integer, parameter :: dp=selected_real_kind(13)
  integer :: i,n
  integer, dimension(12) :: seed
  real :: rnd
  real(kind=dp):: sigma,etot,ekin,epot
  real(kind=dp):: x,x1,x2,x3,x4
  character(len=13), save :: format1 = "(a7,2x,2f9.5)"

  integer :: ios

  x1 = 0.0_dp
  x2 = 0.0_dp
  x3 = 0.0_dp
  x4 = 0.0_dp
  ekin = 0.0_dp
  epot = 0.0_dp

  print*, "seed, n, sigma ="
  read*,   seed
  read*,n,sigma
  call random_seed(put=seed)

  open(unit=1, file='plot.dat', iostat=ios, status="unknown", action="write")
  if ( ios /= 0 ) stop "Error opening file plot.dat"
  
  do i=1,n
     !cccccccccccccccccccccccccc
     call gasdev(rnd)   !
     x=rnd*sigma        ! direct sampling
     !ccccccccccccccccccccccccc!

     !harmonic ocillator
     ! H= -1/2 grad**2 + 1/2 x**2
     ekin = ekin - 0.5_dp * ((x/(2*sigma**2))**2 - 1/(2*sigma**2))
     epot = epot + 0.5_dp * x**2
     !look at the slide for more information slide
     !I should write some extra note for this

     etot = ekin + epot
     x1 = x1 + x
     x2 = x2 + x**2
     x3 = x3 + x**3
     x4 = x4 + x**4
     write(unit=1, fmt=*, iostat=ios) i, epot/i - 5.d-1 * sigma**2, ekin/i - 1 / (8.d0 * sigma**2)&
      , x2/i - sigma ** 2
     if ( ios /= 0 ) stop "Write error in file unit 1"
     
  end do

  close(unit=1, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 1"
  

  write(unit=*,fmt=*)"Results (simulation vs. exact result):"
  write(unit=*,fmt=format1)"etot = ",etot/n,1.0_dp/(8.0_dp*sigma**2)&
       +0.5_dp*sigma**2
  write(unit=*,fmt=format1)"ekin = ",ekin/n,1.0_dp/(8.0_dp*sigma**2)
  write(unit=*,fmt=format1)"epot = ",epot/n,0.5_dp*sigma**2
  write(unit=*,fmt=format1)"<x>  = ",x1/n,0.0_dp
  write(unit=*,fmt=format1)"<x^2>= ",x2/n,sigma**2
  write(unit=*,fmt=format1)"<x^3>= ",x3/n,0.0_dp
  write(unit=*,fmt=format1)"<x^4>= ",x4/n,3.0_dp*sigma**4

end program direct_sampling








