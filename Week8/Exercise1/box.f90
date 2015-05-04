!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! box.f90
!
! simulation of the evolution of a physical system towards equilibrium:
! non interacting particles in a box divided into two parts; 
! at each time step,  one and only one particle (randomly choosen) 
! goes from one side to the other one
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module moduli_box
  use istogram_module
  implicit none
  public :: initial, move
  integer, public :: N,tmax

contains
  
  subroutine initial()
    integer , dimension(12) :: seed  ! on infis: 4 with g95, 8 with gfortran
    print*," total number of particles N >"
    read*,N
    tmax = 10*N  ! we choose the evolution time proportional to N 
    print*," seed (1:4) >"
    read *, seed
    call random_seed(put=seed)
  end subroutine initial

  subroutine move()
    integer :: nl,itime
    real :: r, prob
    
    !my additions
    integer :: supp,n1,n2
    integer :: i
    !initial reset
    n1 = 0
    n2 =0
    !starting position
    !nl = N ! we start with all the particles on the left side
    nl = N/2 ! impose the equilibrium as a starting condition

    open(unit=2,file="box.dat",action="write",status="replace")

    do itime = 1,tmax
       prob = real(nl)/N    ! fraction of particles on the left
                            ! the probability of jumping is 1/N for each particle why? 
                            ! multiplied by the number of particle on that side.
       supp = 0
       do i = 1, 5, 1
        call random_number(r)  
          if (r <= prob) then
              supp = supp - 1
          else
              supp = supp + 1
         end if
       end do
       
       nl = nl + nint(supp/real(5))
       write(unit=2,fmt=*)nl
       call push_ist(dble(nl))
       n1 = n1 + nl
       n2 = n2 + nl**2
    end do
    print*, "sigma^2 = ", n2/dble(tmax) - (n1/dble(tmax))**2
    print*, "<n> = ", n1/dble(tmax)
    print*, "sigma/<n> = ", sqrt(n2/dble(tmax) - (n1/dble(tmax))**2)/n1/dble(tmax)
    close(unit=2)
  end subroutine move
end module moduli_box

program box
  
  use moduli_box
  use istogram_module

  implicit none
  
  !integer, parameter :: dp = selected_real_kind(14,200)

  ! compare a random number with the fraction of particles on the left, nl/N:
  ! if r.le.nl/N we move one particle from tyhe left to the right;
  ! elsewhere from the right to the left
  call initial()
  call init_ist(0d0,dble(N),1d0)
  call move()
  call save_data_ist()
end program box


