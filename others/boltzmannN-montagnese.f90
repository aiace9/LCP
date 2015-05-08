!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! boltzmann.f90
! algoritmo di Metropolis usato come algoritmo di importance sampling:
! per generare microstati con probabilita' secondo la distribuzione di
! Boltzmann, qui per 1 particella classica 1D
! La quantita' di interesse e' la probabilita' P(E)dE che una particella
! abbia energia tra E e E+dE (qui E labella un microstato, a meno del segno)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module common
  implicit none
  public :: initial, Metropolis, data, probability, averages
  real, public :: E,T,del_E,beta,dvmax,accept
  integer, public, dimension(1) :: seed
  integer, public :: nbin,nmcs, nequil, npart
  real, public, dimension(:), allocatable :: P, vel

contains

  subroutine initial(nequil,vcum,ecum,e2cum)
    real, intent(out) :: vcum, ecum,e2cum
    integer, intent(out) :: nequil
    print*," numero di particelle >"
    read *, npart
    print*," numero di steps MC >"
    read *, nmcs
    print*," temperatura assoluta >"
    read *, T
    print*," velocita` iniziale >"
    read *, vel
    print*," massima variaz. nella vel. >"
    read *, dvmax
    print*," seme >"
    read *, seed(1)
    call random_seed(put=seed)
    allocate (vel(npart))

    beta   = 1/T
    nequil = 0.7 * nmcs   ! VERIFICARE LA SCELTA
    E      = 0.5 * sum( vel * vel)
    del_E  = 0.01        ! intervallo di energia per l"istogr. P
    nbin   = int(4*T / del_E)  ! max. numero di canali 
    print *,"# T       :",T
    print *,"# <E0>    :",E
    print *,"# <v0>    :",vel
    print *,"# dvmax   :",dvmax
    print *,"# nMCsteps:",nmcs
    print *,"# deltaE  :",del_E
    print *,"# nbin    :",nbin
    open(unit=9,file="boltzmann.dat",status="replace",action="write")
    open(unit=10,file="temp.dat",position='append',action="write")
    write(unit=9,fmt=*)"# npart   :",npart      
    write(unit=9,fmt=*)"# T       :",T
    write(unit=9,fmt=*)"# <E0>    :",E
    write(unit=9,fmt=*)"# <v0>    :",vel
    write(unit=9,fmt=*)"# dvmax   :",dvmax
    write(unit=9,fmt=*)"# nMCsteps:",nmcs
    write(unit=9,fmt=*)"# deltaE  :",del_E
    write(unit=9,fmt=*)"# nbin    :",nbin
    allocate (P(0:nbin))
    ecum  = 0.0
    e2cum = 0.0
    vcum  = 0.0
    P     = 0.0
    accept= 0.0
  end subroutine initial

  subroutine Metropolis() 
    real, dimension(npart) :: dv,vtrial,de,rnd
    real :: rnd1
    call random_number(rnd)
    dv = (2*rnd - 1) * dvmax              ! variaz. di prova in v
    vtrial = dv + vel                    ! nuova vel. (di prova)
    de = 0.5 * (vtrial*vtrial - vel*vel)  ! variaz. in E
    call random_number(rnd1)
    if (sum(de) >= 0.0) then
       if ( exp(-beta*sum(de)) < rnd1 ) return ! mossa non accettata
    end if
    vel = vtrial
    accept = accept + 1 
    E = E + sum(de)
  end subroutine Metropolis


  subroutine data(vcum,ecum,e2cum)
    real, intent(inout) :: vcum,ecum,e2cum
    Ecum = Ecum + E
    E2cum = E2cum + E*E
    vcum = vcum + sum(vel)
    call probability() 
  end subroutine data


  subroutine probability() 
    integer :: ibin
    ibin = nint(E/(npart*del_E))
    if ( ibin <= nbin )     P(ibin) = P(ibin) + 1
  end subroutine probability


  subroutine averages(vcum,Ecum,E2cum) 
    real, intent(in) :: vcum,Ecum,E2cum
    real :: znorm, Eave, E2ave, vave, sigma2, logP
    integer :: ibin
    znorm  = 1.0/(nmcs)
    accept = accept / (nmcs + nequil) ! acceptance ratio
    Eave   = Ecum * znorm / npart   ! en. media
    E2ave  = E2cum * znorm / npart  ! 
    vave   = vcum * znorm / npart  ! vel. media
    sigma2 = E2ave - Eave*Eave
    print *,"# <E>     :",Eave
    print *,"# <v>     :",vave
    print *,"# accept. :",accept
    print *,"# sigma   :",sqrt(sigma2)
    write(unit=9,fmt=*)"# <E>     :",Eave
    write(unit=9,fmt=*)"# <v>     :",vave
    write(unit=9,fmt=*)"# accept. :",accept
    write(unit=9,fmt=*)"# sigma   :",sqrt(sigma2)
    write(unit=9,fmt=*)"# bin, P(E), log(P(E))"
    do ibin = 0,nbin
 	logP = 0.0
        if ( P(ibin) /= 0.0) logP = log (P(ibin) * znorm)
        write(unit=9,fmt=*) ibin*del_E, P(ibin) * znorm, logP
    end do
    write(unit=10,fmt=*) T, Eave, E2ave
    close(unit=9)
    close(unit=10)
  end subroutine averages
end module common

program Boltzmann
  use common
  real :: vcum, ecum, e2cum
  integer :: imcs
  ! stabilisce i parametri e inizializza le variabili
  call initial(nequil,vcum,ecum,e2cum)
  do  imcs = 1 , nmcs + nequil
     call Metropolis()
     ! accumula dati dopo ogni mossa
     if ( imcs > nequil ) call data(vcum,ecum,e2cum) 
  end do
  call averages(vcum,Ecum,E2cum)
  print*, nequil+nmcs
end program Boltzmann

