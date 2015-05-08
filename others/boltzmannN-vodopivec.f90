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
  real, public :: T,del_E,beta,dvmax,accept
  integer, public, dimension(1) :: seed
  integer, public :: N,nbin,nmcs
  real, public, dimension(:), allocatable :: P
  real, public, dimension(:), allocatable :: vel,E
contains

  subroutine initial(nequil,vcum,ecum,e2cum)
    real, intent(out) :: vcum,ecum,e2cum
    integer, intent(out) :: nequil
    real :: v0
    print*," numero particelle >"
    read *, N
    print*," numero di steps MC >"
    read *, nmcs
    print*," temperatura assoluta >"
    read *, T
    print*," velocita` iniziale >"
    read *, v0
    print*," massima variaz. nella vel. >"
    read *, dvmax
    print*," seme >"
    read *, seed(1)
    call random_seed(put=seed)
    beta   = 1/T
    nequil = 0.1 * nmcs   ! VERIFICARE LA SCELTA
    del_E  = 0.25         ! intervallo di energia per l"istogr. P
    nbin   = int(4*T / del_E)  ! max. numero di canali 
    allocate (vel(N))
    allocate (E(N))
    vel   = v0
    E     = 0.5 * vel * vel
    print *,"# T       :",T
    print *,"# <E0>    :",E(1)
    print *,"# <v0>    :",v0
    print *,"# dvmax   :",dvmax
    print *,"# nMCsteps:",nmcs
    print *,"# deltaE  :",del_E
    print *,"# nbin    :",nbin
    open(unit=9,file="boltzmann.dat",status="replace",action="write")
    write(unit=9,fmt=*)"# T       :",T
    write(unit=9,fmt=*)"# <E0>    :",E(1)
    write(unit=9,fmt=*)"# <v0>    :",v0
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
  
  subroutine deinit
    deallocate(E)
    deallocate(vel)
    deallocate(P)
  end subroutine

  subroutine Metropolis() 
    real :: dv,vtrial,de,rnd
    integer :: i
    
    do i=1,N
	call random_number(rnd)
	dv = (2*rnd - 1) * dvmax              ! variaz. di prova in v
	vtrial = vel(i) + dv                  ! nuova vel. (di prova)
	de = 0.5 * (vtrial*vtrial - vel(i)*vel(i))  ! variaz. in E
	call random_number(rnd)
	if (de >= 0.0 .and. exp(-beta*de) < rnd) then ! mossa non accettata
		cycle
	end if
	vel(i) = vtrial
	accept = accept + 1 
	E(i) = E(i) + de
    end do
  end subroutine Metropolis


  subroutine data(vcum,ecum,e2cum)
    real, intent(inout) :: vcum,ecum,e2cum
    integer :: i
    do i=1,N
    	Ecum = Ecum + E(i)
    	E2cum = E2cum + E(i)*E(i)
    	vcum = vcum + vel(i)
    end do
    call probability() 
  end subroutine data


  subroutine probability() 
    integer :: ibin,i
    real :: Etot
    Etot=0.0
    do i=1,N
    	Etot=Etot+E(i)
    end do
    ibin = int(Etot/(del_E*N))
    if ( ibin <= nbin )     P(ibin) = P(ibin) + 1
  end subroutine probability

  
  subroutine averages(nequil,vcum,Ecum,E2cum) 
    real, intent(in) :: vcum,Ecum,E2cum
    integer, intent(in) :: nequil
    real :: znorm, Eave, E2ave, vave, sigma2, logP
    integer :: ibin
    znorm  = 1.0/nmcs/N
    accept = accept /(nmcs+nequil)/N ! acceptance ratio
    Eave   = Ecum * znorm   ! en. media
    E2ave  = E2cum * znorm  ! 
    vave   = vcum * znorm   ! vel. media
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
!         write(unit=9,fmt=*) ibin,P(ibin) * znorm,logP
        write(unit=9,fmt=*) ibin*del_E,P(ibin)/nmcs
    end do
    close(unit=9)
  end subroutine averages
end module common

program Boltzmann
  use common
  real :: vave, vcum, ecum, e2cum
  real :: Eave, E2ave, sigma2
  integer :: imcs,nequil
  ! stabilisce i parametri e inizializza le variabili
  call initial(nequil,vcum,ecum,e2cum)
  do  imcs = 1 , nmcs + nequil
     call Metropolis()
     ! accumula dati dopo ogni mossa
     if ( imcs > nequil ) then
     	call data(vcum,ecum,e2cum) 
      end if
  end do
  call averages(nequil,vcum,Ecum,E2cum)
  call deinit()
end program Boltzmann

! 1f
! 20, 200, 100, 10, 30, -437
! dvmax=30 c.a. accept ratio 0.5 (provato con diversi seed)
! # <E>     :   46.44997
! # <v>     :  0.1143214

! 20, 100000, 100, 10, 3, -437
! # <E>     :   49.92599
! # <v>     : -2.7030865E-02
!
! 1g 
!  la distribuzione P(E) tende a una gaussiana all'aumentare del numero di particelle:
!  somma di N distribuzioni exp(-x)
! 20, 100000, 100, 10, 30, -437
! 1000, 100000, 100, 10, 30, -437




 
