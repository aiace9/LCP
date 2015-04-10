!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       pi.f90: Calcola pigreco usando l'integrazione Monte Carlo
!C       ("lanciando sassi" in un cerchio)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Program pi

  Implicit none
  integer, dimension(12) :: seed
  real, dimension(2) :: rnd
  Real :: area, x, y
  Integer :: i, max, pigr
    print*,' max numero di punti='
  read*, max
  print*,' seed (1:4) >'
  read*, seed
  call random_seed(put=seed)
  !       apre il file per salvare i dati, inizializza pigr e il generatore
  Open(7, File='pigr.dat', Status='Replace')
  pigr=0
  ! esegue generando coppie di numeri a caso tra -1 e 1
  ! cioe' lanciando sassi nel quadrato di lato 2 e centrato in (0,0):
  ! i sassi che cadono nel cerchio (cioe' con x*x+y*y < o = 1) li conteggiamo 
  ! ai fini del computo di pigr, che e' 4*areacerchio/areaquadrato
  ! cioe' 4*numero sassi caduti nel cerchio/numero sassi lanciati nel quadrato
  Do i=1, max
     call random_number(rnd)
     x = rnd(1)*2-1
     y = rnd(2)*2-1
     If ((x*x + y*y) <= 1) then
        pigr = pigr+1
     Endif
     area = 4.0 * pigr/Real(i)
     if(mod(i,1000)==0)Write(7,*) i, abs(acos(-1.)-area)
  end do
  Close(7)
  Stop 'dati salvati in pigr.dat '
End program pi



















