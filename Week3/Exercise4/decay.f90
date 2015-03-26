!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          Simulation of radioactive decay
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
PROGRAM decay
  IMPLICIT none
  REAL :: lambda
  REAL :: r,dt
  INTEGER :: i, t, nleft1, nleft2, start, sizer, dn1 , dn2
  integer, dimension(:), allocatable :: seed
  !
  call random_seed(sizer)
  allocate(seed(sizer))
  print *,'Here the seed has ',sizer,' components; insert them (or print "/") >'  
  read(*,*)seed
  call random_seed(put=seed)

  !         initial values
  print* ,"probability of each atom to decay in dt:" 
  read *, lambda
  print *,"initial number of nuclei >"
  read *, start
  dt = 1
  t = 1          ! initialize time
  nleft1 = start  ! at the beginning N(t=0)=start
  nleft2 = start  ! sottostima
  ! N(t) nuclei left at time t,
  ! that have a given probability lambda of decay
  ! in the time interval t:t+dt
  !            
  !         open output file
  OPEN(unit=7, FILE="decay.dat", status="replace",action="write")
  WRITE (unit=7,fmt=*) "# t ,       N(t)"
  WRITE (unit=7,fmt=*) "0  ", nleft1, nleft2 !REAL(nleft)/start
  !
  !         Execution
  DO                               ! time loop
     dn1 = 0
     dn2 = 0

     DO  i = 1, nleft1             ! loop on the nuclei left
        call random_number(r)
        IF (r < lambda) THEN
           dn1 = dn1 + 1       ! update the number of nuclei left
        ENDIF
     END DO

     i=1
     do while ( i < nleft2 )
        call random_number(r)
        IF (r < lambda) THEN
           dn2 = dn2 + 1       ! update the number of nuclei left
           i = i + 1
        ENDIF
        i = i + 1
     end do
     
     nleft1 = nleft1 - dn1
     !nleft = nleft - lambda * nleft * dt  
     nleft2 = nleft2 - dn2

     WRITE (unit=7,fmt=*) t , nleft1, nleft2 ! or REAL(nleft)/start
     if (nleft1 <= 0) exit
     t = t + 1
  END DO
  !
  close(7)
  stop
END program decay
