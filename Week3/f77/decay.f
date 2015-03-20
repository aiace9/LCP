CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c          Simulating radioactive decay
c 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM decay
      IMPLICIT none
c
c          Declarations
      REAL*8 r, DRAND48, lambda
      INTEGER i, h, nleft, nloop, seed, start
c
c          Set parameters (decay constant, starting number of atoms, seed)
      lambda=0.01d0
      start=1000 
      seed=11168
c
c   	  Set initial values, seed random number generator     
      h = 1
      nloop = start
      nleft = start
      call SRAND48(seed)
c            
c	  open output file
      OPEN(7, FILE='decayf.dat')
c
c  	  Execution
      DO 20 WHILE (nleft .NE. 0) 
         DO 10 i = 1, nleft
            r = DRAND48()
            IF (r .LE. lambda) THEN
               nloop = nloop -1
            ENDIF
 10      CONTINUE
c
         nleft = nloop
         WRITE (7,*) h, ' ', REAL(nleft)/start
         h = h + 1
 20   END DO
c
      STOP 
      END
