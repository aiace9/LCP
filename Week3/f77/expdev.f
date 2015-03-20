      function expdev(idum)
      IMPLICIT none
      INTEGER idum
c idum qui e' un intero positivo che e' il seed per la ran
      REAL expdev, dum, ran
c
 1         dum = ran(idum)
      if(dum.eq.0) go to 1
      expdev = -log(dum)
      STOP 
      END
