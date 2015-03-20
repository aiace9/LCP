cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  fit.f: Least square fit to decay spectrum                           c
c								       c
c  taken from: "Projects in Computational Physics" by Landau and Paez  c 
c	       copyrighted by John Wiley and Sons, New York            c      
c    			                                               c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Program fit
      Implicit none
c
c declarations
      Integer i, j
      Real*8 s, sx, sy, sxx, sxy, delta, inter, slope
      Real*8 x(12), y(12), d(12)
c
c input value y - exponential fit y > 0      
      Data y /32, 17, 21, 7, 8, 6, 5, 2, 2, 0.1, 4, 1/
c
c input values x
      Do 10 i=1, 12
         x(i)=i*10-5
 10   Continue
c
c input value delta y - estimate
      Do 11 i=1, 12
         d(i)=1.0
 11   Continue
c
c take logs of y values for exponential fit
      Do 20 i=1, 12
         y(i)=log(y(i))
 20   Continue
c
c calculate all the sums
      Do 30 i=1, 12
         s   = s   +         1 / (d(i)*d(i))
         sx  = sx  +      x(i) / (d(i)*d(i))
         sy  = sy  +      y(i) / (d(i)*d(i))
         sxx = sxx + x(i)*x(i) / (d(i)*d(i)) 
         sxy = sxy + x(i)*y(i) / (d(i)*d(i))
 30   Continue
c
c calculate the coefficients
      delta= s*sxx-sx*sx
      slope=  (s*sxy-sx*sy) / delta		
      inter=(sxx*sy-sx*sxy) / delta
      Write(*,*) 'intercept=', inter
      Write(*,*) 'slope=', slope
      Write(*,*) 'correlation=', -sx/sqrt(sxx*s)		

      Stop 'fit'
      End



