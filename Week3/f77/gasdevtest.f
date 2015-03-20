      program gasdevtest
      implicit none
      integer idum,i
      real*8 ran2,gasdev
      write(*,*)'idum (<0)'      
      read(*,*)idum
      do i=1,10
      write(*,*)gasdev(idum)
      end do
      stop
      end
         

      function ran2(idum)
c distribuzione uniforme in (0,1)
      implicit none
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real*8 ran2,am,eps,rnmx
      parameter(im1=2147483563,im2=2147483399,am=1.d0/im1,imm1=im1-1
     &         ,ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211
     &         ,ir2=3791,ntab=32,ndiv=1+imm1/ntab
     &         ,eps=5.6d-17,rnmx=1.d0-eps)
      integer idum2,j,k,iv(ntab),iy
      save iv,iy,idum2
      data idum2/123456789/,iv/ntab*0/,iy/0/
      if(idum.le.0)then
       idum=max(-idum,1)
       idum2=idum
       do j=ntab+8,1,-1
        k=idum/iq1
        idum=ia1*(idum-k*iq1)-k*ir1
        if(idum.lt.0)idum=idum+im1
        if(j.le.ntab)iv(j)=idum
       enddo
       iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0)idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0)idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+imm1
      ran2=min(am*iy,rnmx)
      return
      end

      function gasdev(idum)
c distribuzione gaussiana varianza 1
      implicit none
      integer idum,iset
      real*8 gasdev,fac,gset,rsq,v1,v2,ran2
      save iset,gset
      data iset/0/
      if(iset.eq.0)then
    1  v1=2.d0*ran2(idum)-1.d0
       v2=2.d0*ran2(idum)-1.d0
       rsq=v1**2+v2**2
       if(rsq.ge.1.d0.or.rsq.eq.0.d0)goto 1
       fac=sqrt(-2.d0*log(rsq)/rsq)
       gset=v1*fac
       gasdev=v2*fac
       iset=1
      else
       gasdev=gset
       iset=0
      endif
      return
      end
