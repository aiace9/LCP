      program Boltzmann
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c algoritmo di Metropolis usato come algoritmo di importance sampling:
c per generare microstati con probabilita' secondo la distribuzione di
c Boltzmann, qui per 1 particella classica 1D
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      parameter(maxnbin=500)
      dimension P(0:maxnbin)
c stabilisce i parametri e inizializza le variabili
      open(unit=9,file='boltzmann.dat',status='unknown')
      call initial(nmcs,nequil,beta,vel,E,dvmax,nbin,del_E,idum)
      if (nbin.gt.maxnbin) call error(nbin,maxnbin)
      accept = 0
      do 20 ibin = 0,maxnbin
 20      P(ibin) = 0.
      do 10 imcs = 1 , nmcs + nequil
         call Metropolis(beta,vel,E,dvmax,accept,idum)
c accumula dati dopo ogni mossa
         if(imcs.gt.nequil)call data(E,vel,Ecum,E2cum,vcum,nbin,P,del_E)
 10      continue
      call averages(nmcs,vcum,Ecum,E2cum,nbin,P,accept,beta)
      close(9)
      stop
      end


      subroutine initial(nmcs,nequil,beta,vel,E,dvmax,nbin,del_E,idum)
      print*,' numero di steps MC >'
      read(*,*)nmcs
      print*,' temperatura assoluta >'
      read(*,*)T
      print*,' velocita` iniziale >'
      read(*,*)vel
      print*,' massima variaz. nella vel. >'
      read(*,*)dvmax
      print*,' n. per inizializzare la seq. random (int.<0)>'
      read(*,*)idum
      r = ran2(idum) !n. random di prova
      beta   = 1/T
      nequil = 0.1 * nmcs   ! VERIFICARE LA SCELTA
      E      = 0.5 * vel * vel
      del_E  = 0.05         ! intervallo di energia per l'istogr. P
      nbin   = 4*T / del_E  ! max. numero di canali 
      write(6,*)'# T       :',T
      write(6,*)'# <E0>    :',E
      write(6,*)'# <v0>    :',vel
      write(6,*)'# dvmax   :',dvmax
      write(6,*)'# nMCsteps:',nmcs
      write(6,*)'# deltaE  :',del_E
      write(6,*)'# nbin    :',nbin
      write(9,*)'# T       :',T
      write(9,*)'# <E0>    :',E
      write(9,*)'# <v0>    :',vel
      write(9,*)'# dvmax   :',dvmax
      write(9,*)'# nMCsteps:',nmcs
      write(9,*)'# deltaE  :',del_E
      write(9,*)'# nbin    :',nbin
      return
      end


      subroutine error(n,m)
      print*,' numero di nbin troppo alto! (',n,'.ge.',m,')'
      stop
      end


      subroutine Metropolis(beta,vel,E,dvmax,accept,idum)
      dv = (2*ran2(idum) - 1) * dvmax       ! variaz. di prova in v
      vtrial = vel + dv                     ! nuova vel. (di prova)
      de = 0.5 * (vtrial*vtrial - vel*vel)  ! variaz. in E
      if (de .gt. 0.) then
         if ( exp(-beta*de) .lt. ran2(idum) ) return ! mossa non accettata
      end if
      vel = vtrial
      accept = accept + 1 
      E = E + de
      return
      end


      subroutine data(E,vel,Ecum,E2cum,vcum,nbin,P,del_E)
      dimension P(0:nbin)
      Ecum = Ecum + E
      E2cum = E2cum + E*E
      vcum = vcum + vel
      call probability(E,nbin,P,del_E)
      return
      end


      subroutine probability(E,nbin,P,del_E)
      dimension P(0:nbin)
      ibin = E/del_E
      if ( ibin .gt. nbin ) ibin = nbin
      P(ibin) = P(ibin) + 1
      return
      end


      subroutine averages(nmcs,vcum,Ecum,E2cum,nbin,P,accept,beta)
      dimension P(0:nbin)
      znorm  = 1./nmcs
      accept = accept * znorm ! acceptance ratio
      Eave   = Ecum * znorm   ! en. media
      E2ave  = E2cum * znorm  ! 
      vave   = vcum * znorm   ! vel. media
      sigma2 = E2ave - Eave*Eave
      write(*,*)'# <E>     :',Eave
      write(*,*)'# <v>     :',vave
      write(*,*)'# accept. :',accept
      write(*,*)'# sigma   :',sqrt(sigma2)
      write(9,*)'# <E>     :',Eave
      write(9,*)'# <v>     :',vave
      write(9,*)'# accept. :',accept
      write(9,*)'# sigma   :',sqrt(sigma2)
      write(9,*)'# bin, P(E), log(P(E))'
      do 10 ibin = 0,nbin
         if ( P(ibin) .gt. 0.01*znorm ) then
            prob = P(ibin) * znorm
c            write(*,*)ibin,prob,log(prob)
            write(9,*)ibin,prob,log(prob)
            end if
 10         continue
      return
      end

      function ran2(idum)
c distribuzione uniforme in (0,1) (NO gli endpoints)
c Generatore di n. pseudorandom con un lungo periodo (> 2*10^18)
c metodo di L'Ecuyer + altre sottigliezze
c Inizializzare con IDUM = un intero negativo (poi NON cambiarlo!)
c RNMX approssima il piu' grande floating number minore di 1/
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
      if(idum.le.0)then                   ! inizializza
       idum=max(-idum,1)                  ! si vuole evitare idum=0
       idum2=idum                         
       do j=ntab+8,1,-1                   ! comincia con lo shuffle
        k=idum/iq1                          ! dopo 8 mosse di riscaldamento
        idum=ia1*(idum-k*iq1)-k*ir1
        if(idum.lt.0)idum=idum+im1
        if(j.le.ntab)iv(j)=idum
       enddo
       iy=iv(1)
      endif
      k=idum/iq1                          ! se non inizializza, parte da qui
      idum=ia1*(idum-k*iq1)-k*ir1         ! Calcola idum=mod(IA1*idum,IM1)
      if(idum.lt.0)idum=idum+im1            ! in modo "furbo" per evitare
      k=idum2/iq2                           ! overflow (metodo di Schrage)
      idum2=ia2*(idum2-k*iq2)-k*ir2       ! Calcola idum2=mod(IA2*idum2,IM2)
      if(idum2.lt.0)idum2=idum2+im2
      j=1+iy/ndiv                         ! j e' nell'intervall 1:NTAB
      iy=iv(j)-idum2                      ! idum e' shuffled; idum1,2 combinati
      iv(j)=idum                            ! per generare l'output
      if(iy.lt.1)iy=iy+imm1
      ran2=min(am*iy,rnmx)                ! NO punti estremi
      return
      end



