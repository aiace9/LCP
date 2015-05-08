      program ising2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c algoritmo di Metropolis per calcolare <E>, <M>, nell'insieme canonico
c (fisso T,N,V) con un modello di Ising 2D
c
c Qui e' posto: K_B = 1
c               J   = 1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c il file 'ising2.dat' e' scritto nel formato opprtuno per usare GNUPLOT:
c vedasi i seguenti comandi:
c 
c gnuplot> set view 0,0
c gnuplot> set nosurface
c gnuplot> set contour
c gnuplot> set cntrparam levels discrete -1,1
c gnuplot> splot "ising2.dat"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer spin,L,N,i,j,imcs,maxL,nmcs 
      real*8 ecum,e2cum,xmcum,xm2cum,T,E,M,w,ratio 
      parameter (maxL=32)
      dimension spin(maxL,maxL),w(-4:4)
      ecum = 0.
      e2cum = 0.
      xmcum = 0.
      xm2cum = 0.
      open(unit=8,file='ising2.dat',status='unknown')
      call start(N,L,T,nmcs,spin,E,M,w)
      if(L.gt.maxL)call error(L,maxL)
      print*,' magnetizzazione iniziale netta per spin = ',(1.0*M)/N
      print*,' energia iniziale per spin = ',E/N
      write(8,*)
c     &' imcs, E/(N),M/(N),ecum/(imcs*N),xmcum/(imcs*N)'
      do 2001 j = 1,L
 2001    write(*,'(32i3.0)')(spin(i,j),i=1,L)
      do 100 imcs = 1,nmcs
         call Metrop(N,L,spin,E,M,w,ratio)
         call data(E,M,ecum,e2cum,xmcum,xm2cum,ratio)
c      write(8,*)
c     &'# T,  imcs,  E/N,  M/N,  ecum/(imcs*N),   xmcum/(imcs*N)'
c      write(8,*)
c     &'# 'T,imcs, E/(N),float(M)/(N),ecum/(imcs*N),xmcum/(imcs*N)
 100     continue
         call output(N,nmcs,ecum,e2cum,xmcum,xm2cum,ratio)
      do 201 j = 1,L
 201     write(*,'(32i3.0)')(spin(i,j),i=1,L)
      do 111 j = 1,L
      do 112 i = 1,L
         write(8,*)spin(i,j)
 112  continue
         write(8,*)
 111  continue
         close(8)
         stop
         end

      subroutine start(N,L,T,nmcs,spin,E,M,w)
      implicit none 
      integer N,L,nmcs,spin,idum,iconfig,iup,iright,i,j
      real*8 T,E,M,w,ran2,e4,e8,sum
      dimension spin(32,32),w(-4:4)
      write(6,*)' dimensione lineare del reticolo (un intero, L) = '
      read(5,*) L
      N = L*L ! numero di spin
      write(6,*)' temperatura (T) = '
      read(5,*) T
      write(6,*)' numero di steps Monte Carlo per spin (nmcs) = '
      read(5,*) nmcs
      write(6,*)' n. intero < 0 per generare numero random = '
      read(5,*) idum
c genera la configurazione iniziale
      print*,' inizio: spin down,a caso,up,meta`,antif (-1,0,1,2,3)? >'
      read(*,*)iconfig
      M = 0.
      E = 0.
      if ( iconfig .eq. 0) then
      do 200 j = 1,L
         do 100 i = 1,L
            if( ran2(idum) .lt. 0.5 ) then
               spin(i,j) = 1
               else
               spin(i,j) = -1
               end if
               M = M + spin(i,j)  ! magnetizzazione iniziale netta
 100           continue
 200           continue
      else if (iconfig .eq. -1 .or. iconfig .eq. 1) then   !CONFIG. ORDINATA
      do 201 j = 1,L
         do 101 i = 1,L
            spin(i,j) = iconfig
 101        continue
 201        continue
               M = N * iconfig
      else if (iconfig .eq. 2) then  !ALTRA CONFIG. ORDINATA
      do 211 j = 1,L
         do 111 i = 1,L/2
            spin(i,j) = 1
 111        continue
         do 311 i = L/2+1,L
            spin(i,j) = -1
 311        continue
 211        continue
               M = 0
      else if (iconfig .eq. 3) then  !ALTRA CONFIG. ORDINATA: ANTIFERROM.
      do 212 j = 1,L
         do 112 i = 1,L
            spin(i,j) = 1
            if( mod(i,2)+mod(j,2).eq.1 ) spin(i,j) = -1 !Devono essere due dispari
 112     continue
 212     continue
               M = 0
       else 
                  print*,' caso non previsto'
                  stop
       end if
c calcola l'energia iniziale sommando su tutte le coppie di NN
c (dato uno spin, solo con il NN sopra e il NN destra, 
c e non sotto e a sin., cosi' contiamo le interazioni una volta sola)
      do 300 j = 1,L
         if (j .eq. L) then
            iup = 1
            else
               iup = j + 1
               endif
               do 400 i = 1,L
                  if (i .eq. L) then
                     iright = 1
                     else
                        iright = i + 1
                        endif
                        sum = spin(i,iup) + spin(iright,j)
                        E = E - spin(i,j)*sum ! ENERGIA TOTALE iniziale
 400                    continue
 300                    continue
c calcola le probabilita' di transizione secondo la distrib. di Boltzmann
c (exp(-deltaE/KT)
c Avendo posto la forza delle interazioni (il param. J)=1, 
c le possibili variazioni di energia per uno spin flip sono -8,-4,0,+4,+8: 
c Argomento dell'exp  = cambiamento in energia
      e4 = exp(-4/T)
      e8 = e4 * e4
c argomento di w = isum = somma degli spin 4 NN (puo' essere -4,-2,0,+2,+4)
c Fissiamo il caso in cui lo spin centrale sia +1: allora
c se i 4 NN sono tutti +1, il suo flip corrisponde a deltaE = +8, e
c quindi associamo la prob. exp(-8/T):
      w(4)  = e8    
c e via dicendo per gli altri casi possibili 
c (se isum=0 inutile sprecar tempo, w(0)=1 q.que sia T)
      w(-4) = e8 !1/e8 ! cervellotico, ma tanto poi la usiamo solo per 
c             variaz. di E POSITIVE (le negative son comunque accettate)
      w(2)  = e4
      w(-2) = e4 !1/e4
      return
      end



      subroutine Metrop(N,L,spin,E,M,w,ratio)
      implicit none
      integer i,j,N,L,spin,idum,isum,ispin
      real*8 E,M,w,ratio,ran2
       dimension spin(32,32),w(-4:4)

      do 100 ispin = 1,N
c prende uno spin a caso (i,j sono le coordinate) 
         i = int (L*ran2(idum) + 1)
         j = int (L*ran2(idum) + 1)
c determina i valori degli spin NN usando le PBC
         call periodic(i,j,L,spin,isum)
         if (spin(i,j)*isum .le. 0) then  ! se la somma dei 4 spin NN a
                                          ! (i,j) ha segno opposto a
                                          ! spin(i,j) allora il deltaE
                                          ! per spin flip di (i,j) e' < 0
                                          ! e la mossa e' accettata
            call accept(i,j,M,E,isum,spin,ratio)
            else if ( ran2(idum) .lt. 
     &                w(-spin(i,j)*isum) ) then
               call accept(i,j,M,E,isum,spin,ratio) ! altrimenti, viene
                                                    ! accettata o no
                                                    ! secondo la
                                                    ! probabilita' di
                                                    ! transizione exp(-dE/T)
c cosi' l'argomento di w e' sempre <0 , =-2 o =-4
               end if
 100           continue
      return
      end
               

      subroutine periodic(i,j,L,spin,isum)
      implicit none
      integer i,j,spin,isum,L,left,iright,idown,iup
c determina la somma dei valori sui NN usando le PBC
      dimension spin(32,32)

      if (i .eq. 1) then
         left = spin(L,j)
         else
            left = spin(i-1,j)
            endif
      if (i .eq. L) then
         iright = spin(1,j)
         else
            iright = spin(i+1,j)
            endif
      if (j .eq. 1) then
         idown = spin(i,L)
         else
            idown = spin(i,j-1)
            endif
      if (j .eq. L) then
         iup = spin(i,1)
         else
            iup = spin(i,j+1)
            endif
      isum = left + iright + iup + idown
      return
      end
      


      subroutine accept(i,j,M,E,isum,spin,ratio)
      implicit none
      integer i,j,isum,spin
      real*8 M,E,ratio
c entrano M,E della config. prec. ed escono quelli della nuova config.
      dimension spin(32,32)
      spin(i,j) = -spin(i,j)
      ratio = ratio + 1        ! poi verra' rinormalizzata per
                               ! det. l'acc. ratio
      M = M + 2*spin(i,j)      ! qui ci va il 2 perche' e' la variaz.
      E = E - 2*spin(i,j)*isum ! idem qui
      return
      end


      subroutine data(E,M,ecum,e2cum,xmcum,xm2cum,ratio)
      implicit none
      real*8 E,M,ecum,e2cum,xmcum,xm2cum,ratio
c accumula dati dopo ogni step Monte Carlo per spin
      ecum = ecum + E
      e2cum = e2cum + E*E
      xmcum = xmcum + M
      xm2cum = xm2cum + M*M
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



      subroutine output(N,nmcs,ecum,e2cum,xmcum,xm2cum,ratio)
      implicit none
      real*8 ecum,e2cum,xmcum,xm2cum,ratio,znorm,eave,e2ave,save,s2ave
      integer*4 N,nmcs
      znorm = 1./(1.*nmcs*N)
c medie per spin
      ratio = ratio*znorm
      eave  = ecum*znorm
      e2ave = e2cum*znorm
      save  = xmcum*znorm
      s2ave = xm2cum*znorm
      write(6,*)' acceptance ratio = ',ratio
      write(6,*)' en. media per spin = ',eave
      write(6,*)' en. quadratica media per spin = ',e2ave
      write(6,*)' magnetizzazione media = ',save
      write(6,*)' magnetizzazione quadratica media = ',s2ave
      return
      end


      subroutine error(n,m)
      implicit none
      integer*4 n,m
      print*,' lato troppo grande (troppi spins)! (',n,'.ge.',m,')'
      stop
      end







