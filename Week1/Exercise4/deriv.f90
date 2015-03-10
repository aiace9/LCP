program deriv
!c
!c derivata numerica sin, ds, simmetrica IN SINGOLA PRECISIONE
!c
real :: h(8) 
real :: x, exact
integer :: i, N=8 
data h/0.1, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00001, 0.000001 /
!c
print*, "    h,    f'_ds,      error,      f'_sin,     error,    f'_simm,    error "
x = 1.0                   ! inizializzo le variabili
exact = cos(x)
do i=1,N
   deriv_ds = (sin(x+h(i))-sin(x))/h(i)
   deriv_sin = (sin(x)-sin(x-h(i)))/h(i)
   deriv_simm = (sin(x+h(i))-sin(x-h(i)))/(2*h(i))
   print*, '--------'
   write(*,'(7(f9.6,2x))'),  h(i), deriv_ds, deriv_ds - exact, deriv_sin, deriv_sin - exact, deriv_simm, deriv_simm - exact
   print*, sin(x),sin(x-h(i)),(sin(x)-sin(x-h(i))),h(i)
end do                 

stop
end program deriv
      
