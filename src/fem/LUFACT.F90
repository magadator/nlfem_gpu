   subroutine  LUFACT(a,na,nn,iperm,v,ising) 
!
!  TYPE      AND     DIMENSION 
!
     integer  na, nn, iperm(*), ising 
     double precision a(na,*), v(*), mactol, maceps, x, y 
     integer  i, j, k, l, n 
     data maceps /2.5d-17/
!
!           LOGIC 
!         Initiallzation 
!
!     write(*,*) "1"
     mactol =    8*maceps 
     ising  =    0 
     n  =        iabs(nn) 
     if (n .le. 0) return 
!
!    Calculate the Euclidean norm of each row of A 
!
!     write(*,*) "1a"
     do i =  1,n 
        y  =  0.0 
        do j =  1,n 
           y = y  +  a(i,j)**2 
        end do
        v(i) = 0.0 
        if (y  .ne. 0.d0) v(i) = sqrt(1.0/y) 
     end do
! 
!    Main loop 
! 
!     write(*,*) "2"
     do 4000 k =  1,n 
        if (v(k) .eq. 0.0) then 
           iperm(k) =  k 
           go to 4000
        end if
        l = k 
        x = 0.0 
        do i = k,n 
           y= 0.0 
           do j = 1,k-1
              y  =        y  +  a(i,j) * a(j,k) 
           end do
           a(i,k) = a(i,k) - y 
           y  =     abs(v(i)*a(i,k)) 
           if (y .gt. x) then 
              x = y 
              l = i 
           end if
        end do
! 
!       Interchange rows 
! 
!        write(*,*) "2"
        if (l  .ne. k) then 
           do j =  1,n 
              y      =   a(k,j) 
              a(k,j) =   a(l,j) 
              a(l,j) =   y 
           end do
        end if
        v(l) = v(k) 
        iperm(k) =  l
!
!       Test for singularity 
!
        if  (x  .le. mactol)       then 
           ising =  n - (k-1) 
           return 
        end if
        x =   1.0/a(k,k) 
        a(k,k) = x 
!
!       Calculate elements of strict upper trlangle 
!
        do j = k+1,n 
           y = 0.0 
           do i = 1,k-1 
              y =      y  +  a(k,i) * a(i,j) 
           end do
           a(k,j) = (a(k,j)-y) * x 
        end do
 4000    continue
!      
     return 
   end subroutine LUFACT
