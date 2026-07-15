   subroutine LUSOLV(a,  na, nn, iperm, b, nb, m) 
!
!  TYPE    AND    DIMENSION 
!
     integer        na, nn, iperm(*), nb, m 
     double precision a(na,*), b(nb,*), t, sum 
     logical     transp 
     integer     i,  j,  k, n, nrhs 
!
!    LOGIC 
!    Initlallzation 
!
     nrhs =    iabs(m) 
     transp =  .false. 
     if (m  .lt.  0) transp =  .true. 
     n  =      iabs(nn) 
!
!    Interchange elements of rhs matrix B 
!
     do i = 1,n 
        k  = iperm(i) 
        if  (k.ne.i) then 
           if ( .not. transp) then 
              do j = 1,nrhs 
                 t  =      b(i,j) 
                 b(i,j) =  b(k,j) 
                 b(k,j) =  t 
              end do
           else 
              do j =  1,nrhs 
                 t =  b(j ,i) 
                 b(j,i) =  b(j,k) 
                 b(j,k) =  t 
              end do
           end if
        end if
     end do
! 
!    Cycle over right hand sides 
! 
     do j  =  1,nrhs 
        if (.not. transp) then 
!
!       Solve for rhs  column storage 
!
           b(1,j) = b(1,j) * a(1,1) 
           do i = 1,n-1
              sum =    0.0 
              do k = 1,i
                 sum =    sum  -  a(i+1,k) * b(k,j) 
              end do
              b(i+1,j) = (b(i+1,j) + sum) * a(i+1,i+1) 
           end do
           do i = n-1,1,-1 
              sum = 0.0
              do k  =  1,n-i 
                 sum = sum - a(i,i+k) * b(i+k,j) 
              end do
              b(i, j)  = b(i,j) + sum 
           end do
        else 
! 
!             Solve for rhs' row storage 
! 
           b(j,1) =    b(j,1) * a(1,1) 
           do i = 1,n-1 
              sum =     0.0 
              do k  =  1,i 
                 sum = sum - a(i+1,k) * b(j,k) 
              end do
              b(j,i+1) =  (b(j,i+1) + sum) * a(i+1,i+1) 
           end do
           do i = n-1,1,-1 
              sum = 0.0 
              do k = 1,n-i 
                 sum = sum - a(i,i+k) * b(j,i+k) 
              end do
              b(j,i) =   b(j,i) + sum 
           end do
        end if
     end do
     return 
   end subroutine LUSOLV
