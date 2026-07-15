    subroutine SM3MB(x,y,dm,alpha,f,ls,sm,m) 
! 
!     TYPE  &  DIMENSION 
! 
      implicit none
      !
      integer            m,ls(9) 
      real*8   x(3), y(3),dm(3,3), alpha, f,  p(9,3), smOut(m*m), sm(m,m) 
      real*8   area2,c 
      real*8   d11, d12, d13, d22, d23, d33 
      real*8   x21, x32, x13, y21, y32, y13 
      real*8   x12, x23, x31, y12, y23, y31 
      real*8   s1, s2, s3 
      integer            i,j,k,l,n 
      character*120      status       
!
!                       LOGIC 
!
!      write(*,*) "[I] SM3MB"
      status = ' '
      x21 = x(2) -  x(1)
      x12 =-x21
      x32 = x(3) -  x(2)
      x23 =-x32
      x13 = x(1) -  x(3)
      x31 =-x13
      y21 = y(2) -  y(1)
      y12 =-y21
      y32 = y(3) -  y(2)
      y23 =-y32
      y13 = y(1) -  y(3)
      y31 =-y13

      area2 =     y21*x13 - x21*y13 
      if (area2 .le. 0.0) then 
        status =  'NEGA_AREA'
        if (area2 .eq. 0.0)  status =  'ZERO_AREA' 
        return 
      end if 
      p(1,1) = y23 
      p(2,1) = 0.0 
      p(3,1) = y31 
      p(4,1) = 0.0 
      p(5,1) = y12 
      p(6,1) = 0.0 
      p(1,2) = 0.0 
      p(2,2) = x32 
      p(3,2) = 0.0 
      p(4,2) = x13 
      p(5,2) = 0.0 
      p(6,2) = x21 
      p(1,3) = x32 
      p(2,3) = y23 
      p(3,3) = x13 
      p(4,3) = y31
      p(5,3) = x21 
      p(6,3) = y12 
      n=6 

      if (alpha.ne.0.0) then 
        p(7,1) =  y23*(y13-y21)*alpha/6. 
        p(7,2) =  x32*(x31-x12)*alpha/6.
        p(7,3) =  (x31*y13-x12*y21)*alpha/3. 
        p(8,1) =  y31*(y21-y32)*alpha/6. 
        p(8,2) =  x13*(x12-x23)*alpha/6. 
        p(8,3) =  (x12*y21 -x23*y32)*alpha/3. 
        p(9,1) =  y12*(y32-y13)*alpha/6. 
        p(9,2) =  x21*(x23-x31)*alpha/6. 
        p(9,3) =  (x23*y32-x31*y13)*alpha/3.
        n  = 9 
      end if 

      c   = 0.5*f/area2 
      d11 = c  * dm(1,1) 
      d22 = c  * dm(2,2) 
      d33 = c  * dm(3,3) 
      d12 = c  * dm(1,2) 
      d13 = c  * dm(1,3) 
      d23 = c  * dm(2,3) 

      do j  =  1,n 

         l  =  ls(j) 
         s1 =  d11*p(j,1) +  d12*p(j,2) +  d13*p(j,3) 
         s2 =  d12*p(j,1) +  d22*p(j,2) +  d23*p(j,3) 
         s3 =  d13*p(j,1) +  d23*p(j,2) +  d33*p(j,3) 
         do i =  1,j 

            k = ls(i)
            sm(k,l) =  sm(k,l) +  (s1*p(i,1) + s2*p(i,2) + s3*p(i,3)) 
            sm(l,k) =  sm(k,l) 
         end do
      end do
      
      return 
   end subroutine SM3MB
