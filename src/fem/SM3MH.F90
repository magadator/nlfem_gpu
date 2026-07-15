   subroutine SM3MH(x,y,dm,f,ls,sm,m,status) 
!
!              TYPE  k  DIMENSION 
!
     implicit none
     !
     character*(*) status                                 
     integer       ls(9), m 
     double precision x(3), y(3), dm(3,3), f, sm(m,m) 
     double precision xc(3), yc(3), xm(3), ym(3) 
     double precision bh(3,3), gt(9,9), hh(3,9), sqh(3,3), t(9) 
     double precision qx(3,3), qy(3,3)  
     double precision area, area2 
     double precision a1j, a2j,a3j,b1j,b2j,b3j 
     double precision c, cj, sj, dl, dx, dy 
     double precision d11, d12, d13, d22, d23, d33, jxx, jxy, jyy 
     double precision s1, s2, s3, s4, s5, s6 
     double precision x0, y0, xi, yi 
     integer  i, j, k, l, ising, iperm(9)
!
!    LOGIC 
!
!     write(*,*) "[I] SM3MH"
     !
     status =        '  ' 
     area2  = (y(2)-y(1))*(x(1)-x(3)) -  (x(2)-x(1))*(y(1)-y(3)) 
     if  (area2.le.0.0) then 
        status =  'NEGA_AREA' 
        if  (area2.eq.0.0)  status = 'ZERO_AREA' 
        return 
     end if
!     write(*,*) "A"
     x0  = (x(1)+x(2)+x(3))/3.0 
     y0  = (y(1)+y(2)+y(3))/3.0 
     area = 0.5*area2 
     c  =   1./sqrt(area) 
     xc(1) =   c  *  (x(1)-x0) 
     xc(2) =   c  *  (x(2)-x0) 
     xc(3) =   c  *  (x(3)-x0) 
     yc(1) =   c  *  (y(1)-y0) 
     yc(2) =   c  *  (y(2)-y0) 
     yc(3) =   c  *  (y(3)-y0) 
     xm(1) =   0.5 *  (xc(2)+xc(3)) 
     xm(2) =   0.5 *  (xc(3)+xc(1)) 
     xm(3) =   0.5 *  (xc(1)+xc(2)) 
     ym(1) =   0.5 *  (yc(2)+yc(3)) 
     ym(2) =   0.5 *  (yc(3)+yc(1)) 
     ym(3) =   0.5 *  (yc(1)+yc(2)) 
!
!     Form G'  (G transposed) in GT and initialize HH 
!
!     write(*,*) "B"
     do i = 1,9 
        do j = 1,6 
           gt(j,i) =  0 
        end do
        hh(1,i) =  0. 
        hh(2,i) =  0. 
        hh(3,i) =  0. 
     end do
!
!     write(*,*) "C"
     d11 =   f  * dm(1,1) 
     d22 =   f  * dm(2,2) 
     d33 =   f  * dm(3,3) 
     d12 =   f  * dm(1,2) 
     d13 =   f  * dm(1,3) 
     d23 =   f  * dm(2,3) 
     jxx =   -2.*(xc(1)*xc(2)+xc(2)*xc(3)+xc(3)*xc(1))/3. 
     jxy =       (xc(1)*yc(1)+xc(2)*yc(2)+xc(3)*yc(3))/3. 
     jyy =   -2.*(yc(1)*yc(2)+yc(2)*yc(3)+yc(3)*yc(1))/3. 
!     write(*,*) "D"
     do j =  1,3 
        dx =  xm(j) - xc(j) 
        dy =  ym(j) - yc(j) 
        dl =  sqrt(dx**2 + dy**2) 
        cj =  dx/dl 
        sj =  dy/dl 
        a1j =  -0.5*sj*cj**2 
        a2j =   0.5*cj**3
        b2j =  -0.5*sj**3
        b3j =   0.5*sj**2*cj
        a3j =  -(b2j +  a1j +  a1j)
        b1j =  -(b3j +  b3j +  a2j)
        gt(1,2*j-1) =     1. 
        gt(2,2*j )  =     1. 
        gt(3,2*j-1) =  -yc(j) 
        gt(3,2*j  ) =   xc(j) 
        gt(3,  j+6) =   c 
        gt(4,2*j-1) =   xc(j) 
        gt(6,2*j-1) =   yc(j) 
        gt(5,2*j  ) =   yc(j) 
        gt(6,2*j  ) =   xc(j) 
        hh(j,j+6  ) =   1 
        qx(j,1)     =   a1j 
        qx(j,2)     =   b2j 
        qx(j,3)     =  -2.*b3j 
        qy(j,1)     =   a2j 
        qy(j,2)     =   b3j 
        qy(j,3)     =  -2*a1j 
        s1 =       d11*qx(j,1)+d12*qx(j,2) + d13*qx(j,3) 
        s2 =       d12*qx(j,1)+d22*qx(j,2) + d23*qx(j,3) 
        s3 =       d13*qx(j,1)+d23*qx(j,2) + d33*qx(j,3) 
        s4 =       d11*qy(j,1)+d12*qy(j,2) + d13*qy(j,3) 
        s5 =       d12*qy(j,1)+d22*qy(j,2) + d23*qy(j,3)
        s6 =       d13*qy(j,1)+d23*qy(j,2) + d33*qy(j,3) 
        do i = 1,3 
           xi = xc(i) 
           yi = yc(i) 
           gt(j+6,2*i-1) = a1j*xi*xi + 2 *a2j*xi*yi + a3j*yi*yi 
           gt(j+6,2*i  ) = b1j*xi*xi + 2 *b2j*xi*yi + b3j*yi*yi 
           gt(j+6,i+6  ) =-c*(cj*xi+sj*yi) 
        end do
        do i =  1,j
           sqh(i,j) = jxx *  (qx(i,1)*s1+qx(i,2)*s2+qx(i,3)*s3)&
                   +  jxy *  (qx(i,1)*s4+qx(i,2)*s5+qx(i,3)*s6 &
                   +          qy(i,1)*s1+qy(i,2)*s2+qy(i,3)*s3)&
                   +  jyy *  (qy(i,1)*s4+qy(i,2)*s5+qy(i,3)*s6) 
        end do
     end do
!     write(*,*) "E"
! 
!               Factor G'  and backsolve to obtain H 
!                                                  h 
     call  LUFACT (gt, 9, 9,  iperm, t,  ising) 
     if  (ising .ne.  0) then 
        status =  'SINGULAR_G' 
        return 
     end if
!     write(*,*) "F"
     call  LUSOLV (gt, 9, 9,  iperm, hh, 3, -3) 
!
!    Form physical stiffness and add to incoming SM !
!     write(*,*) "G"
     do j = 1,9 
        l  = ls(j) 
        s1 = sqh(1,1)*hh(1,j) + sqh(1,2)*hh(2,j) + sqh(1,3)*hh(3,j) 
        s2 = sqh(1,2)*hh(1,j) + sqh(2,2)*hh(2,j) + sqh(2,3)*hh(3,j) 
        s3 = sqh(1,3)*hh(1,j) + sqh(2,3)*hh(2,j) + sqh(3,3)*hh(3,j) 
        do i = 1,j 
           k  = ls(i) 
           sm(k,l) =  sm(k,l) +  (s1*hh(1,i) + s2*hh(2,i) + s3*hh(3,i))
           sm(l,k) =  sm(k,l)
        end do
     end do
     return 
   end subroutine SM3MH

