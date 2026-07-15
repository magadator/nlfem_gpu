subroutine SM4M2(x, y, dmIn, alpha, beta, ls, smOut, m)
!                                                                                                                                                                                              
!    TYPE  &  DIMENSION                                                                                                                                                                        
!                                                                                                                                                                                              
  implicit none
  !                                                                                                                                                                                            
     character*120 status 
     integer     iat, m, ls(12) 
     real*8 dmIn(9),smOut(m*m)
     real*8 x(3), y(3), dm(3,3), alpha, beta, sm(m,m) 
     real*8 xt(3), yt(3), dmt(3,3), f 
     integer    i, j,  n, t 
!                                                                                                                                                                                              
!                                 LOGIC                                                                                                                                                        
!               
     sm=0.d0
      do j=1,3
         do i=1,3
            dm(i,j)=dmIn(i+3*(j-1));
         end do
      end do
      status =  ' '
      do i = 1,3
         do j = 1,3
            dmt(i,j) =  dm(i,j)
         end do
      end do
!                                                                                                                                                                                              
      f =    1.0
!                                                                                                                                                                                              
      call  SM3MB (x, y, dmt, alpha, f,  ls,  sm,  m,  status)
      if  (beta .ne. 0.0)  then
         call  SM3MH (x, y, dmt, f*beta, ls, sm, m, status)
      end if
      !if  (status(1:1) .ne.  ' ') return



!  write(*,*) "SM4M"
!  do j=1,m
!     do i=1,m
!        sm(i,j)=j+(i-1)*m;
!     end do
!  end do

  !
!  do n=1,9
!     write(*,'(9E15.7)') sm(n,1:9)
!  end do

! need translate 2D array into 1D array
  do j=1,m
     do i=1,m
        smOut(i+m*(j-1))=sm(j,i);
     end do
  end do
!  write(*,*) smOut
  return
end subroutine SM4M2











