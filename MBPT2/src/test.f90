 program diatest
!---------------!
 implicit none

 integer :: i,n
 real(8), allocatable :: m0(:,:),m1(:,:),m2(:,:),eig(:)

 open(10,file='mat.dat',status='old')
 read(10,*)n
 allocate (m0(n,n))
 allocate (m1(n,n))
 allocate (m2(n,n))
 allocate (eig(n))
 do i=1,n
    read(10,*)m0(:,i)
 enddo
 close(10)

 m1(:,:)=m0(:,:)

 call diasym(m1,eig,n)
   
 open(10,file='dia.dat',status='replace')
 write(10,*)'Eigenvalues:'
 do i=1,n
    write(10,10)i,eig(i)
    10 format(I3,'   ',f14.8)
 enddo
 write(10,*)
 write(10,*)'Eigenvectors:'
 do i=1,n
    write(10,20)i,m1(:,i)
    20 format(i3,'   ',10f14.8)
 enddo
 write(10,*)
 
 m2=matmul(transpose(m1),m0)
 m0=matmul(m2,m1)

 write(10,*)'Transformed matrix (check):'
 do i=1,n
    write(10,30)m0(:,i)
    30 format(10f14.8)
 enddo
 write(10,*)

 close(10)

 deallocate(m0); deallocate(m1); deallocate(m2);  deallocate(eig)

 end program diatest
!-------------------!

!---------------------------------------------------------!
!Calls the LAPACK diagonalization subroutine DSYEV        !
!input:  a(n,n) = real symmetric matrix to be diagonalized!
!            n  = size of a                               !
!output: a(n,n) = orthonormal eigenvectors of a           !
!        eig(n) = eigenvalues of a in ascending order     !
!---------------------------------------------------------!
 subroutine diasym(a,eig,n)
 implicit none

 integer n,l,inf
 real*8  a(n,n),eig(n),work(n*(3+n/2))

 l=n*(3+n/2)
 call dsyev('V','U',n,a,n,eig,work,l,inf)

 end subroutine diasym
!---------------------!
