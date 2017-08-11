program hf

implicit none

double precision,dimension(6,6)::S,invmat,Fprime,P,F,G,Gprime,Fprimed,Cprime,Cmat
double precision,dimension(6)::eigval,eigvalFprime
double precision::E,oldE
double precision,dimension(6,6)::H
double precision,dimension(6,6,6,6)::G2
integer::i,j,mu,nu,lamda,sigma,m
integer::n=6,a
integer::nElec=4
open(unit=4, file='MULBERSHG.DAT',FORM='FORMATTED')
rewind 4
read(4,*) S, H, G2
close(4)
P=0.0
oldE=0.0
m=0

call MatDiag(S,eigval,n)
call inversematrix(S,n,invmat,eigval)
write(*,*)'invertmatrix', invmat

do
  m=m+1
  write(*,*)'SCF=',m

  call Gmat(n,P,G,G2)
  write(*,*) 'Gmatrix',F
  do i=1,6
         do j=1,6
                 F(i,j)=H(i,j)+G(i,j)
         end do
  end do
  write(*,*) 'Fockmatrix',F


  !calculation of Fprime
  Fprime=matmul(transpose(invmat),matmul(F,invmat))
  write(*,*) 'Fprime is=',Fprime
  call diaFprime(Fprime,eigvalFprime,n)
  write(*,*) 'Cprime',Fprime

  Cmat=matmul(invmat,Fprime)
  write(*,*)'Cmat=',Cmat
  !calculation of Pmatrix
  call Pmat(n,P,Cmat)     
  write(*,*)'new Pmat',P

  !calculation for energy
  E=0.0
  do mu=1,6
          do nu=1,6
                 E=E+0.5*(P(mu,nu)*((F(mu,nu))+H(mu,nu)))
          end do
  end do
  write(*,*) 'Energy=',E

  if (abs(oldE-E) .LE. 10.D-7) exit
  oldE=E
  
end do
contains
      subroutine overlapmat(S,eigval,n)
      implicit none
   
      integer n,l,info,lwork,lwmax

      double precision, intent(inout):: S(n,n)
      double precision, intent(out)::eigval(n)
      double precision:: work(n*(3+n/2))
      l=n*(3+n/2)
      lwmax=1000
      lwork=-1
      call dsyev('Vectors','Upper',n,S,n,eigval,work,lwork,l,info)
      lwork=min(LWMAX,int(work(1)))
      call dsyev('Vectors','Upper',n,S,n,eigval,work,lwork,l,info)
      end subroutine overlapmat
  
      subroutine inversematrix(S,n,invmat,eigval)
      implicit none
      integer::i
      integer::n
      double precision, intent(inout):: S(n,n)
      double precision, intent(in)::eigval(n)
      double precision,dimension(6,6),intent(inout)::invmat
     
      invmat=0.d0
      do i=1,6
          invmat(i,i)=(1.d0/dsqrt(eigval(i)))
      end do
      invmat=matmul(S,matmul(invmat,transpose(S)))

      end subroutine inversematrix
      
      subroutine Pmat(n,P,Cmat)
           integer::mu,nu,sigma,lamda,n
           double precision,intent(out)::P(n,n)
           double precision,intent(in)::Cmat(n,n)
           P=0.0
           do mu=1,6
                  do nu=1,6
                     do a=1,2
                          P(mu,nu)=P(mu,nu)+2*(Cmat(mu,a)*Cmat(nu,a))
                     end do
                  end do
           end do
      end subroutine Pmat

      subroutine Gmat(n,P,G,G2)
           integer::mu,nu,sigma,lamda,n
           double precision :: Gprime
           double precision,intent(in)::P(n,n)
           double precision,intent(out)::G(n,n)
           double precision,intent(in)::G2(n,n,n,n) 
           G=0.0
           do mu=1,6
             do nu=1,6
               do lamda=1,6 
                do sigma=1,6
                  G(mu,nu)=G(mu,nu)+(P(lamda,sigma)*(G2(mu,nu,sigma,lamda)-0.5*G2(mu,lamda,sigma,nu)))
                end do
               end do
             end do
           end do
      end subroutine Gmat
 
      subroutine diaFprime(Fprime,eigvalFprime,n)

      implicit none
  
      integer n,l,info,lwork,lwmax
 
      double precision, intent(out)::Fprime(n,n)
      double precision, intent(out)::eigvalFprime(n)
      double precision::work(n*(3+n/2))
      l=n*(3+n/2)
      lwmax=1000
      lwork=-1
      call dsyev('Vectors','Upper',n,Fprime,n,eigvalFprime,work,lwork,l,info)
      lwork=min(LWMAX,int(work(1)))
      call dsyev('Vectors','Upper',n,Fprime,n,eigvalFprime,work,lwork,l,info)
 
      end subroutine diaFprime
    subroutine MatDiag(Matrix,EigenValues,n)
        implicit none
	    double precision  Matrix(n,n),EigenValues(n),work(n*(3+n/2))
        integer::l,INF,n
	    l=n*(3+n/2)
	    call dsyev('V','U',n,Matrix,n,EigenValues,work,l,inf)
    end subroutine MatDiag
end program hf
