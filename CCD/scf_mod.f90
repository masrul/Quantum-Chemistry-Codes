!****************************************************************************************
! By: Zachary Windom
! Date: 7/18/2017
! Description: Program reads a file containing core Hamiltonian matrix, overlap (S) matrix,
!              and two electron integrals. Hartree Fock SCF energy is then calculated
!              using 6 A.O. basis function for a single Be atom. Page 146 of
!              text was used extensively, as well as class notes from 6/13/2017.
!              Calculates MBPT2 correlation energy after transforming AO 2 e-
!              integrals to MO 2 e- integrals.
!***************************************************************************************
!syntax for compiling/linking: 
!gfortran -o fcode scf.f90      -L/usr/local/lib    -llapack -lblas
!syntax found at:
!https://www.math.utah.edu/software/lapack.html
program SCF
IMPLICIT NONE
INTEGER:: alpha,beta,i,j,l,m,k,n,o,p,ierror,ok,lwork,counter,istat
INTEGER, PARAMETER :: LWMAX=17 !syntax found pg. 857
DOUBLE PRECISION :: x, y, z, var,term, energy, energyTerm,orbitalEnergy(6)
DOUBLE PRECISION :: test(6,6),S(6,6), H(6,6), G2(6,6,6,6), diaS(6), work(LWMAX)
DOUBLE PRECISION :: dummy2(LWMAX),gsum,summed,total,quarter_trans(6,6,6,6)
DOUBLE PRECISION :: initP(6,6), initG(6,6), initF(6,6)
DOUBLE PRECISION :: A(6,6), B(6,6), C(6,6), primeF(6,6), diaPrimeF(6),so_energy(12)
DOUBLE PRECISION :: primeC(6,6), Gmatrix(6,6), oldE, newE, diff,energy_den,num
DOUBLE PRECISION :: oldP(6,6), newP(6,6),reg(6,6), rootmeansq,trans(6,6,6,6),trans_so(12,12,12,12)
LOGICAL :: success
integer::q,r,ss,aa,bb
integer::mu,nu,lamda,sigma
INTEGER, PARAMETER :: out_unit=20
OPEN(UNIT=2, FILE='BERSHG.DAT', STATUS='OLD', ACTION='READ', IOSTAT=ierror)
OPEN(UNIT=out_unit, FILE="MOcoeff.dat", STATUS='REPLACE', ACTION='WRITE')
REWIND 2
READ(2,*) S, H, G2
CLOSE(UNIT=2)
!******************************************************
! STEP 1 ::::: Diagonalize S and build S^-1/2 ; Use eigenvectors/values to build
!              X
!******************************************************
! ex. of 'DSYEV' usage found at: 
! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dsyev_ex.f.htm
lwork=-1
call DSYEV('V','U',6,S,6,diaS,work,lwork,ok)
!call PRINT_MATRIX('ORIGINAL S::::',1,6,diaS,1)
lwork=min(LWMAX,int(work(1)))
call DSYEV('V','U',6,S,6,diaS,work,lwork,ok)
!CALL PRINT_MATRIX('EIGENVECTORS:::',6,6,S,6)
!CALL PRINT_MATRIX('EIGENVALUES of S', 1,6, diaS,1)
!*****************************************************
! STEP 2 ::::: Obtain s^-1/2 from eigenvalues of S
DO i=1,6
    DO j=1,6
        if (i==j) THEN
            test(i,j)=1.0/sqrt(diaS(i))
            !write(*,*) test(i,j)
        else
            test(i,j)=0.0
        endif
    END DO
END DO
!CALL PRINT_MATRIX('s^-1/2',6,6,test,6)
!**********************************************
! STEP 3 ::::: Determines transformation matrix X = U*SU
!              from newly formed eigenvectors matrix, 'S',
!              and diagonal, inverted, sq root matrix composed of eigenvalues
!              of S, s^-1/2
reg=MATMUL(S,test)
reg=MATMUL(reg,transpose(S))
!CALL PRINT_MATRIX('X MATRIX:::',6,6,reg,6)
!**********************************************
! STEP 4 ::::: Calculate Fock Matrix from H
!              and G matrices
diff =1
counter=0
success=.FALSE.
DO
!    IF ((diff<=10E-06) .OR. (counter>20)) EXIT
    IF (success .eqv. .TRUE.) EXIT  
    IF (counter==0) THEN
        ! ***** First iteration, G matrix is 0 due to P
        !       being 0
        DO i=1,6
            DO j=1,6
                initF(i,j)=H(i,j)
                newP(i,j)=0
                Gmatrix(i,j)=0    
            END DO
        END DO
    ELSE
        DO i=1,6
            DO j=1,6
                gsum=0
                DO k=1,6
                    DO l=1,6
                        gsum=gsum+newP(k,l)*(G2(i,l,j,k)-(0.5)*G2(i,l,k,j))
                    END DO
                END DO
                Gmatrix(i,j)=gsum
                initF(i,j)=H(i,j)+Gmatrix(i,j)
            END DO
        END DO
    END IF
    counter=counter+1
!    CALL PRINT_MATRIX('G MATRIX::', 6,6, Gmatrix,6)
!    CALL PRINT_MATRIX('FOCK matrix',6,6,initF,6)                        
!**********************************************
! STEP 5 ::::: Calculate Fock PRIME matrix
!              by F'=(s^-1/2)* F s^-1/2
! source:::
! https://software.intel.com/en-us/node/684733
    alpha =1
    beta =1 
    DO i=1,6
        DO j=1,6
            C(i,j)=0
            B(i,j)=0
        END DO
    END DO
    C = MATMUL(initF,reg)
    primeF = MATMUL(transpose(reg),C)
!    CALL PRINT_MATRIX('PRIME FOCK MATRIX',6,6,primeF,6)
!**********************************************
! STEP 6 ::::: Diagonalize primeF to get primeC
!              IE eigenvectors of primeF == primeC
    lwork=-1
    call DSYEV('V','U',6,primeF,6,diaPrimeF,dummy2,lwork,ok)
    lwork=min(LWMAX,int(dummy2(1)))
    call DSYEV('V','U',6,primeF,6,diaPrimeF,dummy2,lwork,ok)
!    CALL PRINT_MATRIX('EIGENVALUES OFPRIME F', 1,6, diaPrimeF,1)
!    CALL PRINT_MATRIX('EIGENVectors OF PRIME F', 6,6, primeF,6)
    do i=1,6
        do j=1,6
            primeC(i,j)=primeF(i,j)
        end do
    end do
!**********************************************
! STEP 7 ::::: Find C matrix by C= s^-1/2 * primeC
    B = MATMUL(reg, primeC)
!    CALL PRINT_MATRIX('real C matrix', 6,6, B,6)
     
!**********************************************
! STEP 8 :::: Form new P matrix from C matrix
!             IE P_munu=2*sum(C_mui*C_nui)
!                   sum from i=1 to 2
! ****** NOTE *******
! Fortran(gfortran) does not like indices other than
! 'i','j','k',...etc .. whY?
    do i=1,6
        do j=1,6
            oldP(i,j)=newP(i,j)
            term = 0
            do k=1,2
                term=term+B(i,k)*B(j,k)
            end do
            term = 2*term
            newP(i,j)=term
        end do 
    end do    
!    CALL PRINT_MATRIX('P MATRIX',6,6,newP,6)
!**********************************************
! STEP 9 ::::: CHECK FOR CONVERGENCE OF  !!!!!!
    CALL ROOT_MEAN_SQUARE(oldP, newP,6,success)
!**********************************************
! STEP 10 ::::: Calculate energy
    if (counter==0) THEN
        oldE=0
    ELSE
        oldE=newE
    endif
    energyTerm=0
    do i=1,6
        do j=1,6
            energyTerm=energyTerm+ newP(j,i)*(H(i,j)+initF(i,j))
        end do
    end do
    energyTerm=(1.0/2)*energyTerm
    newE=energyTerm
!    write(*,*) 'Iteration number is:  ', counter
!    write(*,*) 'energy is: ', newE
    diff = abs(newE-oldE)        
END DO
write(*,*)
write(*,*) 'Final converged energy is:'
write(*,*) '**************************'
write(*,*) newE
write(*,*) 'Energy converged on the',counter,' iteration.'
write(*,*)

              
!trans=0.d0
!DO p=1,6
!  DO q=1,6
!    DO r=1,6
!      DO ss=1,6
!        do mu=1,6
!          do nu=1,6
!            do lamda=1,6
!              do sigma=1,6
!               trans(p,q,r,ss)=trans(p,q,r,ss)+B(p,mu)*B(q,nu)*G2(mu,nu,lamda,sigma)*B(r,lamda)*B(ss,sigma)
!              end do
!            end do
!          end do
!         end do
!       end do
!     end do
!   end do
!end do
!**********************************************
! STEP 11 :::::  Transform from AO to MO

trans=0.d0
DO i=1,6
  DO j=1,6
    DO aa=1,6
      DO bb=1,6
        do mu=1,6
          do nu=1,6
            do lamda=1,6
              do sigma=1,6
               trans(i,j,aa,bb)=trans(i,j,aa,bb)+B(mu,i)*B(nu,aa)*G2(mu,lamda,nu,sigma)*B(lamda,j)*B(sigma,bb)
              end do
            end do
          end do
         end do
       end do
     end do
   end do
end do




!
!!**********************************************
!! STEP 12 :::::  Find MP2 correlation energy 
!!                correction to HF energy
!
!total=0.0
!do i=1,2
!  do j=1,2
!    do aa=3,6
!      do bb=3,6
!        num=trans(i,j,aa,bb)*(2*trans(i,j,aa,bb)-trans(i,j,bb,aa))
!        energy_den=diaPrimeF(i)+diaPrimeF(j)-diaPrimeF(aa)-diaPrimeF(bb)
!        total=total+(num/energy_den)
!      end do
!    end do
!  end do
!end do
!write(*,*)'!!!! MP2 CORRELATION ENERGY IS: !!!!! ', total
!write(*,*)'MBPT(2) energy is: ', total+newE
!call Joshua_Spatial2Spin(trans,trans_so)
call Spatial2Spin(trans,trans_so)
!call another_Spatial2Spin(trans,trans_so)

do i=1,6
    alpha=2*i-1
    beta=alpha+1
    so_energy(alpha)=diaPrimeF(i)
    so_energy(beta)=diaPrimeF(i)
end do
!print*,so_energy

total=0.0
do i=1,4
    do j=1,4
        do aa=5,12
            do bb=5,12
                !num=trans_so(i,j,aa,bb)-trans_so(i,j,bb,aa)
                num=trans_so(i,j,aa,bb)
                num=num*num
                energy_den=so_energy(i)+so_energy(j)-so_energy(aa)-so_energy(bb)
                total=total+0.25*(num/energy_den)
            end do 
        end do 
    end do 
end do 
write(*,*)'!!!! MP2 CORRELATION ENERGY IS: !!!!! ', total
END program SCF













!!! DEF: Prints out matrix elements
SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
  CHARACTER*(*)    DESC
  INTEGER          M, N, LDA
  DOUBLE PRECISION A( LDA, * )
  INTEGER          I, J
  WRITE(*,*)
  WRITE(*,*) DESC
  DO I = 1, M
  WRITE(*,9998) ( A( I, J ), J = 1, N )
  END DO
  9998 FORMAT( 11(:,1X,2E14.4) )
  RETURN
END SUBROUTINE PRINT_MATRIX

!!! DEF: Calculates the RMS of each element of 
!!!      the old and new P (density) matrices
!!!      Returns a bool (T/F) depending on convergence
!!! Source for loop structure :::
!!! https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/297335
SUBROUTINE ROOT_MEAN_SQUARE(prevP, currentP, N,success)
   DOUBLE PRECISION, INTENT(in) :: prevP(N,N), currentP(N,N)
   DOUBLE PRECISION :: var,summed_values
   LOGICAL, INTENT(out):: success
   INTEGER, INTENT(in) :: N
   summed_values=0
   loop1:DO i=1,N
      loop2:DO j=1,N
         var = (prevP(N,N)-currentP(N,N))**2
         summed_values=summed_values+var
      END DO loop2
   END DO loop1
   var = sqrt(summed_values/N**2)
   if (var <= 10E-08) THEN
      success=.TRUE.
   else
      success=.FALSE.
   end if 
!   write(*,*)
!   write(*,*) '**** DENSITY MATRIX, P, HAS CONVERGED? ****'
!   write(*,*) success
   RETURN 
END SUBROUTINE ROOT_MEAN_SQUARE 

subroutine Spatial2Spin(G2,G2_SO)
    double precision::G2(6,6,6,6)
    double precision::G2_SO(12,12,12,12)

    integer::i,j,k,l
    integer::alpha1,beta1    
    integer::alpha2,beta2    
    integer::alpha3,beta3    
    integer::alpha4,beta4
    G2_SO=0.d0
    do i=1,6
        alpha1=2*i-1
        beta1=alpha1+1
        do j=1,6
            alpha2=2*j-1
            beta2=alpha2+1
            do k=1,6
                alpha3=2*k-1
                beta3=alpha3+1
                do l=1,6
                    alpha4=2*l-1
                    beta4=alpha4+1

                    G2_SO(alpha1,alpha2,alpha3,alpha4)=G2(i,j,k,l)-G2(i,j,l,k)
                    G2_SO(alpha1,beta2,alpha3,beta4)=G2(i,j,k,l)
                    G2_SO(alpha1,beta2,beta3,alpha4)=-G2(i,j,l,k)
                    G2_SO(beta1,alpha2,alpha3,beta4)=-G2(i,j,l,k)
                    G2_SO(beta1,alpha2,beta3,alpha4)=G2(i,j,k,l)
                    G2_SO(beta1,beta2,beta3,beta4)=G2(i,j,k,l)-G2(i,j,l,k)
                end do 
            end do 
        end do 
    end do 
end subroutine Spatial2Spin

subroutine another_Spatial2Spin(G2,G2_SO)
    implicit none 
    double precision::G2(6,6,6,6)
    double precision::G2_SO(12,12,12,12)

    integer::p,q,r,s

    G2_SO=0.d0
    do p=1,12
        do r=1,12
            if(mod(p,2) .ne.mod(r,2))continue
            do q=1,12
                do s=1,12
                    if(mod(q,2) .ne.mod(s,2))continue
                    G2_SO(p,q,r,s)=G2_SO(p,q,r,s)+G2((p+1)/2,(r+1)/2,(q+1)/2,(s+1)/2)
                end do 
            end do 
        end do 
    end do 
end subroutine another_Spatial2Spin

subroutine  Joshua_Spatial2Spin(G2,G2_SO)
    implicit none
    double precision::G2(6,6,6,6)
    double precision::G2_SO(12,12,12,12)
    integer::p,q,r,s
    integer::spins(4)
    double precision::value1,value2

    do p=1,12 
        do q=1,12 
            do r=1,12 
                do s=1,12 
                    if(mod(p,2)==mod(r,2))then
                        spins(1)=1
                    else
                        spins(1)=0
                    end if 

                    if(mod(q,2)==mod(s,2))then
                        spins(2)=1
                    else
                        spins(2)=0
                    end if 

                    if(mod(p,2)==mod(s,2))then
                        spins(3)=1
                    else
                        spins(3)=0
                    end if

                    if(mod(q,2)==mod(r,2))then
                        spins(4)=1
                    else
                        spins(4)=0
                    end if 

                    value1 = G2((p+1)/2,(r+1)/2,(q+1)/2,(s+1)/2) * spins(1)*spins(2)
                    value2 = G2((p+1)/2,(s+1)/2,(q+1)/2,(r+1)/2) * spins(3)*spins(4)
                    G2_SO(p,q,r,s)=value1-value2
                end do 
            end do
        end do 
    end do 
end subroutine Joshua_Spatial2Spin

!subroutine Get_DoubleExcitation

!end subroutine Get_DoubleExcitation

!subroutine Update_Intermediates
!! Stanton equation's #3
!do a=5,12
!    do e=5,12
!        Fae(a,e)=(1-KroDelta(a,e))*Fs(a,e)
!        do m=1,4
!            do f=5,12
!                do n=1,4
!                    Fae(a,e)=Fae(a,e)-0.5d0*taus
!                end do 
!            end do 
!        end do 
!    end do 
!end do 
!
!
!end subroutine Update_Intermediates
!
!double function tau(a,b,i,j)
!    integer::i,j,a,b
!    tau()
!end function tau
!
!integer function  KroDelta(i,j)
!    integer::i,j
!    KroDelta=0
!    if(i==j)KroDelta=1
!end function KroDelta
    

