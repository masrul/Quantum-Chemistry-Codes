!****************************************************************************************
! By: Zachary Windom
! Date: 8/1/2017
! Description: Converts Molecular Orbital of Be to Molecular Spin Orbitals
!              The purpose of program is to calculate the CCD energy using
!              several subroutines corresponding to Stanton/Bartlett paper.
!              defined as 'intermediate'. Intermediates are used to iterate and
!              solve RHS of amplitude equation. Interation proceeds until
!              t(12,12,12,12) amplitude matrix has rms value less that 10E-07
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
DOUBLE PRECISION :: A(6,6), B(6,6), C(6,6), primeF(6,6), diaPrimeF(6)
DOUBLE PRECISION :: primeC(6,6), Gmatrix(6,6), oldE, newE, diff,energy_den,num
DOUBLE PRECISION :: oldP(6,6), newP(6,6),reg(6,6), rootmeansq,trans(6,6,6,6)
LOGICAL :: success
integer::q,r,ss,aa,bb,ff,e
integer::mu,nu,lamba,sigma
INTEGER, PARAMETER :: out_unit=20
DOUBLE PRECISION :: energy_SO(12,12),G2_so(12,12,12,12),td(12,12,12,12)
DOUBLE PRECISION :: Fae(12,12),Fmi(12,12),Fme(12,12),tdnew(12,12,12,12)
DOUBLE PRECISION :: Wmnij(12,12,12,12),Wabef(12,12,12,12),Wmbej(12,12,12,12)
DOUBLE PRECISION :: ccd_energy,rmsold,rms,Dijab
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
            write(*,*) test(i,j)
        else
            test(i,j)=0.0
        endif
    END DO
END DO
CALL PRINT_MATRIX('s^-1/2',6,6,test,6)
!**********************************************
! STEP 3 ::::: Determines transformation matrix X = U*SU
!              from newly formed eigenvectors matrix, 'S',
!              and diagonal, inverted, sq root matrix composed of eigenvalues
!              of S, s^-1/2
reg=MATMUL(S,test)
reg=MATMUL(reg,transpose(S))
CALL PRINT_MATRIX('X MATRIX:::',6,6,reg,6)
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
    CALL PRINT_MATRIX('G MATRIX::', 6,6, Gmatrix,6)
    CALL PRINT_MATRIX('FOCK matrix',6,6,initF,6)                        
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
    CALL PRINT_MATRIX('PRIME FOCK MATRIX',6,6,primeF,6)
!**********************************************
! STEP 6 ::::: Diagonalize primeF to get primeC
!              IE eigenvectors of primeF == primeC
    lwork=-1
    call DSYEV('V','U',6,primeF,6,diaPrimeF,dummy2,lwork,ok)
    lwork=min(LWMAX,int(dummy2(1)))
    call DSYEV('V','U',6,primeF,6,diaPrimeF,dummy2,lwork,ok)
    CALL PRINT_MATRIX('EIGENVALUES OFPRIME F', 1,6, diaPrimeF,1)
    CALL PRINT_MATRIX('EIGENVectors OF PRIME F', 6,6, primeF,6)
    do i=1,6
        do j=1,6
            primeC(i,j)=primeF(i,j)
        end do
    end do
!**********************************************
! STEP 7 ::::: Find C matrix by C= s^-1/2 * primeC
    B = MATMUL(reg, primeC)
    CALL PRINT_MATRIX('real C matrix', 6,6, B,6)
     
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
    CALL PRINT_MATRIX('P MATRIX',6,6,newP,6)
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
    write(*,*) 'Iteration number is:  ', counter
    write(*,*) 'energy is: ', newE
    diff = abs(newE-oldE)        
END DO
write(*,*)
write(*,*) 'Final converged energy is:'
write(*,*) '**************************'
write(*,*) newE
write(*,*) 'Energy converged on the',counter,' iteration.'
write(*,*)
!**********************************************
! STEP 9 ::::: Transform AO matrix G2 to MO matrix
!              via quarter transformation
!              (6,6,6,6)=>(6,6,6,6)
!**********************************************
! ***** WHY IS QUARTER TRANSFORMATION FAILING? **********
!quarter_trans=0.d0
!do p=1,6
!  do nu=1,6
!    do lamba=1,6
!      do sigma=1,6
!        summed=0
!        do mu=1,6
!          summed=summed+B(p,mu)*G2(mu,nu,lamba, sigma)
!        end do
!        quarter_trans(p,nu,lamba,sigma)=summed
!      end do
!    end do
!  end do
!end do

!do p=1,6
!  do q=1,6
!    do lamba=1,6
!      do sigma=1,6
!        summed=0
!        do nu=1,6
!          summed=summed+B(q,nu)*quarter_trans(p,nu,lamba,sigma)
!        end do
!        quarter_trans(p,q,lamba,sigma)=summed
!      end do
!    end do
!  end do
!end do

!do p=1,6
!  do q=1,6
!    do r=1,6
!      do sigma=1,6
!        summed=0
!        do lamba=1,6
!          summed=summed+B(r,lamba)*quarter_trans(p,q,lamba,sigma)
!        end do
!        quarter_trans(p,q,r,sigma)=summed
!      end do
!    end do
!  end do
!end do

!do p=1,6
!  do q=1,6
!    do r=1,6
!      do ss=1,6
!        summed=0
!        do sigma=1,6
!          summed=summed+B(ss,sigma)*quarter_trans(p,q,r,sigma)
!        end do
!        quarter_trans(p,q,r,ss)=summed
!      end do
!    end do
!  end do
!end do
DO i=1,6
  DO j=1,6
    DO aa=1,6
      DO bb=1,6
        do mu=1,6
          do nu=1,6
            do lamba=1,6
              do sigma=1,6
               quarter_trans(i,j,aa,bb)=quarter_trans(i,j,aa,bb)+B(mu,i)*B(nu,aa)*G2(mu,lamba,nu,sigma)*B(lamba,j)*B(sigma,bb)
              end do
            end do
          end do
         end do
       end do
     end do
   end do
end do
total=0.0
do i=1,2
  do j=1,2
    do aa=3,6
      do bb=3,6
        num=quarter_trans(i,j,aa,bb)*(2*quarter_trans(i,j,aa,bb)-quarter_trans(i,j,bb,aa))
        energy_den=diaPrimeF(i)+diaPrimeF(j)-diaPrimeF(aa)-diaPrimeF(bb)
        total=total+(num/energy_den)
      end do
    end do
  end do
end do
write(*,*)'!!!! MP2 CORRELATION ENERGY IS: !!!!! ', total
write(*,*)'MBPT(2) energy is: ', total+newE
!**********************************************
!**********************************************
! create diagonal matrix of orbital energies from HF
call diagonal_SOfock(diaPrimeF,energy_SO)
! convert the 2e- MO integrals to MSO integrals
call MOtoMSO(quarter_trans,G2_so)
! calculate the inital t amplitude
! (can be check by plugging into energy expression to recover MBPT(2))
! Energy
call initAmp(G2_so,energy_SO,td)
! Obtain Correlation energy given the inital t amplitude
call get_ECCD(G2_so,td, ccd_energy)
print*,"ccd energy is: ",ccd_energy
counter=0
rmsold=0.0
do while (.True.)
  counter=counter+1
  call get_ECCD(G2_so,td,ccd_energy)
  ! Obtain the intermediates defined in Stanton paper
  call get_Fae(td,energy_SO,G2_so,Fae)
  call get_Fmi(td, energy_SO, G2_so,Fmi)
  call get_Fme(energy_SO,Fme)
  call get_Wmnij(td,G2_so,Wmnij) 
  call get_Wabef(td,G2_so,Wabef)
  call get_Wmbej(G2_so,Wmbej)
  tdnew=0.d0
  do i=1,4
    do j=1,4
      do aa=5,12
        do bb=5,12
          Dijab=energy_SO(i,i)+energy_SO(j,j)-energy_SO(aa,aa)-energy_SO(bb,bb)
          !1st term 
          tdnew(i,j,aa,bb)=G2_so(i,j,aa,bb)
          
          !term 2
          do e=5,12
            tdnew(i,j,aa,bb)=tdnew(i,j,aa,bb)+td(i,j,aa,e)*Fae(bb,e)-td(i,j,bb,e)*Fae(aa,e)
          end do
          !term 3
          do m=1,4
            tdnew(i,j,aa,bb)=tdnew(i,j,aa,bb)+td(i,m,aa,bb)*Fmi(m,j)-td(j,m,aa,bb)*Fmi(m,i)
          end do
          !term 4
          do e=5,12
            do ff=5,12
              tdnew(i,j,aa,bb)=tdnew(i,j,aa,bb)+0.5d0*td(i,j,e,ff)*Wabef(aa,bb,e,ff)
            end do
          end do
          !term 5
          do m=1,4
            do e=5,12
              tdnew(i,j,aa,bb)=tdnew(i,j,aa,bb)+td(i,m,aa,e)*Wmbej(m,bb,e,j)
              tdnew(i,j,aa,bb)=tdnew(i,j,aa,bb)-td(j,m,aa,e)*Wmbej(m,bb,e,i)
              tdnew(i,j,aa,bb)=tdnew(i,j,aa,bb)-td(i,m,bb,e)*Wmbej(m,aa,e,j)
              tdnew(i,j,aa,bb)=tdnew(i,j,aa,bb)+td(j,m,bb,e)*Wmbej(m,aa,e,i)
            end do
            do n=1,4
              tdnew(i,j,aa,bb)=tdnew(i,j,aa,bb)+0.5d0*td(i,j,aa,bb)*Wmnij(m,n,i,j)
            end do
         end do
         tdnew(i,j,aa,bb)=tdnew(i,j,aa,bb)/Dijab
        end do
      end do
    end do
  end do
  call amp_converge(td,tdnew,rms)
  print*, 'iteration',counter, 'with value: ', ccd_energy
  if (abs(rms-rmsold)<10E-07) then
    print*, 'converged at iteration',counter, 'with value: ', ccd_energy
    exit
  else
      td=tdnew
      continue
  end if
end do 

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
   write(*,*)
   write(*,*) '**** DENSITY MATRIX, P, HAS CONVERGED? ****'
   write(*,*) success
   RETURN 
END SUBROUTINE ROOT_MEAN_SQUARE  

SUBROUTINE MOtoMSO(quarter_trans,G2_so)
  DOUBLE PRECISION,INTENT(IN) :: quarter_trans(6,6,6,6)
  DOUBLE PRECISION,INTENT(OUT) :: G2_so(12,12,12,12)
  INTEGER :: i,j,k,l
  INTEGER :: alpha_i, alpha_j, alpha_k, alpha_l
  INTEGER :: beta_i, beta_j, beta_k, beta_l
  G2_so=0.d0
  do i=1,6
    alpha_i=2*i-1
    beta_i=alpha_i+1
    do j=1,6
      alpha_j=2*j-1
      beta_j=alpha_j+1
      do k=1,6
        alpha_k=2*k-1
        beta_k=alpha_k+1
        do l=1,6
          alpha_l=2*l-1
          beta_l=alpha_l+1
          G2_so(alpha_i,alpha_j,alpha_k,alpha_l)=quarter_trans(i,j,k,l)-quarter_trans(i,j,l,k)
          G2_so(alpha_i,beta_j,alpha_k,beta_l)=quarter_trans(i,j,k,l)
          G2_so(alpha_i, beta_j, beta_k, alpha_l)=-quarter_trans(i,j,l,k)
          G2_so(beta_i, alpha_j, alpha_k, beta_l)=-quarter_trans(i,j,l,k)
          G2_so(beta_i,alpha_j,beta_k,alpha_l)=quarter_trans(i,j,k,l)
          G2_so(beta_i,beta_j,beta_k,beta_l)=quarter_trans(i,j,k,l)-quarter_trans(i,j,l,k)
        end do
      end do
    end do
  end do
END SUBROUTINE MOtoMSO

! Purpose: Place MSO energies along diagonal Matrix, found from HF
SUBROUTINE diagonal_SOfock(diaPrimeF,energy_SO)
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(INOUT)::energy_SO(12,12)
  DOUBLE PRECISION, INTENT(IN)::diaPrimeF(6)
  INTEGER:: i,alpha,beta
  energy_SO=0.d0
  do i=1,6
    alpha=2*i-1
    beta=alpha+1
    energy_SO(alpha,alpha)=diaPrimeF(i)
    energy_SO(beta,beta)=diaPrimeF(i)
  end do
END SUBROUTINE diagonal_SOfock

! Purpose: Determines the initial T2 amplitude
SUBROUTINE initAmp(G2_so,energy_SO,td)
  INTEGER::i,j,a,b
  DOUBLE PRECISION, INTENT(IN)::G2_so(12,12,12,12),energy_SO(12,12)
  DOUBLE PRECISION::num,den
  DOUBLE PRECISION, INTENT(OUT)::td(12,12,12,12)
  td=0.d0
  do i=1,4
    do j=1,4
      do a=5,12
        do b=5,12
          num=G2_so(i,j,a,b)
          den=energy_SO(i,i)+energy_SO(j,j)-energy_SO(a,a)-energy_SO(b,b)
          td(i,j,a,b)=td(i,j,a,b)+num/den
        end do
      end do
    end do
  end do
END SUBROUTINE initAmp

! Stanton/Bartlett eq 3
SUBROUTINE get_Fae(td,energy_SO,G2_so,Fae)
  double precision,INTENT(OUT) :: Fae(12,12)
  DOUBLE PRECISION, INTENT(IN)::td(12,12,12,12),energy_SO(12,12),G2_so(12,12,12,12)
  integer::a,e,m,n,f,delta
  do a=5,12
    do e=5,12
      if (a==e) then
        Fae(a,e)=0
      else
          Fae(a,e)=energy_SO(a,e)
      end if
      do m=1,4
        do f=5,12
          do n=1,4
            Fae(a,e)=Fae(a,e)-0.5d0*td(m,n,a,f)*G2_so(m,n,e,f)
          end do
        end do
      end do
    end do
  end do
END SUBROUTINE get_Fae

! ! Stanton/Bartlett eq 4
SUBROUTINE get_Fmi(td, energy_SO, G2_so,Fmi)
  integer ::m,i,e,f,n
  double precision,INTENT(OUT) :: Fmi(12,12)
  DOUBLE PRECISION,INTENT(IN)::td(12,12,12,12),energy_SO(12,12)
  DOUBLE PRECISION, INTENT(IN)::G2_so(12,12,12,12)
  Fmi=0.d0
  do m=1,4
    do i=1,4
      if (m==i) then
        Fmi(m,i)=0
      else
          Fmi(m,i)=energy_SO(m,i)
      end if
      do e=5,12
        do n=1,4
          do f=5,12
            Fmi(m,i)=Fmi(m,i)-0.5d0*td(i,n,e,f)*G2_so(m,n,e,f)
          end do
        end do
      end do
    end do
  end do
END SUBROUTINE get_Fmi

! Stanton/Bartlett eq 5
SUBROUTINE get_Fme(energy_SO,Fme)

  INTEGER::m,e
  double precision,INTENT(OUT)::Fme(12,12)
  DOUBLE PRECISION, INTENT(IN)::energy_SO(12,12)
  Fme=0.d0
  do m=1,4
    do e=5,12
      Fme(m,e)=energy_SO(m,e)
    end do
  end do
END SUBROUTINE get_Fme

! Stanton/Bartlett eq 6
SUBROUTINE get_Wmnij(td,G2_so,Wmnij)
  INTEGER::m,n,i,j,e,f
  DOUBLE PRECISION,INTENT(OUT)::Wmnij(12,12,12,12)
  DOUBLE PRECISION, INTENT(IN)::td(12,12,12,12)
  DOUBLE PRECISION,INTENT(IN)::G2_so(12,12,12,12)
  do m=1,4
    do n=1,4
      do i=1,4
        do j=1,4
          Wmnij(m,n,i,j)=G2_so(m,n,i,j)
          do e=4,12
            do f=4,12
              Wmnij(m,n,i,j)=Wmnij(m,n,i,j)+0.25d0*td(i,j,e,f)*G2_so(m,n,e,f)
            end do
          end do
        end do
      end do
    end do
  end do
END SUBROUTINE get_Wmnij

! Stanton/Bartlett eq 7
SUBROUTINE get_Wabef(td,G2_so,Wabef)
  integer::m,n,e,f,a,b
  double precision,INTENT(OUT):: Wabef(12,12,12,12)
  DOUBLE PRECISION,INTENT(IN)::td(12,12,12,12),G2_so(12,12,12,12)
  do a=5,12
    do b=5,12
      do e=5,12
        do f=5,12
          Wabef(a,b,e,f)=G2_so(a,b,e,f)
          do m=1,4
            do n=1,4
              Wabef(a,b,e,f)=Wabef(a,b,e,f)+0.25d0*td(m,n,a,b)*G2_so(m,n,e,f)
            end do
          end do
        end do
      end do
    end do
  end do
END SUBROUTINE get_Wabef

! Stanton/Bartlett eq 8
SUBROUTINE get_Wmbej(G2_so,Wmbej)
  integer::m,n,e,j,f,b
  DOUBLE PRECISION,INTENT(OUT)::Wmbej(12,12,12,12)
  DOUBLE PRECISION,INTENT(IN)::G2_so(12,12,12,12)
  do m=1,4
    do b=5,12
      do e=5,12
        do j=1,4
          Wmbej(m,b,e,j)=G2_so(m,b,e,j)
          do n=1,4
            do f=5,12
              Wmbej(m,b,e,j)=Wmbej(m,b,e,j)
            end do
          end do
        end do
      end do
    end do
  end do
END SUBROUTINE get_Wmbej

! Purpose: determines the correlation energy given T2 amplitudes and the SO 2 e-
!          matrix 
SUBROUTINE get_ECCD(G2_so,td,ccd_energy)
  INTEGER::i,j,a,b
  double precision, INTENT(OUT)::ccd_energy
  DOUBLE PRECISION,INTENT(IN)::G2_so(12,12,12,12),td(12,12,12,12)
  ccd_energy=0.d0
  do i=1,4
    do j=1,4
      do a=5,12
        do b=5,12
          ccd_energy=ccd_energy+G2_so(i,j,a,b)*td(i,j,a,b)
        end do
      end do
    end do
  end do
  ccd_energy=ccd_energy*0.25
END SUBROUTINE get_ECCD

! Purpose: Determines whether the RMS value of consecutive T amplitudes have
!          fallen beneath 10E-07
SUBROUTINE amp_converge(td,tdnew,rms)
  INTEGER::i,j,k,l
  DOUBLE PRECISION,INTENT(IN)::td(12,12,12,12),tdnew(12,12,12,12)
  DOUBLE PRECISION,intent(out)::rms
  rms=0.d0
  do i=1,12
    do j=1,12
      do k=1,12
        do l=1,12
          rms=rms+(td(i,j,k,l)-tdnew(i,j,k,l))**2
        end do
      end do
    end do
  end do
  rms=dsqrt(rms/12**4)
END SUBROUTINE amp_converge



