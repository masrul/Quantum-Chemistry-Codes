Module HartreeFock
    integer,parameter::K=6,N=4
    double precision,parameter::EnergyTol=10D-8
    double precision,parameter::DensityTol=10D-8
	double precision::S(K,K),H(K,K),G2(K,K,K,K),G_MO(K,K,K,K),GD2(K,K)
	double precision::X(K,K),H_MO(K,K),H_SMO(2*K,2*K),F_SMO(12,12),G_SMO(2*K,2*K,2*K,2*K)
	double precision::G(K,K),F(K,K),C(K,K)
    double precision::D(K,K),pC(K,K)
	double precision::newD(K,K),pF(K,K),D_SMO(2*K,2*K)
    double precision::MOEnergy(K)
    double precision::MBPT2Energy
	integer::iSCF
    double precision::E,newE
    logical::EnergyConverged,DensityConverged
    logical::Converged

    double precision::rmsD,errE

    contains
    
    subroutine  readParams()
        implicit none 
        OPEN (4, FILE='MULBERSHG.DAT', FORM='FORMATTED')
		Rewind 4
      	READ (4,*) S, H, G2
      	CLOSE (4)


        OPEN (4, FILE='BERSHG.DAT', FORM='FORMATTED')
		Rewind 4
      	READ (4,*) S, H, GD2
      	CLOSE (4)
    end subroutine  readParams

	subroutine Diagonalize_Matrix(Matrix,EigenValues,len)
        implicit none
	    integer len,l,inf
	    double precision  Matrix(len,len),EigenValues(len),work(len*(3+len/2))
	    l=len*(3+len/2)
	    call dsyev('V','U',len,Matrix,len,EigenValues,work,l,inf)
  	end subroutine Diagonalize_Matrix

    subroutine Get_XMatrix()
        double precision::DiagInvSmallS(K,K),eig(K)
        integer::i,j
        call Diagonalize_Matrix(S,eig,K)
		DiagInvSmallS=0.d0
		do i=1,K
			DiagInvSmallS(i,i)=1.d0/dsqrt(eig(i))
		end do 
        X=matmul(matmul(S,DiagInvSmallS),transpose(S))
    end subroutine Get_XMatrix

    subroutine Get_transformedCMatrix()
        call Diagonalize_Matrix(pF,MOEnergy,K)
        pC=pF
    end subroutine Get_transformedCMatrix

    
    subroutine Get_GMatrix()
        implicit none 
        integer::mu,nu,lamda,sigma
        G=0.0d0
        do mu=1,K
            do nu=1,K
                do lamda=1,K
                    do sigma=1,K
                        G(mu,nu)=G(mu,nu)+&
                        D(lamda,sigma)*(G2(mu,nu,sigma,lamda)-0.5d0*&
                        G2(mu,lamda,sigma,nu))
                    end do 
                end do 
            end do 
        end do 
    end subroutine Get_GMatrix

    subroutine Get_FockMatrix()
        implicit none 
        integer::i,j

        do i=1,K
            do j=1,K
                F(i,j)=H(i,j)+G(i,j)
            end do 
        end do 
    end subroutine Get_FockMatrix

    subroutine Get_transformedFockMatrix()
	    pF=matmul(matmul(transpose(X),F),X)
    end subroutine Get_transformedFockMatrix


    subroutine Get_CMatrix()
        C=MatMul(X,pC)
    end subroutine Get_CMatrix

    subroutine Get_newDensityMatrix()
        implicit none 
        integer::mu,nu,a
        newD=0.d0

        do mu=1,K
            do nu=1,K
                do a=1,N/2
                    newD(mu,nu)=newD(mu,nu)+C(mu,a)*C(nu,a)
                end do 
                newD(mu,nu)=2*newD(mu,nu)
            end do 
        end do 
    end subroutine Get_newDensityMatrix

    subroutine Init_HF()
        E=0.d0
	    D=0.d0
    end subroutine Init_HF

    subroutine Get_newEnergy()
        integer::mu,nu

        newE=0
        do mu=1,K
            do nu=1,K
                newE=newE+newD(nu,mu)*(H(mu,nu)+F(mu,nu))
            end do 
        end do 
        newE=0.5d0*newE
    end subroutine Get_newEnergy

    subroutine Get_Convergence()
        integer::i,j
        EnergyConverged=.False.
        DensityConverged=.False.
        Converged=.False.
        
        errE=abs(E-newE)
        if(errE<EnergyTol)EnergyConverged=.True.

        rmsD=0.d0
        do i=1,K
            do j=1,K
                rmsD=rmsD+(D(i,j)-newD(i,j))**2
            end do 
        end do 
        rmsD=dsqrt(rmsD/K**2)
        if(rmsD<DensityTol)DensityConverged=.True.

        if(EnergyConverged .and. DensityConverged)Converged=.True.
    end subroutine Get_Convergence

    subroutine PrintResults()

        if(Converged)then
             print*,'**********  Congratulation!!! SCF has Converged  ***********'
             print*,'------------------------------------------------------------'
             print*
        end if 
        print*,'SCF CYCLE NO            :   ',iSCF
        print*,'rmsDensityError         :   ',rmsD
        print*,'DENSITY CONVERGENCE     :   ',DensityConverged
        print*,'absErrorEnergy          :   ',errE
        print*,'ENERGY CONVERGENCE      :   ',EnergyConverged
        print*,'HF ENERGY               :   ',newE
        print*
        print*
    end subroutine PrintResults

    subroutine Get_HFEnergy
        implicit none 
        call readParams()
        call Get_XMatrix()
        call Init_HF()

        iSCF=0
        SCF:do 
            iSCF=iSCF+1
            call Get_GMatrix()
            call Get_FockMatrix()
            call Get_transformedFockMatrix()
            call Get_transformedCMatrix()
            call Get_CMatrix()
            call Get_newDensityMatrix()
            call Get_newEnergy()
            call Get_Convergence()
            !call PrintResults()
            if(Converged)exit SCF
            D=newD
            E=newE
        end do SCF
    end subroutine Get_HFEnergy


    subroutine Convert_AO2MO
        implicit none 
        integer::i,j,a,b
        integer::mu,nu,lamda,sigma
        G2_MO=0.0
        do i=1,K
            do j=1,K
                do a=1,K
                    do b=1,K
                        do  mu=1,K
                            do nu=1,K
                                do lamda=1,K
                                    do sigma=1,K
                                        G2_MO(i,a,j,b)=G2_MO(i,a,j,b)+&
                                        &C(mu,i)*C(nu,a)*C(lamda,j)*C(sigma,b)*G2(mu,nu,lamda,sigma)
                                    end do
                                end do 
                            end do 
                        end do 
                    end do 
                end do 
            end do 
        end do 
        
    end subroutine Convert_AO2MO
    

    subroutine Get_MBPT2Energy
        implicit none 
        integer::i,j,a,b
        double precision::numerator,denominator
        MBPT2Energy=0.0
        do i=1,2
            do j=1,2
                do a=3,6
                    do b=3,6
                        numerator=G2_MO(i,a,j,b)*(2*G2_MO(i,a,j,b)-G2_MO(i,b,j,a))
                        denominator=MOEnergy(i)+MOEnergy(j)-MOEnergy(a)-MOEnergy(b)
                        !print*,denominator
                        MBPT2Energy=MBPT2Energy+numerator/denominator
                    end do 
                end do 
            end do 
        end do
        print*,MBPT2Energy 
    end subroutine Get_MBPT2Energy
	
!    subroutine Get_EMP2
!		implicit none 
!	    integer::i,j,a,b
!		double precision::fs(DK)
!	    EMP2=0.d0
!		fs=0.d0
!		j=0
! 		do i=1,K
!			j=j+1
!			fs(j)=MOEnergy(i)
!			fs(j+1)=MOEnergy(i)
!			j=j+1
!		end do 	
!		print*,fs
    end subroutine Get_MBPT2Energy
    
    subroutine Convert_H6_2_H12
        H_SMO=0.d0
        H_SMO(1:K,1:K)=H_MO
        H_SMO(K+1:2*K,K+1:2*K)=H_MO
    end subroutine Convert_H6_2_H12


    subroutine Convert_G6_2_G12
        G_SMO=0.d0
        !ASSOCIATE(a=>(1:K),b=>(K+1:2*K))
            G_SMO(1:k,1:k,1:k,1:k)=G_MO
            G_SMO(k+1:2*k,k+1:2*k,k+1:2*k,K+1:2*k)=G_MO
            G_SMO(1:k,k+1:2*k,1:k,k+1:2*k)=G_MO
            G_SMO(k+1:2*k,1:k,k+1:2*k,1:k)=G_MO
        !END ASSOCIATE
    end subroutine Convert_G6_2_G12

    subroutine Convert_P6_2_P12
        D_SMO=0.d0
        D_SMO(1,1)=1.d0
        D_SMO(2,2)=1.d0
        D_SMO(7,7)=1.d0
        D_SMO(8,8)=1.d0
    end subroutine Convert_P6_2_P12

    subroutine Get_F_SMO
        do mu=1,2*K
            do nu=1,2*K
                do lamda=1,2*K
                    do sigma=1,2*K
                        F_SMO(mu,nu)=F_SMO(mu,nu)+H_SMO(mu,nu)+&
                        &D_SMO(lamda,sigma)*(G_SMO(mu,nu,sigma,lamda)-0.5d0*&
                        &G_SMO(mu,lamda,sigma,nu))
                    end do
                end do 
            end do 
        end do 
      do mu=1,12
           Print*,F_SMO(mu,1),
      end do

    end subroutine Get_F_SMO





End Module HartreeFock
