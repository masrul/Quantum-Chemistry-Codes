Module HartreeFock
    integer,parameter::K=6,N=4
    double precision,parameter::EnergyTol=10D-8
    double precision,parameter::DensityTol=10D-8
	double precision::S(K,K),H(K,K),G2(K,K,K,K)
	double precision::X(K,K)
	double precision::G(K,K),F(K,K),C(K,K)
    double precision::D(K,K),pC(K,K)
	double precision::newD(K,K),pF(K,K)
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
        double precision::eig(K)
        call Diagonalize_Matrix(pF,eig,K)
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
End Module HartreeFock
