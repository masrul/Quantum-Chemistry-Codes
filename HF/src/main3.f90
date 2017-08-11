program main
	implicit none
    integer,parameter::K=6,N=4
    double precision,parameter::EnergyTol=10e-10 
    double precision,parameter::DensityTol=10e-8
	double precision::S(K,K),H(K,K),G2(K,K,K,K)
	double precision::X(K,K)
	double precision::G(K,K),F(K,K),C(K,K)
    double precision,pointer::D(:,:),pC(:,:)
	double precision,target::newD(K,K),pF(K,K)
	integer::i
    double precision::E,newE
    logical::EnergyConverged,DensityConverged

	call readParams()
    call Get_XMatrix()
    call Init_HF()
    i=0

    SCF:do 
        i=i+1

        write(4,*) 'SCF CYCLE:',i
        call Get_GMatrix()
        write(4,*) 'GMatrix'
        call Print_Matrix(G)

        call Get_FockMatrix()
        write(4,*) 'FockMatrix'
        call Print_Matrix(F)

        call Get_transformedFockMatrix()
        write(4,*) 'pFMatrix'
        call Print_Matrix(pF)

        call Get_transformedCMatrix()
        write(4,*) 'pCMatrix'
        call Print_Matrix(pC)

        call Get_CMatrix()
        write(4,*) 'CMatrix'
        call Print_Matrix(C)

        call Get_newDensityMatrix()
        write(4,*) 'PMatrix'
        call Print_Matrix(newD)

        call Get_newEnergy()
        write(4,*) 'Energy:',newE
        write(4,*) 
        write(4,*) 

        if(Converged())then
            print*,'**********  Congratulation!!! SCF has Converged  ***********'
            print*,'------------------------------------------------------------'
            print*
            call PrintResults()
            exit SCF
        else
            call PrintResults()
        end if 
        D=>newD
        E=newE
    end do SCF

    contains
    
    subroutine  readParams()
        implicit none 
        OPEN (5, FILE='MULBERSHG.DAT', FORM='FORMATTED')
		Rewind 5
      	READ (5,*) S, H, G2
        call Print_Matrix(S)
        call Print_Matrix(H)
        call Print_Matrix(G2)
      	CLOSE (5)
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
        write(4,*)'DiagInvSmallS'
        call Print_Matrix(DiagInvSmallS) 
        X=matmul(matmul(S,DiagInvSmallS),transpose(S))
        write(4,*)'X Matrix'
        call Print_Matrix(X) 
    end subroutine Get_XMatrix

    subroutine Get_transformedCMatrix()
        double precision::eig(K)
        call Diagonalize_Matrix(pF,eig,K)
        pC=>pF
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
        allocate(D(K,K))
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

    function Converged()
        logical::Converged
        double precision::rmsD
        integer::i,j
        EnergyConverged=.False.
        DensityConverged=.False.
        Converged=.False.

        if(abs(E-newE)<EnergyTol)EnergyConverged=.True.

        rmsD=0.d0
        do i=1,K
            do j=1,K
                rmsD=rmsD+(D(i,j)-newD(i,j))**2
            end do 
        end do 
        rmsD=dsqrt(rmsD)
        if(rmsD<DensityTol)DensityConverged=.True.

        if(EnergyConverged .and. DensityConverged)Converged=.True.
    end function Converged

    subroutine PrintResults()
        print*,'SCF CYCLE NO            :   ',i
        print*,'DENSITY CONVERGENCE     :   ',DensityConverged
        print*,'ENERGY CONVERGENCE      :   ',EnergyConverged
        print*,'HF ENERGY               :   ',newE
        print*
        print*
    end subroutine PrintResults

    subroutine Print_Matrix(Matrix)
        double precision::Matrix(K,K)
        integer::i
        do i=1,K
            write(4,'(6F15.8)') Matrix(i,:)
        end do 
        write(4,*)

    end subroutine Print_Matrix
end program main
