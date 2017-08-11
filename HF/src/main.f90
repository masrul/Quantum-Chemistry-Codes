program main
	implicit none
    integer::K=6,N=4 
	double precision::S(6,6),H(6,6),G2(6,6,6,6)
	double precision::X(6,6)
	double precision::D(6,6),G(6,6),F(6,6),pF(6,6),C(6,6)
	double precision::newD(6,6)
	integer::i
    double precision::E,newE
    logical::EnergyConverged,DensityConverged

	call readParams()
    call Get_XMatrix()
    call Init_HF()
    i=0

    SCF:do 
        i=i+1
        call Get_GMatrix()
        call Get_FockMatrix()
        call Get_transformedFockMatrix()
        call Get_transformedCMatrix()
        call Get_CMatrix()
        call Get_newDensityMatrix()
        call Get_newEnergy()
        if(Converged())then
            print*,'**********  Congratulation SCF has Converged  ***********'
            print*,'---------------------------------------------------------'
            print*
            call PrintResults()
            exit SCF
        else
            call PrintResults()
        end if 
        D=newD
        E=newE
    end do SCF

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
        C=MatMul(X,pF)
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

    function Converged()
        logical::Converged
        double precision::rmsD
        integer::i,j
        EnergyConverged=.False.
        DensityConverged=.False.
        Converged=.False.

        if(abs(E-newE)<10e-9)EnergyConverged=.True.

        rmsD=0.d0
        do i=1,K
            do j=1,K
                rmsD=rmsD+(D(i,j)-newD(i,j))**2
            end do 
        end do 
        rmsD=dsqrt(rmsD)
        if(rmsD<10e-5)DensityConverged=.True.

        if(EnergyConverged .and. DensityConverged)Converged=.True.
    end function Converged

    subroutine PrintResults()
        print*,'SCF cycle no       :',i
        print*,'Density Convergence:',DensityConverged
        print*,'Energy Convergence :',EnergyConverged
        print*,'HF energy          :',newE
        print*
        print*
    end subroutine PrintResults
end program main
