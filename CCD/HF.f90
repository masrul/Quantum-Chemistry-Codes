Module HF
    implicit none 
    integer,parameter::nO=6,nElec=4,nSO=12
    double precision,parameter::EnergyTol=10D-8
    double precision,parameter::DensityTol=10D-8
    double precision::S(nO,nO),H(nO,nO),G2(nO,nO,nO,nO)
    double precision::X(nO,nO)
    double precision::G(nO,nO),F(nO,nO),C(nO,nO)
    double precision::D(nO,nO),pC(nO,nO)
    double precision::newD(nO,nO),pF(nO,nO)
    integer::iSCF
    double precision::MOEnergy(nO)
    double precision::EHF,newEHF
    logical::EnergyConverged,DensityConverged
    logical::Converged
    double precision::rmsD,errE

contains

subroutine  readParams()
    implicit none 
    OPEN (4, FILE='BERSHG.DAT', FORM='FORMATTED')
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
    double precision::DiagInvSmallS(nO,nO),eig(nO)
    integer::i,j
    call Diagonalize_Matrix(S,eig,nO)
    DiagInvSmallS=0.d0
    do i=1,nO
        DiagInvSmallS(i,i)=1.d0/dsqrt(eig(i))
    end do 
    X=matmul(matmul(S,DiagInvSmallS),transpose(S))
end subroutine Get_XMatrix


subroutine Get_transformedCMatrix()
    call Diagonalize_Matrix(pF,MOEnergy,nO)
    pC=pF
end subroutine Get_transformedCMatrix

subroutine Get_GMatrix()
    implicit none 
    integer::mu,nu,lamda,sigma
    G=0.0d0
    do mu=1,nO
        do nu=1,nO
            do lamda=1,nO
                do sigma=1,nO
                    G(mu,nu)=G(mu,nu)+&
                    D(lamda,sigma)*(G2(mu,sigma,nu,lamda)-0.5d0*&
                    G2(mu,sigma,lamda,nu))
                end do 
            end do 
        end do 
    end do 
end subroutine Get_GMatrix

subroutine Get_FockMatrix()
    implicit none 
    integer::i,j

    do i=1,nO
        do j=1,nO
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

    do mu=1,nO
        do nu=1,nO
            do a=1,nElec/2
                newD(mu,nu)=newD(mu,nu)+C(mu,a)*C(nu,a)
            end do 
            newD(mu,nu)=2*newD(mu,nu)
        end do 
    end do 
end subroutine Get_newDensityMatrix

subroutine Init_HF()
    EHF=0.d0
    D=0.d0
end subroutine Init_HF

subroutine Get_newEHFnergy()
    integer::mu,nu

    newEHF=0
    do mu=1,nO
        do nu=1,nO
            newEHF=newEHF+newD(nu,mu)*(H(mu,nu)+F(mu,nu))
        end do 
    end do 
    newEHF=0.5d0*newEHF
end subroutine Get_newEHFnergy

subroutine Get_Convergence()
    integer::i,j
    EnergyConverged=.False.
    DensityConverged=.False.
    Converged=.False.
    
    errE=abs(EHF-newEHF)
    if(errE<EnergyTol)EnergyConverged=.True.

    rmsD=0.d0
    do i=1,nO
        do j=1,nO
            rmsD=rmsD+(D(i,j)-newD(i,j))**2
        end do 
    end do 
    rmsD=dsqrt(rmsD/nO**2)
    if(rmsD<DensityTol)DensityConverged=.True.

    if(EnergyConverged .and. DensityConverged)Converged=.True.
end subroutine Get_Convergence

subroutine Get_EHF
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
        call Get_newEHFnergy()
        call Get_Convergence()
        if(Converged)exit SCF
        D=newD
        EHF=newEHF
    end do SCF
end subroutine Get_EHF

End Module HF
