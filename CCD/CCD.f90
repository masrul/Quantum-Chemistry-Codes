Module CCD
    use HF
    double precision::G2_MO(nO,nO,nO,nO)
    double precision::EMP2,ECCD,OLDCCD
    double precision::G2_SO(nSO,nSO,nSO,nSO)
    double precision::td(nSO,nSO,nSO,nSO)
    double precision::tdnew(nSO,nSO,nSO,nSO)
    double precision::ESO(nSO,nSO)
    double precision,parameter::rmsT2tol=1d-7
    logical::ccdConverged
    
contains
    
subroutine Convert_AOtoMO
    implicit none 
    integer::i,j,a,b
    integer::mu,nu,lamda,sigma
    double precision::muINT,nuINT,lamdaINT,sigmaINT
    G2_MO=0.0
    do i=1,nO
        do j=1,nO
            do a=1,nO
                do b=1,nO
                    sigmaINT=0.d0
                    do  sigma=1,nO
                        lamdaINT=0.d0
                        do lamda=1,nO
                            nuINT=0.d0
                            do nu=1,nO
                                muINT=0.d0
                                do mu=1,nO
                                    muINT=muINT+C(mu,i)*G2(mu,lamda,nu,sigma)
                                end do
                                nuINT=nuINT+C(nu,a)*muINT
                            end do 
                            lamdaINT=lamdaINT+C(lamda,j)*nuINT
                        end do 
                        sigmaINT=sigmaINT+C(sigma,b)*lamdaINT
                    end do
                    G2_MO(i,j,a,b)=G2_MO(i,j,a,b)+sigmaINT 
                end do 
            end do 
        end do 
    end do 
     
end subroutine Convert_AOtoMO

subroutine Convert_Spatial2Spin
    implicit none 
    integer::i,j,k,l
    integer::alpha1,beta1    
    integer::alpha2,beta2    
    integer::alpha3,beta3    
    integer::alpha4,beta4
    G2_SO=0.d0
    do i=1,nO
        alpha1=2*i-1
        beta1=alpha1+1
        do j=1,nO
            alpha2=2*j-1
            beta2=alpha2+1
            do k=1,nO
                alpha3=2*k-1
                beta3=alpha3+1
                do l=1,nO
                    alpha4=2*l-1
                    beta4=alpha4+1

                    G2_SO(alpha1,alpha2,alpha3,alpha4)=G2_MO(i,j,k,l)-G2_MO(i,j,l,k)
                    G2_SO(alpha1,beta2,alpha3,beta4)=G2_MO(i,j,k,l)
                    G2_SO(alpha1,beta2,beta3,alpha4)=-G2_MO(i,j,l,k)
                    G2_SO(beta1,alpha2,alpha3,beta4)=-G2_MO(i,j,l,k)
                    G2_SO(beta1,alpha2,beta3,alpha4)=G2_MO(i,j,k,l)
                    G2_SO(beta1,beta2,beta3,beta4)=G2_MO(i,j,k,l)-G2_MO(i,j,l,k)
                end do 
            end do 
        end do 
    end do 
end subroutine Convert_Spatial2Spin

subroutine  Init_CCD
    integer::i,j,a,b
    double precision::numerator,denominator
    td=0.d0
    do i=1,nElec
        do j=1,nElec
            do a=nElec+1,nSO
                do b=nElec+1,nSO
                    numerator=G2_SO(i,j,a,b)
                    denominator=ESO(i,i)+ESO(j,j)-ESO(a,a)-ESO(b,b)
                    td(i,j,a,b)=td(i,j,a,b)+numerator/denominator
                end do 
            end do 
        end do 
    end do 
end subroutine Init_CCD


subroutine Convert_EMO_2_ESO
    integer::i,alpha,beta
    do i=1,nO
        alpha=2*i-1
        beta=alpha+1
        ESO(alpha,alpha)=MOEnergy(i)
        ESO(beta,beta)=MOEnergy(i)
    end do
end subroutine Convert_EMO_2_ESO

subroutine Get_ECCD
    implicit none 
    integer::i,j,a,b
    ECCD=0.d0
    do i=1,nElec
        do j=1,nElec
            do a=nElec+1,nSO
                do b=nElec+1,nSO
                    ECCD=ECCD+0.25*G2_SO(i,j,a,b)*td(i,j,a,b)
                end do 
            end do 
        end do 
    end do 
end subroutine Get_ECCD

subroutine Check_rmsT2Convergence
    integer::i,j,k,l
    double precision::rmsT2
    rmsT2=0.d0
    ccdConverged=.False.
    do i=1,nSO
        do j=1,nSo
            do k=1,nSO
                do l=1,nSO
                    rmsT2=rmsT2+(td(i,j,k,l)-tdnew(i,j,k,l))**2
                end do 
            end do 
        end do 
    end do 
    rmsT2=dsqrt(rmsT2/nSO**4)
    if(rmsT2<rmsT2tol)ccdConverged=.True.
end subroutine Check_rmsT2Convergence


subroutine Update_tijab

    implicit none 
    integer::i,j,k,l
    integer::a,b,c,d
    double precision:: Dijab
 
    tdnew=0.d0
    do i=1,nElec
        do j=1,nElec
            do a=nElec+1,nSO
                do b=nElec+1,nSO
                    Dijab=ESO(i,i)+ESO(j,j)-ESO(a,a)-ESO(b,b)

                    ! 1st term
                    tdnew(i,j,a,b)=tdnew(i,j,a,b)+G2_SO(a,b,i,j)

                    ! 2nd term
                    do c=nElec+1,nSO
                        do d=nElec+1,nSO
                            tdnew(i,j,a,b)=tdnew(i,j,a,b)+0.5d0*G2_SO(a,b,c,d)*td(i,j,c,d)
                        end do 
                    end do 

                    ! 3rd term
                    do k=1,nElec
                        do l=1,nElec
                            tdnew(i,j,a,b)=tdnew(i,j,a,b)+0.5d0*G2_SO(i,j,k,l)*td(k,l,a,b)
                        end do 
                    end do 

                    ! 4th term 
                    do k=1,nElec
                        do c=nElec+1,nSO
                            tdnew(i,j,a,b)=tdnew(i,j,a,b)-G2_SO(b,k,c,j)*td(i,k,a,c)
                            tdnew(i,j,a,b)=tdnew(i,j,a,b)+G2_SO(b,k,c,i)*td(j,k,a,c)
                            tdnew(i,j,a,b)=tdnew(i,j,a,b)+G2_SO(a,k,c,j)*td(i,k,b,c)
                            tdnew(i,j,a,b)=tdnew(i,j,a,b)-G2_SO(a,k,c,i)*td(j,k,b,c)
                        end do 
                    end do 


                    ! 5th term
                    do k=1,nElec
                        do l=1,nElec
                            do c=nElec+1,nSO
                                do d=nElec+1,nSO
                                    tdnew(i,j,a,b)=tdnew(i,j,a,b)+0.25d0*G2_SO(k,l,c,d)*td(i,j,c,d)*td(k,l,a,b)

                                    tdnew(i,j,a,b)=tdnew(i,j,a,b)-0.50d0*G2_SO(k,l,c,d)*td(i,j,a,c)*td(k,l,b,d)
                                    tdnew(i,j,a,b)=tdnew(i,j,a,b)-0.50d0*G2_SO(k,l,c,d)*td(i,j,b,d)*td(k,l,a,c)

                                    tdnew(i,j,a,b)=tdnew(i,j,a,b)-0.50d0*G2_SO(k,l,c,d)*td(i,k,a,b)*td(j,l,c,d)
                                    tdnew(i,j,a,b)=tdnew(i,j,a,b)-0.50d0*G2_SO(k,l,c,d)*td(i,k,c,d)*td(j,l,a,b)

                                    tdnew(i,j,a,b)=tdnew(i,j,a,b)+G2_SO(k,l,c,d)*td(i,k,a,c)*td(j,l,b,d)
                                    tdnew(i,j,a,b)=tdnew(i,j,a,b)+G2_SO(k,l,c,d)*td(i,k,b,d)*td(j,l,a,c)

                                end do 
                            end do 
                        end do 
                    end do 

                    ! Divided by denominator
                    tdnew(i,j,a,b)=tdnew(i,j,a,b)/Dijab
                end do 
            end do 
        end do 
    end do 
    call Check_rmsT2Convergence
    td=tdnew
end subroutine Update_tijab

End Module CCD
