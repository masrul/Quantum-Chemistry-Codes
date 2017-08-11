Module CCD
    use HF
	double precision::G2_MO(nO,nO,nO,nO)
    double precision::EMP2
    double precision::G2_SO(nSO,nSO,nSO,nSO)
    double precision::tau(nSO,nSO,nSO,nSO)
    double precision::ESO(nSO)
    double precision::Fae(nSO,nSO)
    double precision::Fmi(nSO,nSO)
    double precision::Fme(nSO,nSO)
    double precision::Wmnij(nSO,nSO,nSO,nSO)
    double precision::Wabef(nSO,nSO,nSO,nSO)
    double precision::Wmbej(nSO,nSO,nSO,nSO)
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


subroutine Get_EMP2
    implicit none 
    integer::i,j,a,b
    double precision::numerator,denominator
    EMP2=0.0
    do i=1,nElec/2
        do j=1,nElec/2
            do a=nElec/2+1,nO
                do b=nElec/2+1,nO
                    numerator=G2_MO(i,j,a,b)*(2*G2_MO(i,j,a,b)-G2_MO(i,j,b,a))
                    denominator=MOEnergy(i)+MOEnergy(j)-MOEnergy(a)-MOEnergy(b)
                    EMP2=EMP2+numerator/denominator
                end do 
            end do 
        end do 
    end do 
    print*,EMP2 
end subroutine Get_EMP2


subroutine Convert_Spatial2Spin(G2,G2_SO)
    double precision::G2(nO,nO,nO,nO)
    double precision::G2_SO(nSO,nSO,nSO,nSO)

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
end subroutine Convert_Spatial2Spin

subroutine  Get_initTau
    integer::i,j,a,b
    double precision::numerator,denominator
    tau=0.d0
    do i=1,nElec
        do j=1,nElec
            do a=nElec+1,nSO
                do b=nElec+1,nSO
                    numerator=G2_SO(i,j,a,b)
                    denominator=ESO(i)+ESO(j)-ESO(a)-ESO(b)
                    tau(i,j,a,b)=tau(i,j,a,b)+numerator/denominator
                end do 
            end do 
        end do 
    end do 
end subroutine Get_initTau


subroutine Convert_EMO_2_ESO
    integer::i,alpha,beta
    do i=1,nO
        alpha=2*i-1
        beta=alpha+1
        ESO(alpha)=MOEnergy(i)
        ESO(beta)=MOEnergy(i)
    end do
end subroutine Convert_EMO_2_ESO

subroutine Get_Fae   !#eq-3
    integer::a,e,m,f,n
    do a=nElec+1,nSO
        do e=nElec+1,nSO
            Fae(a,e)=(1-KroDelta(a,e))*ESO(a)
            do m=1,nElec
                do f=nElec+1,nSO
                    do n=1,nElec
                        Fae(a,e)=Fae(a,e)-0.5d0*tau(m,n,a,f)*G2_SO(m,n,e,f)
                    end do 
                end do 
            end do 
        end do 
    end do 
end subroutine Get_Fae

subroutine Get_Fmi   !#eq-4
    integer::m,i,e,n,f

    do m=1,nElec
        do i=1,nElec
            Fmi(m,i)=(1-KroDelta(m,i))*ESO(m)
            do e=1,nElec
                do n=nElec+1,nSO
                    do f=1,nElec
                        Fmi(m,i)=Fmi(m,i)-0.5d0*tau(i,n,e,f)*G2_SO(m,n,e,f)
                    end do 
                end do 
            end do 
        end do 
    end do 
end subroutine Get_Fmi

subroutine Get_Fme   !#eq-5
    implicit none 
    integer::m,e 
    Fme=0.d0
    do m=nElec+1,nSO
        do e=nElec+1,nSO
            Fme(m,e)=KroDelta(m,e)*ESO(m)
        end do 
    end do 
end subroutine Get_Fme


subroutine Get_Wmnij !#eq-6
    implicit none 
    integer::m,n,i,j,e,f
    do m=1,nElec
        do n=1,nElec
            do i=1,nElec
                do j=1,nElec
                    Wmnij(m,n,i,j)=G2_SO(m,n,i,j)
                    do e=nElec+1,nSO
                        do f=nElec+1,nSO
                            Wmnij(m,n,i,j)=Wmnij(m,n,i,j)+.25d0*tau(i,j,e,f)*G2_SO(m,n,e,f)
                        end do 
                    end do 
                end do 
            end do 
        end do 
    end do 

end subroutine Get_Wmnij

subroutine Get_Wabef !#eq-7
    implicit none 
    integer::m,n,e,f,a,b
    do a=nElec+1,nSO
        do b=nElec+1,nSO
            do e=nElec+1,nSO
                do f=nElec+1,nSO
                    Wabef(a,b,e,f)=G2_SO(a,b,e,f)
                    do m=1,nElec
                        do n=1,nElec
                            Wabef(a,b,e,f)=Wabef(a,b,e,f)+0.25d0*tau(m,n,a,b)*G2_SO(m,n,e,f)
                        end do 
                    end do 
                end do 
            end do 
        end do 
    end do 
end subroutine Get_Wabef

subroutine Get_Wmbej !#eq-8
    implicit none 
    integer::m,n,e,j,f,b
    do m=1,nElec
        do b=nElec+1,nSO
            do e=nElec+1,nSO
                do j=1,nElec
                    Wmbej(m,b,e,j)=G2_SO(m,b,e,j)
                    do n=nElec+1,nSO
                        do f=1,nElec
                            Wmbej(m,b,e,j)=Wmbej(m,b,e,j)-0.5d0*tau(j,n,f,b)*G2_SO(m,n,e,f)
                        end do 
                    end do 
                end do 
            end do 
        end do 
    end do 
end subroutine Get_Wmbej

subroutine Get_DoubleExcitation
    do i=1,nElec
        do j=1,nElec
            do a=nElec+1,nSO
                do b=nElec+1,nSO
                    tdnew(i,j,a,b)=G2_SO(i,j,a,b)
                    
                    do e=nElec+1,nSO
                        tdnew(i,j,a,b)=t(i,j,a,b)*




        
end subroutine Get_DoubleExcitation

integer function  KroDelta(i,j)
    integer::i,j
    KroDelta=0
    if(i==j)KroDelta=1
end function KroDelta

End Module CCD
