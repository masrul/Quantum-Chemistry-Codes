Module CCD
    use HF
	double precision::G2_MO(nO,nO,nO,nO)
    double precision::EMP2,ECCD,OLDCCD
    double precision::G2_SO(nSO,nSO,nSO,nSO)
    double precision::td(nSO,nSO,nSO,nSO)
    double precision::tdnew(nSO,nSO,nSO,nSO)
    double precision::ESO(nSO,nSO)
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
                    td(a,b,i,j)=td(a,b,i,j)+numerator/denominator
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

subroutine Get_Fae   !#eq-3
    integer::a,e,m,f,n
    Fae=0.d0
    do a=nElec+1,nSO
        do e=nElec+1,nSO
            Fae(a,e)=(1-KroDelta(a,e))*ESO(a,e)
            do m=1,nElec
                do f=nElec+1,nSO
                    do n=1,nElec
                        Fae(a,e)=Fae(a,e)-0.5d0*td(a,f,m,n)*G2_SO(m,n,e,f)
                    end do 
                end do 
            end do 
        end do 
    end do
!    print*,'Fae' 
!    print*,Fae 
end subroutine Get_Fae

subroutine Get_Fmi   !#eq-4
    integer::m,i,e,n,f
    Fmi=0.d0
    do m=1,nElec
        do i=1,nElec
            Fmi(m,i)=(1-KroDelta(m,i))*ESO(m,i)
            do e=1,nElec
                do n=nElec+1,nSO
                    do f=1,nElec
                        Fmi(m,i)=Fmi(m,i)-0.5d0*td(e,f,i,n)*G2_SO(m,n,e,f)
                    end do 
                end do 
            end do 
        end do 
    end do 
!    print*,'Fmi'
!    print*,Fmi
end subroutine Get_Fmi

subroutine Get_Fme   !#eq-5
    implicit none 
    integer::m,e 
    Fme=0.d0
    do m=nElec+1,nSO
        do e=nElec+1,nSO
            Fme(m,e)=ESO(m,e)
        end do 
    end do 
!    print*,'Fme'
!    print*,Fme
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
                            Wmnij(m,n,i,j)=Wmnij(m,n,i,j)+.25d0*td(e,f,i,j)*G2_SO(m,n,e,f)
                        end do 
                    end do 
                end do 
            end do 
        end do 
    end do 
!    print*,'Wmnij'
!    print*,Wmnij
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
                            Wabef(a,b,e,f)=Wabef(a,b,e,f)+0.25d0*td(a,b,m,n)*G2_SO(m,n,e,f)
                        end do 
                    end do 
                end do 
            end do 
        end do 
    end do 
!    print*,'Wabef'
!    print*,Wabef
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
                            Wmbej(m,b,e,j)=Wmbej(m,b,e,j)-0.5d0*td(f,b,j,n)*G2_SO(m,n,e,f)
                        end do 
                    end do 
                end do 
            end do 
        end do 
    end do 
!    print*,'Wmbej'
!    print*,Wmbej
end subroutine Get_Wmbej

subroutine Update_DoubleExcitation
    implicit none 
    integer::i,j,a,b
    integer::e,m,n,f
    double precision:: Dijab
    tdnew=0.d0
    call Get_Fae
    call Get_Fmi
    call Get_Fme

    call Get_Wmnij
    call Get_Wabef
    call Get_Wmbej
    do i=1,nElec
        do j=1,nElec
            do a=nElec+1,nSO
                do b=nElec+1,nSO
                    Dijab=ESO(i,i)+ESO(j,j)-ESO(a,a)-ESO(b,b)

                    !*** First term ***!
                    tdnew(a,b,i,j)=G2_SO(i,j,a,b)

                    !*** 2nd Term ***!
                    do e=nElec+1,nSO
                        tdnew(a,b,i,j)=tdnew(a,b,i,j)+td(a,e,i,j)*Fae(b,e)-td(b,e,i,j)*Fae(a,e)
                    end do 
                    !*** 3rd Term ***!
                    do m=1,nElec 
                        tdnew(a,b,i,j)=tdnew(a,b,i,j)+td(a,b,i,m)*Fmi(m,j)-td(a,b,j,m)*Fmi(m,i)
                    end do 

                    !*** 4th Term ***!
                    do e=nElec+1,nSO
                        do f=nElec+1,nSO
                            tdnew(a,b,i,j)=tdnew(a,b,i,j)+0.5d0*td(e,f,i,j)*Wabef(a,b,e,f)
                        end do 
                    end do 

                    !*** 5th Term ***!
                    do m=1,nElec
                        do e=nElec+1,nSO
                            tdnew(a,b,i,j)=tdnew(a,b,i,j)+td(a,e,i,m)*Wmbej(m,b,e,j)
                            tdnew(a,b,i,j)=tdnew(a,b,i,j)-td(a,e,j,m)*Wmbej(m,b,e,i)
                            tdnew(a,b,i,j)=tdnew(a,b,i,j)-td(b,e,i,m)*Wmbej(m,a,e,j)
                            tdnew(a,b,i,j)=tdnew(a,b,i,j)+td(b,e,j,m)*Wmbej(m,a,e,i)
                        end do 
                    end do

                   !*** 6th Term ***! 
                   do n=1,nElec
                        tdnew(a,b,i,j)=tdnew(a,b,i,j)+0.5d0*td(a,b,i,j)*Wmnij(m,n,i,j)
                   end do 
                   tdnew(a,b,i,j)=tdnew(a,b,i,j)/Dijab
               end do 
            end do 
        end do 
    end do
    do i=1,12
        do j=1,12
            do a=1,12
                do b=1,12
                    if(td(i,j,a,b)>0.000000000001)then
                        !print*,td(i,j,a,b),tdnew(i,j,a,b)
                    end if 
                end do 
            end do
        end do 
    end do  


    td=tdnew
end subroutine Update_DoubleExcitation

integer function  KroDelta(i,j)
    integer::i,j
    KroDelta=0
    if(i==j)KroDelta=1
end function KroDelta

subroutine Get_ECCD
    implicit none 
    integer::i,j,a,b
    ECCD=0.d0
    do i=1,nElec
        do j=1,nElec
            do a=nElec+1,nSO
                do b=nElec+1,nSO
                    ECCD=ECCD+0.25*G2_SO(i,j,a,b)*td(a,b,i,j)
                end do 
            end do 
        end do 
    end do 
end subroutine Get_ECCD

End Module CCD
