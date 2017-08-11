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

subroutine Convert_AOtoMO_2
    implicit none 
    integer::p,q,r,s
    integer::mu,nu,lamda,sigma
    double precision::muINT,nuINT,lamdaINT,sigmaINT
	double precision::G2_Inter(nO,nO,nO,nO)

    G2_MO=0.d0
    G2_Inter=0.d0
    do p=1,nO
        do nu=1,nO
            do lamda=1,nO
                do sigma=1,nO
                    do mu=1,nO
                        G2_Inter(p,nu,lamda,sigma)=G2_Inter(p,nu,lamda,sigma)+C(p,mu)*G2(mu,nu,lamda,sigma)
                    end do 
                end do 
            end do 
        end do 
    end do 

    do nu=1,nO
    do p=1,nO
        do q=1,nO
            do lamda=1,nO
                do sigma=1,nO
                    G2_MO(p,q,lamda,sigma)=G2_MO(p,q,lamda,sigma)+C(q,nu)*G2_Inter(p,nu,lamda,sigma)
                end do 
            end do 
        end do 
    end do
    end do 

    !G2_Inter=0.d0
    do p=1,nO
        do q=1,nO
            do r=1,nO
                do sigma=1,nO
                    do lamda=1,nO
                        G2_Inter(p,q,r,sigma)=G2_Inter(p,q,r,sigma)+C(r,lamda)*G2_MO(p,q,lamda,sigma)
                    end do 
                end do 
            end do 
        end do 
    end do

    do p=1,nO
        do q=1,nO
            do r=1,nO
                do s=1,nO
                    do sigma=1,nO
                        G2_MO(p,q,r,s)=G2_MO(p,q,r,s)+C(s,sigma)*G2_Inter(p,q,r,sigma)
                    end do 
                end do 
            end do 
        end do 
    end do 
end subroutine Convert_AOtoMO_2




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

subroutine Get_Fae   !#eq-3
    implicit none 
    integer::a,e,m,f,n
    Fae=0.d0
    do a=nElec+1,nSO
        do e=nElec+1,nSO
            do m=1,nElec
                do f=nElec+1,nSO
                    do n=1,nElec
                        Fae(a,e)=Fae(a,e)-0.5d0*td(m,n,a,f)*G2_SO(m,n,e,f)
                    end do 
                end do 
            end do 
        end do 
    end do
end subroutine Get_Fae

subroutine Get_Fmi   !#eq-4
    implicit none 
    integer::m,i,e,n,f
    Fmi=0.d0
    do m=1,nElec
        do i=1,nElec
            do e=nElec+1,nSO
                do n=1,nElec
                    do f=nElec+1,nSO
                        Fmi(m,i)=Fmi(m,i)+0.5d0*td(i,n,e,f)*G2_SO(m,n,e,f)
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
    do m=1,nElec
        do e=nElec+1,nSO
            Fme(m,e)=ESO(m,e)
        end do 
    end do 
end subroutine Get_Fme


subroutine Get_Wmnij !#eq-6
    implicit none 
    integer::m,n,i,j,e,f
    Wmnij=0.d0
    do m=1,nElec
        do n=1,nElec
            do i=1,nElec
                do j=1,nElec
                    Wmnij(m,n,i,j)=G2_SO(m,n,i,j)
                    do e=nElec+1,nSO
                        do f=nElec+1,nSO
                            Wmnij(m,n,i,j)=Wmnij(m,n,i,j)+0.25d0*td(i,j,e,f)*G2_SO(m,n,e,f)
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
    Wabef=0.d0
    do a=nElec+1,nSO
        do b=nElec+1,nSO
            do e=nElec+1,nSO
                do f=nElec+1,nSO
                    Wabef(a,b,e,f)=G2_SO(a,b,e,f)
                    do m=1,nElec
                        do n=1,nElec
                            Wabef(a,b,e,f)=Wabef(a,b,e,f)+0.25d0*td(m,n,a,b)*G2_SO(m,n,e,f)
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
    Wmbej=0.d0
    do m=1,nElec
        do b=nElec+1,nSO
            do e=nElec+1,nSO
                do j=1,nElec
                    Wmbej(m,b,e,j)=G2_SO(m,b,e,j)
                    do n=1,nElec
                        do f=nElec+1,nSO
                            !Wmbej(m,b,e,j)=Wmbej(m,b,e,j)
                            Wmbej(m,b,e,j)=Wmbej(m,b,e,j)-0.5d0*td(j,n,f,b)*G2_SO(m,n,e,f)
                        end do 
                    end do 
                end do 
            end do 
        end do 
    end do 
end subroutine Get_Wmbej

subroutine Update_DoubleExcitation
    implicit none 
    integer::i,j,a,b
    integer::e,m,n,f
    double precision:: Dijab

    call Get_Fae
    call Get_Fmi
    call Get_Fme

    call Get_Wmnij
    call Get_Wabef
    call Get_Wmbej

    tdnew=0.d0
    do i=1,nElec
        do j=1,nElec
            do a=nElec+1,nSO
                do b=nElec+1,nSO
                    Dijab=ESO(i,i)+ESO(j,j)-ESO(a,a)-ESO(b,b)

                    ! 1st term
                    tdnew(i,j,a,b)=G2_SO(i,j,a,b)

                    ! 2nd Term 
                    do e=nElec+1,nSO
                        tdnew(i,j,a,b)=tdnew(i,j,a,b)+td(i,j,a,e)*Fae(b,e)
                        tdnew(i,j,a,b)=tdnew(i,j,a,b)-td(i,j,b,e)*Fae(a,e)

                        do f=nElec+1,nSO
                            tdnew(i,j,a,b)=tdnew(i,j,a,b)+0.5d0*Wabef(a,b,e,f)*td(i,j,e,f)
                        end do 
                    end do 

                    ! 3rd Term 
                    do m=1,nElec 
                        tdnew(i,j,a,b)=tdnew(i,j,a,b)-td(i,m,a,b)*Fmi(m,j)
                        tdnew(i,j,a,b)=tdnew(i,j,a,b)+td(j,m,a,b)*Fmi(m,i)                    
                     end do

                    ! 4th Term 
                    do m=1,nElec
                        do e=nElec+1,nSO
                            tdnew(i,j,a,b)=tdnew(i,j,a,b)+td(i,m,a,e)*Wmbej(m,b,e,j)
                            tdnew(i,j,a,b)=tdnew(i,j,a,b)-td(j,m,a,e)*Wmbej(m,b,e,i)
                            tdnew(i,j,a,b)=tdnew(i,j,a,b)-td(i,m,b,e)*Wmbej(m,a,e,j)
                            tdnew(i,j,a,b)=tdnew(i,j,a,b)+td(j,m,b,e)*Wmbej(m,a,e,i)
                        end do
                    end do

                   tdnew(i,j,a,b)=tdnew(i,j,a,b)/Dijab
               end do 
            end do 
        end do 
    end do
    call Check_rmsT2Convergence
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


                    ! 7th term
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
