Module MBPT
    use HartreeFock
	double precision::G2_MO(K,K,K,K)
    double precision::EMP2
    contains
    subroutine Convert_AOtoMO
        implicit none 
        integer::i,j,a,b
        integer::mu,nu,lamda,sigma
        double precision::muINT,nuINT,lamdaINT,sigmaINT
        G2_MO=0.0
        do i=1,K
            do j=1,K
                do a=1,K
                    do b=1,K
                        sigmaINT=0.d0
                        do  sigma=1,K
                            lamdaINT=0.d0
                            do lamda=1,K
                                nuINT=0.d0
                                do nu=1,K
                                    muINT=0.d0
                                    do mu=1,K
                                        muINT=muINT+C(mu,i)*G2(mu,nu,lamda,sigma)
                                    end do
                                    nuINT=nuINT+C(nu,a)*muINT
                                end do 
                                lamdaINT=lamdaINT+C(lamda,j)*nuINT
                            end do 
                            sigmaINT=sigmaINT+C(sigma,b)*lamdaINT
                        end do
                        G2_MO(i,a,j,b)=G2_MO(i,a,j,b)+sigmaINT 
                    end do 
                end do 
            end do 
        end do 
         
    end subroutine Convert_AOtoMO


    subroutine Convert_AOtoMO
        implicit none 
        integer::i,j,a,b
        integer::mu,nu,lamda,sigma
        double precision::muINT,nuINT,lamdaINT,sigmaINT
        G2_MO=0.0
        

         
    end subroutine Convert_AOtoMO
    

    subroutine Get_EMP2
		implicit none 
	    integer::i,j,a,b
        double precision::numerator,denominator
        EMP2=0.0
        do i=1,2
            do j=1,2
                do a=3,6
                    do b=3,6
                        numerator=G2_MO(i,a,j,b)*(2*G2_MO(i,a,j,b)-G2_MO(i,b,j,a))
                        denominator=MOEnergy(i)+MOEnergy(j)-MOEnergy(a)-MOEnergy(b)
                        EMP2=EMP2+numerator/denominator
                    end do 
                end do 
            end do 
        end do 
		print*,EMP2 
    end subroutine Get_EMP2
End Module MBPT
