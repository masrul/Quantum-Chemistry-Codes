program spin
    integer::i,j,k,l
    integer::alpha1,beta1    
    integer::alpha2,beta2    
    integer::alpha3,beta3    
    integer::alpha4,beta4

    do i=1,6
        alpha1=2*i-1
        beta1=alpha1+1
        do j=1,6
            alpha2=2*j-1
            beta2=alpha2+1
            do k=1,6
                alpha3=2*k-1
                beta3=alpha3+1
                do l=1,6
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
end program spin
