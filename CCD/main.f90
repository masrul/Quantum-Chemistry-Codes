program main
    use HF
    use CCD
    implicit none 
    integer::iter
    call Get_EHF
    call Convert_AOtoMO
    call Convert_EMO_2_ESO
    call Convert_Spatial2Spin
    call Init_CCD

    iter=0
    do while (.not.ccdConverged)
        iter=iter+1
        call Get_ECCD
        call Update_DoubleExcitation
        print*,'Iteration number: ',iter,'-------> ECCD: ',ECCD
    end do 
end program main
