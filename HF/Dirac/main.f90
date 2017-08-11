program main
    use HartreeFock
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
        call Get_newEnergy()
        call Get_Convergence()
        call PrintResults()
        if(Converged)exit SCF
        D=newD
        E=newE
    end do SCF

end program main
