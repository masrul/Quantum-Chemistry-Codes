program main
    use HartreeFock
    use MBPT
    call Get_HFEnergy
    call Convert_AOtoMO
    call Get_EMP2
    print*,MOEnergy(:)
    print*,'G2:',G2_MO(3,4,5,6)
    print*,'C:',C(1,:)
end program main
