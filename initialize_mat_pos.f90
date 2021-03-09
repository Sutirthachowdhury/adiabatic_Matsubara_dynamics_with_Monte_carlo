subroutine initialize_mat_pos(Q_mats)

    use global

    real*8, intent(out) :: Q_mats(ndof,lb_m:ub_m)

    integer :: i,j,k

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Initialize positions 
    do j = 1,ndof
        do i=lb_m,ub_m
            Q_mats(j,i) = 0.
        enddo
    enddo   

end subroutine initialize_mat_pos