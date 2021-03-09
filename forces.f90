subroutine forces(Q_mats,f_mats)

    use global

    implicit none


    real*8, intent(in) :: Q_mats(ndof,lb_m:ub_m)    
    real*8, intent(out) :: f_mats(ndof,lb_m:ub_m)
    integer :: i,j,k,l
    real*8 :: pot
  

    real*8 :: q(ndof,nb),f_local(ndof,nb)
    
    call back_transform_matrix(Q_mats,q)

    !--- compute the potential in bead representation------

    pot = 0.

    do l = 1,nb
      pot = pot + 0.25*q(1,l)**4
    enddo 

    pot = pot/real(nb)
    !----------------------------------------
    !--------- force in bead rrepresentation-------
    f_local = 0.0d0

    do l = 1,nb
      f_local(1,l) = -q(1,l)**3
    enddo 
    !-------------------------------------

    call transform_matrix(f_local,f_mats)


    !!! Compute the quartic force term

    !f_mats(:,:) = 0.0d0

    !do i = lb_m,ub_m
    !  do j = lb_m,ub_m
    !    do k = lb_m,ub_m
    !      do l = lb_m,ub_m
    !        f_mats(:,i) = f_mats(:,i) - Q_mats(:,j)*Q_mats(:,k)*Q_mats(:,l)*a_4(i,j,k,l)*4.
    !      enddo 
    !    enddo 
    !  enddo 
    !enddo 

    
  
  return

end subroutine forces   