subroutine monte_carlo(Q_mats)

    use global
    implicit none

    integer :: tot_mode,imode
    integer :: init_nmc
    integer :: i,j,iup,idn,k,init_mc,l


    real*8 :: iacc_x,iatt_x,rnum,rand,pot
    real*8 :: energy,energy_old,del_e,wt
    real*8 :: q(ndof,nb)

    real*8, allocatable :: Q_new(:,:)

    real*8, intent(inout) :: Q_mats(ndof,lb_m:ub_m)

    allocate(Q_new(ndof,lb_m:ub_m))

  !initialize all zeroes

  iatt_x = 0.d0
  iacc_x = 0.d0

  do init_mc=1,nmc

    !-----randomly pick a mode to move-----                                 
                                                
    call random_number(rand)
    tot_mode = int(rand*mats_mode)

    imode = int(tot_mode-1)
   

    !--------------------------------
    
    !move individual modes
    !generating trial moves
        
    do j=1,ndof

      Q_new(j,:) = Q_mats(j,:)    
      call random_number(rand)
      rnum = rand
      Q_new(j,imode) = Q_mats(j,imode) + (rnum-0.5d0)*step(1)  ! here need to check with gran as well
      
   end do

   !------------- compute new energy----------------
   call back_transform_matrix(Q_new,q)

     !--- compute the potential in bead representation------

   pot = 0.

   do l = 1,nb
     pot = pot + 0.25*q(1,l)**4
   enddo 

   pot = pot/real(nb)


  energy = pot
  !------------------------------------------------------
  iatt_x=iatt_x+1.0

 
  !-----------computing old energy------------------
  
  call back_transform_matrix(Q_mats,q)

  !--- compute the potential in bead representation------

  pot = 0.

  do l = 1,nb
    pot = pot + 0.25*q(1,l)**4
  enddo 

  pot = pot/real(nb)

  energy_old = pot
  !------------------------------------------

  ! calculating the energy difference
  del_e = (energy-energy_old)

  !accepting/rejecting step
  wt = dexp(-beta*del_e)
     
     
  call random_number(rand)
  rnum=rand

  if(rnum.lt.wt) then
           
    Q_mats(:,imode) = Q_new(:,imode)
    iacc_x=iacc_x+1.d0
       
  endif

  enddo 

end subroutine monte_carlo  