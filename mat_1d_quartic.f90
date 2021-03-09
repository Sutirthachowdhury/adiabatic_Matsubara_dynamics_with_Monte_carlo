program mat_1d_quartic

    use global
    implicit none

    integer :: istep,iup,idn,isample,ieq     
    integer,allocatable :: seed(:)
    integer seed_dimension,time


    integer i,j,k,l,m,n,ncount

    !======= premitive bead variables=============

    !real*8, allocatable :: xe_mc(:,:),xe(:,:) ! nuclear bead positions 
    !real*8, allocatable :: ve(:,:) ! nuclear beads momentum
    !=======================================


    real*8, allocatable:: Q_mats(:,:),P_mats(:,:),f_mats(:,:) !"Matsubara" positions, momentum and force
    real*8, allocatable:: Q_mats_dyn(:,:),P_mats_dyn(:,:),f_mats_dyn(:,:) ! this variables will use in dynamics
    real*8, allocatable :: matsfreq(:)
    real*8, allocatable :: corr_xx(:) ! this is RR auto-correlation

    real*8 :: dp,gran,rand,theta 
    real*8 :: corr_xx0 ! initial value of correlation function

    real*8 :: partition ! weighting the cosine term

    real*8 :: q(ndof,nb)

    ! getting the inputs
    call parameter_input
 
    !-------- random seed-------------------
    CALL RANDOM_SEED(size=seed_dimension)
    ALLOCATE (seed(seed_dimension))
    do i=1,seed_dimension
      seed(i) =time()+3*i-1 + ourseed
    end do
    CALL RANDOM_SEED(PUT=seed)
    !-----------------------------------------

    !allocate(xe_mc(ndof,nb),xe(ndof,nb),ve(ndof,nb))
    allocate(Q_mats(ndof,lb_m:ub_m),P_mats(ndof,lb_m:ub_m),f_mats(ndof,lb_m:ub_m))
    allocate(Q_mats_dyn(ndof,lb_m:ub_m),P_mats_dyn(ndof,lb_m:ub_m),f_mats_dyn(ndof,lb_m:ub_m))
    allocate(matsfreq(lb_m:ub_m))
    
    allocate(corr_xx(nstep))

    !--- initialization of correlation function 
    corr_xx = 0.0
    corr_xx0 = 0.0
    !----------------------------

    !------- initializing the partition function----

    partition = 0.0

    !---- initialize the matsubara pos,momentum and potential----
    call initialize_mat_pos(Q_mats)
    call momentum_sample(P_mats) ! for initializing momentum from Maxwell-Boltzmann distribution
    call init_pot() ! initialize potential variables (coeff for infinity N limit tranformation)
    !===========================================================

       
    !----- run thermalization trajectory---
    do j = 1,neq
      !if (mod(ieq,nskip)==0) then
      !  call momentum_sample(P_mats)
      !end if  
      call monte_carlo(Q_mats)
      !call evolve(Q_mats,P_mats)

    enddo 
    !----------------------------------------


    do j = 1,ntraj
     
      !do isample = 1,nsample
  
      !   call evolve(Q_mats,P_mats)
         !============================
 
      ! enddo

      call monte_carlo(Q_mats)
       
       !store sample configuration
       Q_mats_dyn = Q_mats
  
       call momentum_sample(P_mats)

      !------------- calculation of phase term----------
      do i = lb_m,ub_m
        matsfreq(i) = (2.0*pi*i)/beta
      enddo  

      theta = 0.
      do i=lb_m,ub_m
        theta = theta + matsfreq(i)*Q_mats_dyn(1,-i)*P_mats(1,i)
      enddo

      !---------- calculating partiton function (accumulation of real(phase))-----

      partition = partition + cos(beta*theta)


      !----------- C_xx(t=0)----------------------

      corr_xx0 = corr_xx0 + cos(beta*theta)*Q_mats_dyn(1,0)*Q_mats_dyn(1,0)
  
      !================dynamics step ==========================================
     
      do istep = 1,nstep
      

        call evolve(Q_mats_dyn,P_mats)

        call back_transform_matrix(Q_mats_dyn,q)

        write(201,222) real(istep*dt), (q(1,k),k=1,nb)


        corr_xx(istep) = corr_xx(istep) + cos(beta*theta)*Q_mats(1,0)*Q_mats_dyn(1,0)

      enddo 

      write(201,222)
      write(201,222)

    enddo 

    corr_xx = corr_xx/real(partition)
    corr_xx0 = corr_xx0/real(partition)

    !write(3001,*) corr_xx0

    !do ncount=1,nstep
    !  write(200,222) ncount*dt,corr_xx(ncount)
    !end do

222      format(300(e13.6,2x))
    stop
end program mat_1d_quartic