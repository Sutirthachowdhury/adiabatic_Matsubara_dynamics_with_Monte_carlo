subroutine parameter_input
    use global 
    
    implicit none

    integer i
    
    !=========================================
    open(10,file='params.in',status='old')
    read(10,*) nstep
    read(10,*) ntraj
    read(10,*) neq
    read(10,*) nmc
    read(10,*) beta
    read(10,*) mnuc
    read(10,*) dt
    read(10,*) ourseed
    close(10)
    !==========================================

    allocate(step(nstate+1))
    allocate(a_4(lb_m:ub_m,lb_m:ub_m,lb_m:ub_m,lb_m:ub_m))

    open(10,file='steps.in',status='old')
    read(10,*) step(1), step(2), step(3)
    close(10)

    !--- freq of HO-----
    omega(1) = 1.0    
    !--------------------
   
    ! kronicker delta for mapping hamiltonian (we will use later)---- 

     del=0.

     do i=1,nstate
        del(i,i)=1.
     end do

    !----------------------------------- 

end subroutine parameter_input
