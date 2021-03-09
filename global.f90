module global
    implicit none
  
    integer, parameter:: nosc=1
    integer, parameter:: nstate=2
    integer,parameter:: ndof=nosc
   
  
    real*8 omega(ndof) ! harmonic oscillator frequency

    integer, parameter :: nb = 115  ! 5*M
    integer,parameter :: mats_mode = 3
    integer, parameter :: ub_m   = +(mats_mode-1)/2
    integer, parameter :: lb_m   = -(mats_mode-1)/2
    
    real*8 mnuc,beta,betap,dt,delta  !delta is the diabatic coupling
    
    integer nstep,ntraj,nskip,nsample,neq,nmc
  
    integer :: ourseed
  
    real*8  del(nstate,nstate) ! delta symbol

    !* monte carlo step size
    real*8,allocatable :: step(:)

    ! this is the normal mode coefficients, this needs to be called once
    real*8, allocatable :: a_4(:,:,:,:)
  
    complex*16,parameter :: eye=(0.,1.)
  
    real,parameter :: pi=acos(-1.)
    real*8, parameter:: hbar=1.
  end module global
  