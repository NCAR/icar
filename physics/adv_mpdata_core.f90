! This is the core 1D implementation of mpdata
! It is stored in a separate file so that the same code can
! be used for advect_u, advect_v, and advect_w after l, r, etc. are set up

    ! First Calculate the standard upwind advection
    call flux1(l,r,U2,f)
    
    ! q1 is a temporary array to store the data while working with it here. 
    q1(1)=l(1)
    q1(n)=r(n-1)
    q1(2:n-1) = l(2:n-1) + (f(:n-2) - f(2:n-1))

    ! This is the MPDATA diffusion correction term for 1D flow
    ! U is defined on the mass grid for the pseudo-velocities?
    ! left and right pseudo (diffusive) velocities
    
    ! we will copy the q1 data into r to potentially minimize aliasing problems 
    ! for the compiler, and improve memory alignment for vectorization
    r  = q1(2:n) 
    ! l  = q1(1:n-1) ! no need to copy these data over again
    
    ! In MPDATA papers (r-l)/(r+l) is usually refered to as "A"
    ! compute the denomenator first so we can check that it is not zero
    denom=(r + q1(1:n-1))
    where(denom==0) denom=1e-10
    ! U2 is the diffusive pseudo-velocity
    U2 = abs(U2) - U2**2
    U2 = U2 * (r-q1(1:n-1)) / denom
    
    ! now calculate the MPDATA flux term
    call flux1(q1(1:n-1),r,U2,f)
    
! Fluxes are added to the original scalar field in the advect_u and advect_v subroutines
