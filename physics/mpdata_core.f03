    ! This is the core 1D implementation of mpdata
    ! It is stored in a separate file so that the same code can
    ! be used for advect_u, advect_v, and advect_w after l, r, etc. are set up
    call flux1(l,r,U2,f)

    q1(1)=l(1)
    q1(n)=r(n-1)
    q1(2:n-1) = l(2:n-1) + (f(:n-2) - f(2:n-1))

    ! This is the MPDATA diffusion correction term for 1D flow
    ! l, c, and r are now len(n-2) defined with q from the first update
    l2  = q1(1:n-2)
    c2  = q1(2:n-1)
    r2  = q1(3:n)
    ! U is defined on the mass grid for the pseudo-velocities?
    U3 = (U2(1:n-2) + U2(2:n-1)) / 2
    U3 = abs(U3)-U3**2
    ! left and right pseudo (diffusive) velocities
    denom=(q1(1:n-2)+q1(2:n-1))
    where(denom==0) denom=1e-10
    Vl = U3 * (q1(1:n-2)-q1(2:n-1)) / denom

    denom=(q1(2:n-1)+q1(3:n))
    where(denom==0) denom=1e-10
    Vr = U3 * (q1(3:n)  -q1(2:n-1)) / denom

    call flux1(l2,c2,Vl,fl)
    call flux1(c2,r2,Vr,fr)
