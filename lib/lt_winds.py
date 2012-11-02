from __future__ import print_function
import sys
import numpy as np
from numpy import fft


def rotate_winds(rot_matrix,windfield):
    '''apply a matrix multiple to each element in an XxYxZx3 array
    
    WARNING, because of the nested for loops, this is likely to be slow
    but it is simple enough to convert to inline C or f2py fortran fairly easily.
    '''
    sz=windfield.shape
    for i in range(1,sz[0]-1):
        for j in range(1,sz[1]-1):
            for k in range(sz[2]):
                windfield[i,j,k,:]=np.dot(rot_matrix[i,j,k,...],windfield[i,j,k,:])

def generate_rotation_matrix(z,dx=4000.0):
    '''Calculate rotation matrices for each x,y,z model grid cell
    
    This is also likely to be painfully slow, and could be speed up
    However, this should only need to be run once at the begining, so
    it probably isnt worth spending much time on for now'''
    sz=z.shape
    rot_matrix=np.empty((sz[0],sz[1],sz[2],3,3),dtype='f')
    identity_matrix=np.array([[1,0,0],
                              [0,1,0],
                              [0,0,1]],dtype='f')
    print("Calculating wind Rotation Matrix... this will also take a while")
    ref_inc=sz[0]/10.0
    ref=ref_inc
    for i in range(1,sz[0]-1):
        if i>=ref:
            ref+=ref_inc
            print(str(np.round(float(i)/sz[0]*100))[0:4],end="% ")
            sys.stdout.flush()
        for j in range(1,sz[1]-1):
            for k in range(sz[2]):
                #normal vector 1 = cross product between x and y vectors on 3d grid
                nvec1=np.cross(np.array([dx,0,z[i-1,j,k]-z[i+1,j,k]],dtype='f'),
                             np.array([0,dx,z[i,j-1,k]-z[i,j+1,k]],dtype='f'))
                #normal vector 2 = (0,0,1) ?
                nvec2 = np.array([0,0,1],dtype='f')
                
                # order of vectors is important to get axis pointed in the correct direction
                # you can change the order, but then you have to change the sign of the angle theta
                rot_axis=np.cross(nvec2,nvec1) #the rotation axis is the cross product of the two normal vectors
                rot_axis/=np.sqrt(np.sum(rot_axis**2)) #normalize to a unit vector
                
                #angle = arccos( dot product / magnitudeA*magnitudeB )
                vdot=np.dot(nvec1,nvec2)
                mag1=np.sqrt(np.sum(nvec1**2))
                mag2=np.sqrt(np.sum(nvec2**2))
                # print(nvec1,nvec2,vdot,mag1,mag2)
                theta = np.arccos(vdot/(mag1*mag2))
                
                #3D rotation from the origin and about an arbitrary unit vector
                # from http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/
                # c1=1-cos(theta)
                # c=cos(theta)
                # s=sin(theta)
                # x1 = uu*c1*x + uv*c1*y + uw*c1*z + c*x + v*s*z-w*s*y
                # y1 = vu*c1*x + vv*c1*y + vw*c1*z + c*y + w*s*x-u*s*z
                # z1 = wu*c1*x + wv*c1*y + ww*c1*z + c*z + u*s*y-v*s*x
                # I've done more simplification, so be careful
                # x1 = x(uu*c1+c) + y(uv*c1-w*s) + z(uw*c1+v*s)
                # y1 = x(vu*c1+w*s) + y(vv*c1+c) + z(vw*c1-u*s)
                # z1 = x(wu*c1-v*s) + y(wv*c1+u*s) + z(ww*c1+c)
                c=np.cos(theta)
                c1=1-c
                s=np.sin(theta)
                u=rot_axis[0]
                v=rot_axis[1]
                w=rot_axis[2]
                rot_matrix[i,j,k,...]=np.array([[(u*u*c1+c),   (u*v*c1-w*s), (u*w*c1+v*s)],
                                     [(v*u*c1+w*s), (v*v*c1+c),   (v*w*c1-u*s)],
                                     [(w*u*c1-v*s), (w*v*c1+u*s), (w*w*c1+c)]])
                if not np.isfinite(rot_matrix[i,j,k,:,:].mean()):
                    rot_matrix[i,j,k,:,:]=identity_matrix
    print("100%")            
    return rot_matrix
    
def test_topo():
    # Make a theoretical gaussian mountain to test the code with. 
    # setup the model physical domain
    Lx    = 512.0*100     # % length of domain 
    Ly    = 512.0*100
    
    dx    = 100.0       # % grid spacing
    dy    = 100.0
    
    hm    = 1000.0         # % mnt height
    xm    = Lx/2        # % mountain location in domain 
    ym    = Ly/2
    am    = 5000.0       # % mountain half-width (m)
    
    Nx    = np.round(Lx/dx)+1   # % no of grid pts
    Ny    = np.round(Ly/dy)+1
    
    # % 2D distance arrays
    x     = np.linspace(0,Lx,Nx).reshape((1,Nx)).repeat(Ny,axis=0)  
    y     = np.linspace(0,Ly,Ny).reshape((Ny,1)).repeat(Nx,axis=1)
    
    # %----------------------------------------------------------------
    # % Make the mountain (theoretical)
    xyratio = 1. # ratio of x to y mountain half-width in domain
    zs=hm*np.exp(-(((x-xm)*xyratio)**2+(y-ym)**2)/am**2)  # % Gaussian
    
    Fzs   = fft.fftshift(fft.fft2(zs))/(Nx*Ny)
    zs=zs[np.newaxis,:,:].repeat(3,axis=0)
    zs[1,...]+=1000
    zs[2,...]+=2000
    return (zs,Fzs)

def linear_winds(Fzs, U,V,z,dx=4000.0,dy=4000.0,Ndsq  = 1E-5):
    # see Appendix A of Barstad and Gronas (2006) Tellus,58A,2-18
    # % -------------------------------------------------------------------------------
    # Ndsq  = 0.003**2      # % dry BV freq sq. was 0.01 initially Idar suggested 0.005
                            # should be calculated from the real atm profile
    f  = 0.938e-4           # rad/s Coriolis frequency for 40deg north
    # %---------------------------------------------------------------------------------
    
    (Ny,Nx)=Fzs.shape
    # % Compute 2D k and l wavenumber fields (this could be done once not on every pass)
    # (meshgrid may be faster?)
    k    = np.linspace(-np.pi/dx,np.pi/dx,Nx).reshape((1,Nx)).repeat(Ny,axis=0)
    l    = np.linspace(-np.pi/dy,np.pi/dy,Ny).reshape((Ny,1)).repeat(Nx,axis=1)
    i=1j #np.sqrt(np.complex(-1))
    
    sig  = U*k+V*l
    denom = sig**2-f**2
    kl = k**2+l**2
    kl[kl==0]=1E-30
    sig[sig==0]=1E-30
    
    msq = (Ndsq/denom * kl).astype('complex')          # % vertical wave number, hydrostatic
    mimag=np.zeros((Ny,Nx)).astype('complex')
    mimag.imag=(np.sqrt(-msq)).real
    m=np.where(msq>=0, (np.sign(sig)*np.sqrt(msq)).astype('complex'), mimag)

    m2 = np.sqrt(((Ndsq-sig**2)/denom * kl).astype('complex'))          # % vertical wave number, hydrostatic
    mimag.imag=(np.sqrt(-msq)).real
    m2=np.where(msq>=0, (np.sign(sig)*np.sqrt(msq)).astype('complex'), mimag)
    
    neta=Fzs*np.exp(i*m2*z)
    w_hat=i*sig*neta
    u_hat=(-m2*(sig*k-i*l*f)*i*neta)/kl
    v_hat=(-m2*(sig*l+i*k*f)*i*neta)/kl
    
    w_hat=Ny*Nx*np.real(fft.ifft2(fft.ifftshift(w_hat)))
    u_hat=Ny*Nx*np.real(fft.ifft2(fft.ifftshift(u_hat)))
    v_hat=Ny*Nx*np.real(fft.ifft2(fft.ifftshift(v_hat)))
    return(u_hat,v_hat,w_hat)
    
def update_winds(z,Fzs,U,V,W,dx=4000.0,Ndsq=1E-5,r_matrix=None,padx=0,pady=0,rotation=True):
    
    if r_matrix==None and rotation:
        r_matrix=np.array(generate_rotation_matrix(z.transpose((1,2,0)),dx=dx))
    
    sz=U.shape
    windfield=np.empty((sz[1],sz[2],sz[0],3),dtype="f")
    endx=np.choose(padx==0,[-padx,None])
    endy=np.choose(pady==0,[-pady,None])
    for i in range(U.shape[0]):
        uh,vh,wh=linear_winds(Fzs,U[i,:,:].mean(),V[i,:,:].mean(),(z[i,:,:]-z[0,:,:]).mean(),dx=dx,dy=dx,Ndsq=Ndsq)
        windfield[:,:,i,0]=U[i,:,:]+uh[pady:endy,padx:endx]
        windfield[:,:,i,1]=V[i,:,:]+vh[pady:endy,padx:endx]
        windfield[:,:,i,2]=W[i,:,:]+wh[pady:endy,padx:endx]
    if rotation:
        rotate_winds(r_matrix,windfield)
    for i in range(U.shape[0]):
        U[i,:,:]=windfield[:,:,i,0]
        V[i,:,:]=windfield[:,:,i,1]
        W[i,:,:]=windfield[:,:,i,2]
    return r_matrix