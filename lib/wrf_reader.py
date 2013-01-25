import glob
import numpy as np

from bunch import Bunch
import swim_io
import fast_mean # inline C to do a running spatial mean

R=8.3144621 # J/mol/K
cp=29.19 # J/mol/K   =1.012 J/g/K
g=9.81 # m/s^2

def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)    
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)               
    x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
    g = np.exp(-(x**2/float(size)+y**2/float(sizey)))
    g=g[size/2:-size/2,sizey/2:-sizey/2]
    return g / g.sum()


def match_xy(lat1,lon1,lat2,lon2):
    N=lat2.shape
    x=np.zeros(N)
    y=np.zeros(N)
    winhalfsize=2

    for i in range(N[1]):
        j=0
        dists=(lat1-lat2[j,i])**2 + (lon1-lon2[j,i])**2
        (lasty,lastx)=np.unravel_index(dists.argmin(), dists.shape)
        (prevx,prevy)=(lastx,lasty)
        latwin=lat1[lasty-winhalfsize:lasty+winhalfsize,lastx-winhalfsize:lastx+winhalfsize]
        lonwin=lon1[lasty-winhalfsize:lasty+winhalfsize,lastx-winhalfsize:lastx+winhalfsize]
        for j in range(N[0]):
            if (prevx!=lastx) or (prevy!=lasty):
                latwin=lat1[lasty-winhalfsize:lasty+winhalfsize,lastx-winhalfsize:lastx+winhalfsize]
                lonwin=lon1[lasty-winhalfsize:lasty+winhalfsize,lastx-winhalfsize:lastx+winhalfsize]
                (prevx,prevy)=(lastx,lasty)
            dists=(latwin-lat2[j,i])**2+(lonwin-lon2[j,i])**2
            (newy,newx)=np.unravel_index(dists.argmin(), dists.shape)
            
            # print(j,i,newx,newy)
            lastx=(newx-winhalfsize)+lastx
            lasty=(newy-winhalfsize)+lasty
            x[j,i]=lastx
            y[j,i]=lasty
        
    return (x,y)
 
def bilin_weights(yi,y,xi,x):
    x0=np.abs((xi-x[0])/(x[1]-x[0]))
    x1=1-x0 # np.abs((xi-x[1])/(x[1]-x[0]))
    x2=np.abs((xi-x[2])/(x[3]-x[2]))
    x3=1-x2 # np.abs((xi-x[3])/(x[3]-x[2]))
    y5=y[0]*x1+y[1]*x0
    y6=y[2]*x3+y[3]*x2
    f1=(yi-y5)/(y6-y5)
    f2=1-f1# (y6-yi)/(y6-y5)
    return np.array([x1*f2,x0*f2,x3*f1,x2*f1])
    
    
def match_xy_bilin(lat1,lon1,lat2,lon2):
    N1=lat1.shape
    N=lat2.shape
    out=np.zeros((N[0],N[1],4,3))
    x=np.zeros(4).astype('i')
    y=np.zeros(4).astype('i')
    winhalfsize=3
    
    dxinc=np.sign(lon1[1,1]-lon1[0,0]).astype('i')
    dyinc=np.sign(lat1[1,1]-lat1[0,0]).astype('i')

    for i in range(N[1]):
        j=0
        dists=(lat1-lat2[j,i])**2 + (lon1-lon2[j,i])**2
        (lasty,lastx)=np.unravel_index(dists.argmin(), dists.shape)
        (prevx,prevy)=(lastx,lasty)
        ymin=lasty-winhalfsize
        if ymin<0:ymin=0
        ymax=lasty+winhalfsize
        if ymax>=N1[0]:ymax=N1[0]-1
        xmin=lastx-winhalfsize
        if xmin<0:xmin=0
        xmax=lastx+winhalfsize
        if xmax>=N1[1]:xmax=N1[1]-1
        latwin=lat1[ymin:ymax,xmin:xmax]
        lonwin=lon1[ymin:ymax,xmin:xmax]
        for j in range(N[0]):
            if (prevx!=lastx) or (prevy!=lasty):
                ymin=lasty-winhalfsize
                if ymin<0:ymin=0
                ymax=lasty+winhalfsize
                if ymax>=N1[0]:ymax=N1[0]-1
                xmin=lastx-winhalfsize
                if xmin<0:xmin=0
                xmax=lastx+winhalfsize
                if xmax>=N1[1]:xmax=N1[1]-1
                latwin=lat1[ymin:ymax,xmin:xmax]
                lonwin=lon1[ymin:ymax,xmin:xmax]
                (prevx,prevy)=(lastx,lasty)
            dists=(latwin-lat2[j,i])**2+(lonwin-lon2[j,i])**2
            (newy,newx)=np.unravel_index(dists.argmin(), dists.shape)
            lastx=xmin+newx
            lasty=ymin+newy
            x[0]=newx
            y[0]=newy
            if latwin[newy,newx]<lat2[j,i]:
                if lonwin[newy,newx]<lon2[j,i]:
                    x[1]=newx+dxinc
                    x[2]=newx
                    x[3]=newx+dxinc
                    y[1]=newy
                    y[2]=newy+dyinc
                    y[3]=newy+dyinc
                else:
                    x[1]=newx-dxinc
                    x[2]=newx
                    x[3]=newx-dxinc
                    y[1]=newy
                    y[2]=newy+dyinc
                    y[3]=newy+dyinc
            else:
                if lonwin[newy,newx]<lon2[j,i]:
                    x[1]=newx+dxinc
                    x[2]=newx
                    x[3]=newx+dxinc
                    y[1]=newy
                    y[2]=newy-dyinc
                    y[3]=newy-dyinc
                else:
                    x[1]=newx-dxinc
                    x[2]=newx
                    x[3]=newx-dxinc
                    y[1]=newy
                    y[2]=newy-dyinc
                    y[3]=newy-dyinc
            x[x<0]=0
            y[y<0]=0
            x[x>=(xmax-xmin)]=xmax-xmin-1
            y[y>=(ymax-ymin)]=ymax-ymin-1
                    
            # bilinear interpolation for an arbitrary grid spacing 
            # (must be a grid, but this will handle a lot of irregularities)
            weights=bilin_weights(lat2[j,i],latwin[y,x],lon2[j,i],lonwin[y,x])

            out[j,i,:,0]=y+ymin
            out[j,i,:,1]=x+xmin
            out[j,i,:,2]=weights
        
    return out
              
class WRF_Reader(object):
    _filenames=None
    _sfc_files=None
    _curfile=0
    curpos=0
    x=-99
    y=-99
    _bilin=False
    _nn=True
    _geoLUT=None
    topo=None
    # atm variables in WRF 3d files
    timevar="Times"      # time...                                  1D
    hvar='HGT'           #terran height(?)                  [m]     2D
    tvar='T'             #potential temperature (-300)      [K]     3D
    gphvar_base="PHB"    #base geopotential height          [m]     3D
    gphvar_pert="PH"     #geopotential height perturbation  [m]     3D
    pvar_base="PB"       #base pressure                     [pa]    3D
    pvar_pert="P"        # pressure perturbation            [pa]    3D
    qcvar='QCLOUD'       # cloud water mixing ratio         [kg/kg] 3D
    qvvar='QVAPOR'       # water vapor mixing ratio         [kg/kg] 3D
    qivar='QICE'         # cloud ice mixing ratio           [kg/kg] 3D
    wvar='W'             # vertical winds                   [pa/s]  3D NOT USED
    uvar='U'             # U wind component                 [m/s]   3D 
    vvar='V'             # V wind component                 [m/s]   3D 
    # land surface variables in WRF files
    pblvar='PBLH'        # PBL height                       [m]     2D
    lhvar="LH"           # Latent heat flux from surface    [W/m^2] 2D
    shvar="HFX"          # Sensible heat flux from surface  [W/m^2] 2D

    def make_model_domain(self,nlevels,maxheight,topo):
        sz=topo.shape
        
        slopes=(maxheight-topo)/nlevels
        if self.use_linear_winds:
            slopes=slopes.mean()
            self.hgt3d=(topo[np.newaxis,...].repeat(nlevels,axis=0)+
                    slopes[np.newaxis,...]*np.arange(nlevels)[:,np.newaxis,np.newaxis])
        elif self.use_wrf_winds:
            self.hgt3d=(swim_io.read_nc(self._filenames[0],var=self.gphvar_base)
                    .data[0,self.usepressures,self.sub_y0:self.sub_y1,self.sub_x0:self.sub_x1]/9.81)
            if topo.shape[0]<self.hgt3d.shape[1]:
                halfk=(self.hgt3d.shape[1]-topo.shape[0])/2
                self.hgt3d=self.hgt3d[:,halfk:-halfk,halfk:-halfk]
                
        else:
            self.hgt3d=(topo[np.newaxis,...].repeat(nlevels,axis=0)+
                    slopes[np.newaxis,...].repeat(nlevels,axis=0)*
                    np.arange(nlevels)[:,np.newaxis].repeat(sz[0]*sz[1],axis=1).reshape((nlevels,sz[0],sz[1])))
        
    
    def init_xy(self,driverfilename='forcing/wrfout_d01_2000-10-01_00:00:00',domain=None,options=None):
        if options.wrfres==2000:
            wrffilename=options.baseline_file
        elif options.wrfres==4000:
            wrffilename=options.baseline_file
        halfk=options.halfk
        sub_y0,sub_y1,sub_x0,sub_x1=options.subset
        self.sub_y0=sub_y0
        self.sub_y1=sub_y1
        self.sub_x0=sub_x0
        self.sub_x1=sub_x1
        d=swim_io.Dataset(driverfilename, 'r')
        nlat=d.variables['XLAT'][0,:,:]
        nlon=d.variables['XLONG'][0,:,:]
        nulat=d.variables['XLAT_U'][0,:,:]
        nulon=d.variables['XLONG_U'][0,:,:]
        nvlat=d.variables['XLAT_V'][0,:,:]
        nvlon=d.variables['XLONG_V'][0,:,:]
        hgt=d.variables['HGT'][0,:,:]
        self.nlevels=options.nlevels
        self.usepressures=np.arange(self.nlevels)
        self.base_pressure=d.variables[self.pvar_base][0,...][self.usepressures,...]
        self.base_gph=d.variables[self.gphvar_base][0,...][self.usepressures,...]
        maxheight=self.base_gph[-1,:,:]/g
        d.close()
        wlat=swim_io.read_nc(wrffilename,var='XLAT').data[0,sub_y0:sub_y1,sub_x0:sub_x1]
        wlon=swim_io.read_nc(wrffilename,var='XLONG').data[0,sub_y0:sub_y1,sub_x0:sub_x1]
        if domain!=None:
            hires_topo=domain.topo
        else:
            hires_topo=swim_io.read_nc(wrffilename,var='HGT').data[0,sub_y0:sub_y1,sub_x0:sub_x1][halfk:-halfk,halfk:-halfk]
        print('Calculating XY match lookup table, this may take a while.')
        if self._nn:
            wlat=wlat[halfk:-halfk,halfk:-halfk]
            wlon=wlon[halfk:-halfk,halfk:-halfk]
            (x,y)=match_xy(nlat,nlon,wlat,wlon)
            (xu,yu)=match_xy(nulat,nulon,wlat,wlon)
            (xv,yv)=match_xy(nvlat,nvlon,wlat,wlon)
            self.x=x.astype('i')
            self.y=y.astype('i')
            self.xu=xu.astype('i')
            self.yu=yu.astype('i')
            self.xv=xv.astype('i')
            self.yv=yv.astype('i')
            self.topo=hgt[self.y,self.x]
            maxheight=maxheight[self.y,self.x]
            self.base_gph=self.base_gph[:,self.y,self.x]
            self.base_pressure=self.base_pressure[:,self.y,self.x]
            self.make_model_domain(self.nlevels,maxheight,hires_topo)
            
        elif self._bilin:
            self._geoLUT=match_xy_bilin(nlat,nlon,wlat,wlon)
            self._geoLUTu=match_xy_bilin(nulat,nulon,wlat,wlon)
            self._geoLUTv=match_xy_bilin(nvlat,nvlon,wlat,wlon)
            curx=(self._geoLUT[:,:,:,1]).astype('i')
            cury=(self._geoLUT[:,:,:,0]).astype('i')
            w=self._geoLUT[:,:,:,2]
            self.topo=np.zeros(curx.shape[0:2],dtype=np.float32)
            for i in range(4):self.topo+=hgt[cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            self.topo=self.topo[halfk:-halfk,halfk:-halfk]
            MH=np.zeros(curx.shape[0:2],dtype=np.float32)
            for i in range(4):MH+=maxheight[cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            MH=MH[halfk:-halfk,halfk:-halfk]
            self.make_model_domain(self.nlevels,MH,hires_topo)
            
        
    def __init__(self, file_search,sfc=None,nn=True, bilin=False, windonly=False, domain=None, options=None,*args, **kwargs):
        super(WRF_Reader,self).__init__(*args, **kwargs)
        self.use_linear_winds=options.use_linear_winds
        self.use_wrf_winds=options.use_wrf_winds
        self.halfk=options.halfk
        self._filenames=np.sort(glob.glob(file_search))
        self.io_ratio=options.io_ratio
        self.windonly=windonly
        if sfc!=None:
            self._sfc_files=glob.glob(sfc)
        if nn:
            self._nn=True
            self._bilin=False
        elif bilin:
            self._nn=False
            self._bilin=True
        self.init_xy(self._filenames[0],domain=domain,options=options)
        self.curpos=options.start_position
        if windonly:
            self.curpos=0
    
    
    # we are our own iterator...
    def __iter__(self):
        return self

    def nextNN(self):
        curfile=self._curfile
        if curfile>=len(self._filenames):
            raise StopIteration
        minx=self.x.min()
        maxx=self.x.max()+1
        miny=self.y.min()
        maxy=self.y.max()+1
        curx=self.x-minx
        cury=self.y-miny

        minxu=self.xu.min()
        maxxu=self.xu.max()+1
        minyu=self.yu.min()
        maxyu=self.yu.max()+1
        curxu=self.xu-minxu
        curyu=self.yu-minyu

        minxv=self.x.min()
        maxxv=self.x.max()+1
        minyv=self.y.min()
        maxyv=self.y.max()+1
        curxv=self.x-minxv
        curyv=self.y-minyv
        plevels=self.usepressures
        d=swim_io.Dataset(self._filenames[curfile], 'r')
        if self.windonly:
            wind_u=d.variables[self.uvar][self.curpos,:,minyu:maxyu,minxu:maxxu][plevels,...]
            wind_u=wind_u[:,curyu,curxu]
            wind_v=d.variables[self.vvar][self.curpos,:,minyv:maxyv,minxv:maxxv][plevels,...]
            wind_v=wind_v[:,curyv,curxv]
            hgt=d.variables[self.gphvar_pert][self.curpos,:,miny:maxy,minx:maxx][plevels,...]
            hgt=(hgt[:,cury,curx]+self.base_gph)/9.8
            temperature=d.variables[self.tvar][self.curpos,:,miny:maxy,minx:maxx][plevels,...]
            temperature=temperature[:,cury,curx]+300
            specific_humidity=d.variables[self.qvvar][self.curpos,:,miny:maxy,minx:maxx][plevels,...]
            specific_humidity=specific_humidity[:,cury,curx]
            pressure=d.variables[self.pvar_pert][self.curpos,:,miny:maxy,minx:maxx][plevels,...]
            pressure=pressure[:,cury,curx]+self.base_pressure
            # temperature=None
            # pressure=None
            # specific_humidity=None
            relative_humidity=None
        else:
            temperature=d.variables[self.tvar][self.curpos,:,miny:maxy,minx:maxx][plevels,cury,curx]
            pressure=d.variables[self.pvar_pert][self.curpos,:,miny:maxy,minx:maxx][plevels,cury,curx]
            specific_humidity=d.variables[self.qvvar][self.curpos,:,miny:maxy,minx:maxx][plevels,cury,curx]
            relative_humidity=specific_humidity.copy()
            wind_u=d.variables[self.uvar][self.curpos,:,minyu:maxyu,minxu:maxxu][plevels,curyu,curxu]
            wind_v=d.variables[self.vvar][self.curpos,:,minyv:maxyv,minxv:maxxv][plevels,curyv,curxv]
            hgt=(d.variables[self.gphvar_pert][self.curpos,:,miny:maxy,minx:maxx][plevels,cury,curx]+self.base_gph)/9.8

        datestr=''.join(d.variables[self.timevar][self.curpos])
        d.close()


        self.time_inc()
        N=list(wind_v.shape)
        N[2]+=1 #v is staggered in y direction
        return Bunch(th=temperature, p=pressure,hgt=hgt,
                     qv=specific_humidity/(1-specific_humidity), 
                     rh=relative_humidity, 
                     u=wind_u, v=wind_v,w=np.zeros(N,dtype=np.float32,order="F"),
                     date=datestr)
    

    def _next_surface(self,curpos):
        geoLUT=self._geoLUT
        x=geoLUT[:,:,:,1]
        y=geoLUT[:,:,:,0]
        w=geoLUT[:,:,:,2]
        halfk=self.halfk
        offset=np.int(np.ceil(halfk/self.io_ratio))
        minx=x.min().astype('i')-offset 
        maxx=x.max().astype('i')+offset
        miny=y.min().astype('i')-offset
        maxy=y.max().astype('i')+offset
        curx=(x-minx).astype('i')
        cury=(y-miny).astype('i')
        d=swim_io.Dataset(self._sfc_files[self._curfile], 'r')
        Nxy=curx.shape[0:2]
        sensible_heat=np.zeros(Nxy,dtype=np.float32)
        latent_heat=np.zeros(Nxy,dtype=np.float32)
        pblh=np.zeros(Nxy,dtype=np.float32)
        # longwave_up=np.zeros(Nxy,dtype=np.float32)
        # albedo=np.zeros(Nxy,dtype=np.float32)
        
        gaussian=gauss_kern(2)
        
        curdata=d.variables[self.shvar][curpos,miny:maxy,minx:maxx]
        for i in range(4):sensible_heat+=np.float32(curdata[cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        curdata=d.variables[self.lhvar][curpos,miny:maxy,minx:maxx]
        for i in range(4):latent_heat+=np.float32(curdata[cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        curdata=d.variables[self.pblvar][curpos,miny:maxy,minx:maxx]
        for i in range(4):pblh+=np.float32(curdata[cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        # curdata=convolve(d.variables[self.LWuvar][miny:maxy,minx:maxx],gaussian,mode="same")
        # for i in range(4):longwave_up+=np.float32(curdata[cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        # curdata=convolve(d.variables[self.albvar][miny:maxy,minx:maxx],gaussian,mode="same")
        # for i in range(4):albedo+=np.float32(curdata[cury[:,:,i],curx[:,:,i]]*w[:,:,i])
        d.close()
        sensible_heat=sensible_heat[halfk:-halfk,halfk:-halfk]/6
        latent_heat=latent_heat[halfk:-halfk,halfk:-halfk]/6
        pblh=pblh[halfk:-halfk,halfk:-halfk]
        # print(sensible_heat.max())
        # print(latent_heat.max(),sensible_heat.max())
        return Bunch(sensible_heat=sensible_heat,latent_heat=latent_heat,pblh=pblh)
        # longwave_up=longwave_up[halfk:-halfk,halfk:-halfk]
        # albedo=albedo[halfk:-halfk,halfk:-halfk]
        # return Bunch(sensible_heat=-sensible_heat,latent_heat=-latent_heat,
        #              longwave_up=longwave_up,albedo=albedo/100.0)

    def nextBilin(self):
        curfile=self._curfile
        if curfile>=len(self._filenames):
            raise StopIteration
        x=self._geoLUT[:,:,:,1]
        y=self._geoLUT[:,:,:,0]
        w=self._geoLUT[:,:,:,2]
        xu=self._geoLUTu[:,:,:,1]
        yu=self._geoLUTu[:,:,:,0]
        wu=self._geoLUTu[:,:,:,2]
        xv=self._geoLUTv[:,:,:,1]
        yv=self._geoLUTv[:,:,:,0]
        wv=self._geoLUTv[:,:,:,2]
        halfk=self.halfk
        offset=np.int(np.ceil(halfk/self.io_ratio))
        minx=x.min().astype('i')-offset
        maxx=x.max().astype('i')+offset
        miny=y.min().astype('i')-offset
        maxy=y.max().astype('i')+offset
        curx=(x-minx).astype('i')
        cury=(y-miny).astype('i')
        minux= xu.min().astype('i')-offset
        maxux= xu.max().astype('i')+offset
        minuy= yu.min().astype('i')-offset
        maxuy= yu.max().astype('i')+offset
        curux=(xu-minux).astype('i')
        curuy=(yu-minuy).astype('i')
        minvx= xv.min().astype('i')-offset
        maxvx= xv.max().astype('i')+offset
        minvy= yv.min().astype('i')-offset
        maxvy= yv.max().astype('i')+offset
        curvx=(xv-minvx).astype('i')
        curvy=(yv-minvy).astype('i')

        d=swim_io.Dataset(self._filenames[curfile], 'r')
        
        datestr=''.join(d.variables[self.timevar][self.curpos])
        usepressures=self.usepressures
        
        Nz=usepressures.size
        Nxy=curx.shape[0:2]
        N=[Nz,Nxy[0],Nxy[1]]
        
        wind_u=np.zeros(N,dtype=np.float32,order="F")
        wind_v=np.zeros(N,dtype=np.float32,order="F")
        
        if self.windonly:
            curdata=fast_mean.fast_smooth(d.variables[self.uvar][self.curpos,...][:,minuy:maxuy,minux:maxux][usepressures,:,:],offset)
            for i in range(4):wind_u+=curdata[:,curuy[:,:,i],curux[:,:,i]]*wu[:,:,i]
            # curdata=d.variables[self.vvar][self.curpos,...][:,minvy:maxvy,minvx:maxvx][usepressures,:,:]
            curdata=fast_mean.fast_smooth(d.variables[self.vvar][self.curpos,...][:,minvy:maxvy,minvx:maxvx][usepressures,:,:],offset)
            for i in range(4):wind_v+=curdata[:,curvy[:,:,i],curvx[:,:,i]]*wv[:,:,i]
            # we also want geopotential height so we can properly interpolate vertically if necessary
            curdata=(d.variables[self.gphvar_pert][self.curpos,...][usepressures,:,:]+self.base_gph)[:,miny:maxy,minx:maxx]
            for i in range(4):hgt+=curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            hgt=hgt[:,halfk:-halfk,halfk:-halfk]/9.8
            potential_temperature=None
            qc=None
            specific_humidity=None
            pressure=None
            qr=None
            nr=None
            qs=None
            qi=None
            ni=None
            qg=None
        else:
            potential_temperature=np.zeros(N,dtype=np.float32,order="F")
            hgt=np.zeros(N,dtype=np.float32,order="F")
            qc=np.zeros(N,dtype=np.float32,order="F")
            # qi=np.zeros(N,dtype=np.float32) #for now qi is added to qc so we don't have to worry about ni
            specific_humidity=np.zeros(N,dtype=np.float32,order="F")
            pressure=np.zeros(N,dtype=np.float32,order="F")
            
            # Note, this currently reads all pressure levels from disk then subsets to "usepressures"
            curdata=d.variables[self.pvar_pert][self.curpos,...]
            curdata=(curdata[usepressures,:,:]+self.base_pressure)[:,miny:maxy,minx:maxx]
            for i in range(4):pressure+=curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            curdata=d.variables[self.tvar][self.curpos,...][:,miny:maxy,minx:maxx][usepressures,:,:]+300
            for i in range(4):potential_temperature+=curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            curdata=(d.variables[self.gphvar_pert][self.curpos,...][usepressures,:,:]+self.base_gph)[:,miny:maxy,minx:maxx]
            for i in range(4):hgt+=curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            curdata=d.variables[self.qcvar][self.curpos,...][:,miny:maxy,minx:maxx][usepressures,:,:]
            for i in range(4):qc+=curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            curdata=d.variables[self.qivar][self.curpos,...][:,miny:maxy,minx:maxx][usepressures,:,:]
            for i in range(4):qc+=curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            curdata=d.variables[self.qvvar][self.curpos,...][:,miny:maxy,minx:maxx][usepressures,:,:]
            for i in range(4):specific_humidity+=curdata[:,cury[:,:,i],curx[:,:,i]]*w[:,:,i]
            # curdata=d.variables[self.uvar][self.curpos,...][:,minuy:maxuy,minux:maxux][usepressures,:,:]
            curdata=fast_mean.fast_smooth(d.variables[self.uvar][self.curpos,...][:,minuy:maxuy,minux:maxux][usepressures,:,:],offset)
            for i in range(4):wind_u+=curdata[:,curuy[:,:,i],curux[:,:,i]]*wu[:,:,i]
            # curdata=d.variables[self.vvar][self.curpos,...][:,minvy:maxvy,minvx:maxvx][usepressures,:,:]
            curdata=fast_mean.fast_smooth(d.variables[self.vvar][self.curpos,...][:,minvy:maxvy,minvx:maxvx][usepressures,:,:],offset)
            for i in range(4):wind_v+=curdata[:,curvy[:,:,i],curvx[:,:,i]]*wv[:,:,i]

            potential_temperature=potential_temperature[:,halfk:-halfk,halfk:-halfk]
            pressure=pressure[:,halfk:-halfk,halfk:-halfk]
            specific_humidity=specific_humidity[:,halfk:-halfk,halfk:-halfk]
            hgt=hgt[:,halfk:-halfk,halfk:-halfk]/9.8
            qc=qc[:,halfk:-halfk,halfk:-halfk]
            
            N2=[Nz,Nxy[0]-halfk*2,Nxy[1]-halfk*2]
            qr=np.zeros(N2,dtype=np.float32,order="F")
            nr=np.zeros(N2,dtype=np.float32,order="F")
            qs=np.zeros(N2,dtype=np.float32,order="F")
            qi=np.zeros(N2,dtype=np.float32,order="F")
            ni=np.zeros(N2,dtype=np.float32,order="F")
            qg=np.zeros(N2,dtype=np.float32,order="F")
            
        
        self.curpos+=1
        kernelsize=np.int(np.floor(halfk))
        #        wind_u=wind_u[:,halfk:-halfk,halfk:-halfk]
        #        wind_v=wind_v[:,halfk:-halfk,halfk:-halfk]
        wind_u=fast_mean.fast_smooth(wind_u,halfk)[:,halfk:-halfk,halfk:-halfk]
        wind_v=fast_mean.fast_smooth(wind_v,halfk)[:,halfk:-halfk,halfk:-halfk]
        N=qc.shape
        #        print(wind_u.max(),wind_v.max())
        if self._sfc_files!=None:
            sfc=self._next_surface(self.curpos)
            sfc.pblh+=hgt[0,:,:]
            # convert the pblh into a vertical index into the 3D arrays
            pblindex=np.argmin(np.abs(sfc.pblh[np.newaxis,:,:]-hgt),axis=0)
            pblindex[pblindex<2]=2
            pblindex[pblindex>=17]=16
            sfc.pblh=pblindex+1
        else:
            sfc=None
        
        d.close()
        self.time_inc()
        return Bunch(p=pressure,th=potential_temperature,
                     sh=specific_humidity, hgt=hgt,qc=qc,
                     qv=specific_humidity/(1-specific_humidity), 
                     u=wind_u, v=wind_v,w=np.zeros(N2,dtype=np.float32,order="F"),
                     date=datestr,sfc=sfc,
                     qr=qr,qi=qi,qs=qs,qg=qg,ni=ni,nr=nr
                     )

    def time_inc(self):
        self.curpos+=1
        d=swim_io.Dataset(self._filenames[self._curfile], 'r')
        npos=d.variables[self.vvar].shape[0]-1
        if self.curpos>=npos:
            self._curfile+=1
            self.curpos=0
        d.close()
        
        
    def next(self):
        if self._geoLUT==None:
            return self.nextNN()
        else:
            return self.nextBilin()
            
    def close(self):
        pass
    def __enter__(self):
        return self
    
    def __exit__(self):
        self.close()
    
    def __del__(self):
        self.close()
