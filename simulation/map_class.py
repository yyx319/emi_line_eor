#from DefineOptions import *
import numpy as np
#from icecream import ic
from matplotlib.colors import LogNorm
# Commencing to implement new map
class mapseed():
    ######################################################################################
    def __init__(self):
        print("Should this be a dataclass?")
        # Image properties
        self.pf=None           # Map values container
        self.complexmap=False  # Flag to determine whether map requires complex processing
        self.nchannel=0        # Number of information channels within map
        self.display_norm=None # Normalisation type for image display
        self.has_stream=False  # Determine whether this map has streamlines to be printed
        self.has_quiver=False  # Determine whether this map has quivers to be printed        
        
        # Complex images raw information
        self.valextrema=False  # Extrema values
        self.source_maps=None  # Source maps from which complex maps are generated
        self.stream_data=None  # Source map used to generate streamlines print
        self.quiver_data=None  # Source map used to generate quivers print
        print("TODO: Implement streamlines and quivers data")
        return
    ######################################################################################
    def set_map(self,pf):
        """ Sets the map value containers, and extracts additional map information """
        self.pf=pf
        print(pf.max())
        # Determine whether map is a multichannel image
        if (len(self.pf.shape)>2):
            self.nchannel=self.pf.shape[2]
            self.complexmap=True
        # Determine map extrema
        if (self.complexmap):
            minlist=[]; maxlist=[]
            for idim in range(self.nchannel):
                minlist.append(pf[:,:,idim].min())
                maxlist.append(pf[:,:,idim].max())
        else:
            minlist=[pf.min()]
            maxlist=[pf.max()]            
        self.pfmins=np.array(minlist)
        self.pfmaxs=np.array(maxlist)
        return
    ######################################################################################
    def square_map(self):
        """ Converts the map to a square shape by reducing the longest dimension """
        xdim_pix=self.pf.shape[0]
        ydim_pix=self.pf.shape[1]
        # If map already square, all done
        if (xdim_pix == ydim_pix): return
        cut_dim_y=False
        if (ydim_pix > xdim_pix): cut_dim_y=True
        cut_pix_left=int(0.5*abs(ydim_pix-xdim_pix))
        cut_pix_right=max(xdim_pix,ydim_pix)-(abs(ydim_pix-xdim_pix)-cut_pix_left)
        if (cut_dim_y):
            self.pf=self.pf[:,cut_pix_left:cut_pix_right]
        else:
            self.pf=self.pf[cut_pix_left:cut_pix_right]            
        return
    ######################################################################################
    def announce_extrema(self,pret):
        """ Announces the map seed max-min values """
        print(pret+"Map min=",self.pfmins,", max=",self.pfmaxs)
        if (self.valextrema): print(pret+"Map min=",self.valmins,", max=",self.valmaxs)
        return
    ######################################################################################
    def setup_RGB_map_from_maps(self,maplist):
        """ Combines multiple simple maps into a new complex map seed """
        logtype=type(LogNorm())
        self.nchannel=len(maplist)
        self.complexmap=True
        self.source_maps={}

        # Save information from source maps
        self.valextrema=True        
        minlist=[]; maxlist=[]; normlist=[]
        for imap, source_map in enumerate(maplist):
            minlist.append(np.nanmin(source_map.pf))
            maxlist.append(np.nanmax(source_map.pf))
            normlist.append(source_map.display_norm)
            self.source_maps[imap]=source_map
        self.valmins=np.array(minlist)
        self.valmaxs=np.array(maxlist)
        
        # Default RGB images ---------
        RGBmatrix="Identity"
        colscheme=[[255,0,0] ,[0,255,0],[0,0,255]]   # Red - Green - Blue print
        attenuation=[1.0,1.0,1.0]
        # WRB image for inflow - outflow -----
        #self.valmins[:]=1.e-27; self.valmaxs[:]=1.e-23
        self.valmins[:]=1.e-22; self.valmaxs[:]=5.e-18
        colscheme=[[150,150,150],[250,0,0],[0,0,250]]  # White - Red - Blue print
        

        # Superseed extrema if provided
        #if (len(plot_mivals)==self.nchannel): self.valmins=np.array(plot_mivals)
        #if (len(plot_mavals)==self.nchannel): self.valmaxs=np.array(plot_mavals)        
        if (RGBmatrix=="Identity"): RGBmatrix=[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
        colscheme=np.array(colscheme)
        RGBmatrix=np.array(RGBmatrix)
        attenuation=np.array(attenuation)

        # Generate channels data
        channels_data=[]
        for imap, source_map in enumerate(maplist):        
            data = source_map.pf
            # Normalise data according to received display norm
            if (type(source_map.display_norm) is logtype):
                # Limit negative values for log10
                data[data<0.0]=self.valmins[imap]
                data[data<0.0]=0.0
                # Convert to log10 values                
                data=np.log10(data)
                self.valmins[imap]=np.log10(self.valmins[imap])
                self.valmaxs[imap]=np.log10(self.valmaxs[imap])
            # Applying plot limits to map
            data[data!=data]=self.valmins[imap]
            data[data<self.valmins[imap]]=self.valmins[imap]
            data[data>self.valmaxs[imap]]=self.valmaxs[imap]
            # Normalise channel contribution to [0, 1]
            data = (data-self.valmins[imap])/(self.valmaxs[imap]-self.valmins[imap])            
            channels_data.append(data)

        # Process each channel according to colour and attenuation
        nx = channels_data[0].shape[0]
        ny = channels_data[0].shape[1]
        rgbArray = np.zeros((nx,ny,3), 'float32')
        for imap, channel_data in enumerate(channels_data):
            filter_values=RGBmatrix[imap][:]*attenuation[imap]
            for i in range(0,nx):
                for j in range(0,ny):
                    rgbArray[i][j][:]=rgbArray[i][j][:]+filter_values[:]*channel_data[i][j]

        # Colour scheme rotation + colour saturation cap
        self.pf = np.zeros((nx,ny,3), 'uint8')
        for i in range(0,nx):
            for j in range(0,ny):
                self.pf[i][j][0]=min(255,
                                     rgbArray[i][j][0]*colscheme[0][0]+
                                     rgbArray[i][j][1]*colscheme[1][0]+
                                     rgbArray[i][j][2]*colscheme[2][0])
                self.pf[i][j][1]=min(255,
                                     rgbArray[i][j][0]*colscheme[0][1]+
                                     rgbArray[i][j][1]*colscheme[1][1]+
                                     rgbArray[i][j][2]*colscheme[2][1])
                self.pf[i][j][2]=min(255,
                                     rgbArray[i][j][0]*colscheme[0][2]+
                                     rgbArray[i][j][1]*colscheme[1][2]+
                                     rgbArray[i][j][2]*colscheme[2][2])
            
        return
    ######################################################################################
