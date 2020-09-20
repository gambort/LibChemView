import matplotlib.patches as pa

import numpy as np

# Basic vector stuff
def length(x): return np.sqrt(np.sum(x**2))
def norm(x): return x/length(x)

# Basic matrix stuff
def Rot1(ang):
    c = np.cos(ang*np.pi/180.)
    s = np.sin(ang*np.pi/180.)
    return np.array([[1,0,0],[0,c,-s],[0,s,c]])

def Rot2(ang):
    c = np.cos(ang*np.pi/180.)
    s = np.sin(ang*np.pi/180.)
    return np.array([[c,0,s],[0,1,0],[-s,0,c]])

def Rot3(ang):
    c = np.cos(ang*np.pi/180.)
    s = np.sin(ang*np.pi/180.)
    return np.array([[c,-s,0],[s,c,0],[0,0,1.]])
       

# Color mixing
def CMix(c1,c2,alpha):
    a,b=alpha,1.-alpha
    return (c1[0]*a+c2[0]*b,
            c1[1]*a+c2[1]*b,
            c1[2]*a+c2[2]*b)

# Helper function
def zorder(z): return int(np.round(z*1000.))

# Chemical pallette system for pseudo-3D ball and stick
# models. May eventually be extended
class ChemPallette:
    def __init__(self,
                 XY0=(0.,0.),
                 ang1=0, ang2=0, ang3=0,
                 zDecay = 40., zOrigin=0.,
                 TextSize=12):
        # Initialise the pallette and (optionally)
        # set some defaults
        
        self.Patches=[]
        self.Text=[]

        self.xyz0 = np.array([XY0[0],XY0[1],0.])
        self.SetProjection(ang1, ang2, ang3)
        self.Scale = np.eye(3)
        
        self.Effects = {'zDecay': zDecay, 'zOrigin': zOrigin}

        self.TextSize = TextSize
    #######################################################
        
    def SetProjectionMatrix(self, M):
        self.Projection = M

    def SetProjection(self, ang1=0, ang2=0, ang3=0):
        # Define a projection for the system
        # must be done at the start

        M2 = Rot1(-ang1)
        M1 = Rot2(ang2)
        M3 = Rot3(ang3)
        
        self.Projection = M1.dot(M2).dot(M3)

    def ViewAxes(self, ID):
        if ID.lower() in ('xz' or 'zx'):
            self.SwapAxes('yz')
        if ID.lower() in ('yz' or 'zy'):
            self.SwapAxes('xz')

    def SwapAxes(self, ID):
        Indx={'x':0, 'y':1, 'z':2}
        I1 = Indx[ID[0].lower()]
        I2 = Indx[ID[1].lower()]

        if (I1==I2): return

        Swap = np.eye(3)
        Swap[I1,I2] = 1.
        Swap[I2,I1] =-1.
        Swap[I1,I1] = 0.
        Swap[I2,I2] = 0.

        self.Projection = np.dot(Swap, self.Projection)
        
    def TransformPoint(self, xyz):
        return np.dot(xyz,self.Projection).dot(self.Scale) + self.xyz0

    def SetEffects(self, EffectDict):
        for d in EffectDict:
            self.Effects[d] = EffectDict[d]

    #######################################################
    def Add(self, P,z,ID=""):
        # Add a patch
        self.Patches += [{'P':P,'z':z,'ID':ID}]

    def AddText(self,Label,xy,angle=0):
        # Add text
        self.Text += [{'Label':Label, 'xy':xy, 'angle':angle}]

    def Show(self, ax=None, alpha=None):
        # Show the patches in ax or gca()
        #
        # zDecay makes further atoms/bonds darker

        if ax is None: ax = plt.gca()

        zAll = np.array([x['z'] for x in self.Patches])
        ii = np.argsort(zAll)

        if self.Effects['zOrigin'] is None:
            zOrigin = zAll.mean()
        else: zOrigin = self.Effects['zOrigin']

        zDecay = self.Effects['zDecay']

        for i in ii:
            #print("%20s z = %6d"%(self.Patches[i]['ID'],
            #                      self.Patches[i]['z']))
            P = self.Patches[i]['P']
            if zDecay>0. or not(alpha is None):
                if alpha is None: alpha=1.
                # Darken the patch, if needed
                dz = self.Patches[i]['z'] - zOrigin
                beta = min(1., np.exp(dz/zorder(zDecay)))
                C = P.get_facecolor()
                P.set_facecolor((beta*C[0],beta*C[1],beta*C[2],alpha*C[3]))
            ax.add_patch(P)

        if self.TextSize>0:
            for T in self.Text:
                ax.text(T['xy'][0], T['xy'][1], T['Label'],
                        rotation=T['angle'],
                        horizontalalignment='center',
                        verticalalignment='center',
                        fontsize=self.TextSize)
            
    #######################################################
    def Atom(self, xyz_global, r=1,col=(1,0,0),
             ID="Atom",
             ShowHighlight=False):
        
        # Convert the atom position into the projection
        xyz = self.TransformPoint(xyz_global)
        
        sdk = 0.2
        slt = 0.1
        dkcol = CMix(col,(0,0,0),sdk)
        ltcol = CMix(col,(1,1,1),slt)

        ZP = zorder(xyz[2])

        xy0, r0 = xyz[:2], r
        if ShowHighlight:
            r1 = r*0.8
            d1 = 0.707
        else:
            r1 = r*0.7
            d1 = 0.5
        xy1 = [xyz[0]+(r0-r1)*d1,xyz[1]+(r0-r1)*d1]

        C0 = pa.Circle(xy0,r0,ec=None,fc=dkcol)
        C1 = pa.Circle(xy1,r1,ec=None,fc=col)

        self.Add(C0, ZP+0, ID = ID)
        self.Add(C1, ZP+1, ID = ID+"_hi")

        if ShowHighlight:
            r2, dr2 = r*0.2, r*0.4
            xy2 = [xyz[0]+dr2*0.707,xyz[1]+dr2*0.707]
            C2 = pa.Circle(xy2,r2,ec=None,fc=ltcol)
            self.Add(C2, ZP+2)


    def Bond(self,xyz_g1,xyz_g2,r1=1,r2=1,
             t=0.2,
             ID="Bond",
             Label=None,
             col=(0.5,0.5,0.5)):
        xyz1 = self.TransformPoint(xyz_g1)
        xyz2 = self.TransformPoint(xyz_g2)

        # dark color
        sdk = 0.2
        dkcol = CMix(col,(0,0,0),sdk)

        # Normal ray
        N = norm(np.array([1.,-1.,1.]))

        # Diagonal ray
        Q1 = norm(xyz2-xyz1)

        xyz1 += Q1*r1*0.95
        xyz2 -= Q1*r2*0.95

        #print("Bond from %8.3f %8.3f %8.3f -> %8.3f %8.3f %8.3f"\
        #      %(xyz1[0],xyz1[1],xyz1[2],xyz2[0],xyz2[1],xyz2[2]))

        NQ = norm(Q1 - np.dot(Q1,N)*N)

        S = 1.
        if Q1[0]<=0.: S*=-1.
        
        e1 = norm(Q1[:2])*S
        e2 = np.array([e1[1],-e1[0]])

        dw = e2*t/2.
        B1 = pa.Polygon([xyz1[:2]-dw, xyz1[:2]+dw,
                         xyz2[:2]+dw, xyz2[:2]-dw],
                        closed=True,
                        facecolor=dkcol,
                        linewidth=0,
        )
        
        zMiddle = (xyz2[2]+xyz1[2])/2.
        self.Add(B1, zorder(zMiddle), ID)

        # Add the highlight
        x = (1. - np.dot(N,NQ)**2)
        if x>0.1:
            dw1 =  e2*t/2.5
            dw2 =  e2*t/2.5*(1.-x)+t/3.0*e1*Q1[2]**2
            B2 = pa.Polygon([xyz1[:2]-dw1, xyz1[:2]+dw2,
                             xyz2[:2]+dw2, xyz2[:2]-dw1],
                            closed=True,
                            facecolor=col,
                            linewidth=0,
            )

            self.Add(B2, zorder(zMiddle)+1, ID+"_hi")


        # Add end caps
        if (Q1[2]>0.):
            w = t
            h = t*Q1[2]**2
            angle = -np.arctan2(Q1[0],Q1[1])*180./np.pi
            T = pa.Ellipse(xyz1[:2],w,h,angle=angle,
                           facecolor=dkcol)
            self.Add(T, zorder(xyz1[2])-1, ID+"_cap")
        else:
            w = t
            h = t*Q1[2]**2
            angle = -np.arctan2(Q1[0],Q1[1])*180./np.pi
            T = pa.Ellipse(xyz2[:2],w,h,angle=angle,
                           facecolor=dkcol)
            self.Add(T, zorder(xyz2[2])-1, ID+"_cap")

        if not(Label is None):
            # Add label if required
            xyc = (xyz2+xyz1)[:2]/2.
            if e2[1]>=0.:
                Lxy = xyc - (t+0.15)*e2
            else:
                Lxy = xyc + (t+0.15)*e2
            angle = -np.arctan2(e2[0],e2[1])*180./np.pi
            angle = ((angle+90.)%180.)-90.

            self.AddText(Label, Lxy, angle=angle)

if __name__=="__main__":
    import matplotlib.pyplot as plt
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,4))

    for ax,ang in ((ax1,0),(ax2,45)):
        CP = ChemPallette()
        CP.SetProjection(30,ang,0)
        X1,r1 = [ 0,0,0],0.3
        X2,r2 = [ 0,0,1],0.3
        X3,r3 = [ 0,1,0],0.3
        X4,r4 = [ 1,0,0],0.3
        
        CP.Atom(X1,r=r1,col=(0.5,0.5,0.5))
        CP.Atom(X2,r=r2,col=(0,0,1))
        CP.Atom(X3,r=r3,col=(0,1,0))
        CP.Atom(X4,r=r4,col=(1,0,0))
        CP.Bond(X1,X2,r1,r2)
        CP.Bond(X1,X3,r1,r3)
        CP.Bond(X1,X4,r1,r4)

        for Ang in (0.,60.,120.,180.,240.,300.):
            A1 = Ang
            A2 = Ang+60.
            c1,s1 = np.cos(A1/180.*np.pi), np.sin(A1/180.*np.pi)
            c2,s2 = np.cos(A2/180.*np.pi), np.sin(A2/180.*np.pi)

            X1, r1 = [3.*c1, 0., 3.*s1], 0.3
            X2, r2 = [3.*c2, 0., 3.*s2], 0.3
            CP.Atom(X1, r1, col=(1,1,0))
            CP.Bond(X1, X2, r1, r2, col=(1,0,1))

            X1, r1 = [6.*c1, 0., 6.*s1], 0.6
            X2, r2 = [6.*c2, 0., 6.*s2], 0.6
            CP.Atom(X1, r1, col=(1,1,0))
            CP.Bond(X1, X2, r1, r2, col=(1,0,1))
        CP.Show(ax)
    ax1.axis('equal')
    ax2.axis('equal')
    plt.show()
    
