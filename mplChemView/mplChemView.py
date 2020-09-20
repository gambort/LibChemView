from mplChemView.mplChemViewBackEnd import ChemPallette
from NiceColours import *
import numpy as np
from matplotlib.transforms import Bbox

SpeciesColours = {
    'H': 'White',
    'AM': 'Lavender',
    'AEM': 'Green|Black',
    'Boron': 'Pink',
    'TM': 'Pink',
    'C': 'Black|Black|Grey',
    'N': 'Blue',
    'O': 'Red',
    'F': 'Green|Yellow',
    'Cl': 'Lime',
    'Br': 'Red|Black',
    'I': 'Lavender|Black',
    'NG': 'Cyan',
    'P': 'Orange',
    'S': 'Yellow',
    'Ti': 'Grey',
    'Cu': 'Brown',
    'Hg': 'Grey|White',
    'Unknown': 'Beige',
}

SpeciesClass = { # This needs to be filled in
    'AM': ('Li','Na','K','Rb','Cs','Fr'),
    'AEM': ('Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra'),
    'TM': ('Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn',
           'Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
           'La','Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg',
           'Ac','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cm'),
    'NG': ('He', 'Ne', 'Ar', 'Kr', 'Xe'),
}

SpeciesRadii = {
    'H': 0.38, 'He': 0.32, 'Li': 1.34, 'Be': 0.90, 'B': 0.82 ,
    'C': 0.77, 'N': 0.75, 'O': 0.73, 'F': 0.71, 'Ne': 0.69 ,
    'Na': 1.54, 'Mg': 1.30, 'Al': 1.18, 'Si': 1.11, 'P': 1.06 ,
    'S': 1.02, 'Cl': 0.99, 'Ar': 0.97, 'K': 1.96, 'Ca': 1.74 ,
    'Sc': 1.44, 'Ti': 1.36, 'V': 1.25, 'Cr': 1.27, 'Mn': 1.39 ,
    'Fe': 1.25, 'Co': 1.26, 'Ni': 1.21, 'Cu': 1.38, 'Zn': 1.31 ,
    'Ga': 1.26, 'Ge': 1.22, 'As': 1.19, 'Se': 1.16, 'Br': 1.14 ,
    'Kr': 1.10, 'Rb': 2.11, 'Sr': 1.92, 'Y': 1.62, 'Zr': 1.48 ,
    'Nb': 1.37, 'Mo': 1.45, 'Tc': 1.56, 'Ru': 1.26, 'Rh': 1.35 ,
    'Pd': 1.31, 'Ag': 1.53, 'Cd': 1.48, 'In': 1.44, 'Sn': 1.41 ,
    'Sb': 1.38, 'Te': 1.35, 'I': 1.33, 'Xe': 1.30, 'Cs': 2.25 ,
    'Ba': 1.98, 'La': 1.69, 'Ce': 1.65, 'Pr': 0.11, 'Nd': 0.20 ,
    'Pm': 1.50, 'Sm': 1.50, 'Eu': 1.50, 'Gd': 1.50, 'Tb': 1.50 ,
    'Dy': 1.50, 'Ho': 1.50, 'Er': 1.50, 'Tm': 1.50, 'Yb': 1.50 ,
    'Lu': 1.60, 'Hf': 1.50, 'Ta': 1.38, 'W': 1.46, 'Re': 1.59 ,
    'Os': 1.28, 'Ir': 1.37, 'Pt': 1.28, 'Au': 1.44, 'Hg': 1.49 ,
    'Tl': 1.48, 'Pb': 1.47, 'Bi': 1.46, 'Po': 1.50, 'At': 1.50 ,
    'Rn': 1.45, 'Fr': 1.50, 'Ra': 1.50, 'Ac': 1.50, 'Th': 1.50 ,
    'Pa': 1.50, 'U': 1.50, 'Np': 1.50, 'Pu': 1.50, 'Am': 1.50 ,
    'Cm': 1.50, 'Bk': 1.50, 'Cf': 1.50, 'Es': 1.50, 'Fm': 1.50 ,
    'Md': 1.50, 'No': 1.50, 'Lr': 1.50, 'Rf': 1.50, 'Db': 1.50 ,
    'Sg': 1.50, 'Bh': 1.50, 'Hs': 1.50, 'Mt': 1.50, 'Ds': 1.50 ,
    'Rg': 1.50, 'Cn': 1.50, 'Nh': 1.50, 'Fl': 1.50, 'Mc': 1.50 ,
    'Lv': 1.50, 'Ts': 1.50, 'Og': 1.50 ,
}

class MoleculeDrawer:
    def __init__(self, xyzUnits="Ang"):
        self.xyzScale = 1.0 # Internally works in Angstrom
        if xyzUnits.upper() in ("A0", "HA", "BOHR"):
            self.xyzScale = 0.53
        elif xyzUnits.upper() in ("PM"):
            self.xyzScale = 0.01
        elif xyzUnits.upper() in ("NM"):
            self.xyzScale = 10.

        self.Clear()

    def Clear(self):
        self.Atoms, self.Bonds = [], []
    def ClearAtoms(self): self.Atoms = []
    def ClearBonds(self): self.Bonds = []
        
    def AddAtom(self, Element, xyz,
                AtomRadius=0.3):
        # Add Element to position xyz
        
        # Get the element's colour
        if Element in SpeciesColours:
            ColDesc = SpeciesColours[Element]
        else:
            ColDesc = SpeciesColours['Unknown']
            for Class in SpeciesClass:
                if Element in SpeciesClass[Class]:
                    ColDesc = SpeciesColours[Class]

        # Average the colours
        CC = ColDesc.split('|')
        col = [0.,0.,0.]
        for k in range(len(CC)):
            TC = NiceColour(CC[k])
            col[0]+=TC[0]
            col[1]+=TC[1]
            col[2]+=TC[2]
        col[0]/=len(CC)
        col[1]/=len(CC)
        col[2]/=len(CC)
        col = tuple(col)


        # Get the element's bond information
        if Element in SpeciesRadii:
            BondRadius = SpeciesRadii[Element]
        else:
            BondRadius = 1.5
            
        # This should also have bond information
        self.Atoms += [{'xyz':np.array(xyz)*self.xyzScale,
                        'col':col,
                        'Element':Element,
                        'rA':AtomRadius,
                        'rB':BondRadius}]

    def AddBonds(self, fudge=0.05, thickness=0.1):
        NAtoms = len(self.Atoms)
        for K1 in range(NAtoms):
            for K2 in range(K1+1,NAtoms):
                A1 = self.Atoms[K1]
                A2 = self.Atoms[K2]
                xyz1 = A1['xyz']
                xyz2 = A2['xyz']

                R = np.sqrt(np.sum((xyz2-xyz1)**2))
                R12 = A1['rB']+A2['rB']

                if R<(R12*(1.+fudge)):
                    self.Bonds += [{'K12':(K1,K2),
                                    'xyz1':xyz1, 'xyz2':xyz2,
                                    'r1':A1['rA'], 'r2':A2['rA'],
                                    't':thickness,
                                    'Label':None,
                                    'col':(.2,.2,.2),
                    }]

    def RemoveAtom(self, K):
        if K<len(self.Atoms):
            self.Atoms.pop(K)

    def FindBond(self, K1,K2):
        BList = []
        for KB in range(len(self.Bonds)):
            B = self.Bonds[KB]
            if B['K12']==(K1,K2) or B['K12']==(K2,K1):
                BList += [KB]
        return BList
        
    def RemoveBond(self, K1,K2):
        Remove = self.FindBond(K1,K2)
        for R in Remove:
            self.Bonds.pop(R)
            
    def LabelBond(self, K1,K2,Label=None):
        BList = self.FindBond(K1,K2)
        if Label is None:
            Label = "%s-%s"%(self.Atoms[K1]['Element'],
                             self.Atoms[K2]['Element'])

        for kB in BList:
            self.Bonds[kB]['Label']=Label

    def ReadXYZFile(self, FileName):
        F = open(FileName)
        N = int(F.readline())
        Comment = F.readline()
        for k in range(N):
            try:
                T = F.readline().split()
            except:
                print("File does not contain enough entries")
                break
            El = T[0]
            xyz = np.array([float(x) for x in T[1:4]])
            self.AddAtom(El, xyz)
        F.close()
        self.AddBonds()

    def AddMolecule(self, M, xyz=[0,0,0]):
        for A in M.Atoms:
            AC = dict(A)
            AC['xyz']=A['xyz']+np.array(xyz)
            self.Atoms += [AC]
            
        for B in M.Bonds:
            BC = dict(B)
            BC['xyz1']=B['xyz1']+np.array(xyz)
            BC['xyz2']=B['xyz2']+np.array(xyz)
            self.Bonds += [BC]
        
    def DrawToPallette(self, ChemPallette):
        for A in self.Atoms:
            xyz = A['xyz']
            col = A['col']
            r = A['rA']
            ChemPallette.Atom(xyz, r, col=col)

        for B in self.Bonds:
            xyz1 = B['xyz1']
            xyz2 = B['xyz2']
            r1 = B['r1']
            r2 = B['r2']
            col = B['col']
            t = B['t']
            Label = B['Label']

            ChemPallette.Bond(xyz1, xyz2, r1, r2,
                              t=t, col=col, Label=Label)

    def DrawToAxes(self, ax=None,
                   Inset=None,
                   XY0=(0,0),
                   ViewAng=(30,45,0),
                   ViewAxes='xy',
                   TextSize=12,
                   alpha=1.,
    ):
        if not(Inset is None):
            # add as a new inset
            fig = ax.figure
            Coords = Bbox.from_bounds(Inset[0], Inset[1], Inset[2], Inset[3])
            AxCoords = ax.transAxes.transform(Coords)
            FigCoords = fig.transFigure.inverted().transform(AxCoords)
            
            ax2 = fig.add_axes(Bbox(FigCoords))
            ax = ax2

        CP = ChemPallette(XY0=XY0, TextSize=TextSize)
        CP.SetProjection(ViewAng[0],ViewAng[1],ViewAng[2])
        CP.ViewAxes(ViewAxes)

        self.DrawToPallette(CP)
        CP.Show(ax, alpha=alpha)

        ax.axis('equal')
        if not(Inset is None):
            ax.axis('off')

        return ax
    
if __name__=="__main__":
    import matplotlib.pyplot as plt
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,4))

    #OutputNiceColourTable()

    # ME is e1,e2,e3 as a molecule
    ME = MoleculeDrawer()
    ME.AddAtom('C',[0,0,0])
    ME.AddAtom('H',[0.75,0,0])
    ME.AddAtom('O',[0,1.4,0])
    ME.AddAtom('V',[0,0,2.0])
    ME.AddBonds() # Fills in bonds using covalent radii
    ME.LabelBond(0,2,"C-O")
    ME.LabelBond(0,3,"C-V")

    # A real molecule
    M2 = MoleculeDrawer()
    M2.ReadXYZFile("TestMolecule.xyz")
    # And embed ME
    M2.AddMolecule(ME, xyz=(2,2,2))

    # Add a LiH plot to go with the molecule
    R = np.linspace(1.0,6)
    x = 2.*(R/1.4-1.)
    E = -4.*(1.+x)*np.exp(-x)
    X0 = 1.4
    Y0 = -0.5
    ax1.plot(R,E)
    
    # Add LiH
    M1 = MoleculeDrawer()
    M1.AddAtom('Li',[0,0,0])
    M1.AddAtom('H' ,[1.4,0,0])
    M1.AddBonds()
    M1.LabelBond(0,1,"Li-H")
    
    for ax,M in ((ax1,M1), (ax2,M2)):
        if M==M1: ang=(0,0,0)
        else: ang=(30,45,0)
        M.DrawToAxes(ax, XY0=[X0,Y0], ViewAng=ang)

        ax.axis('equal')

    M2.DrawToAxes(ax1, Inset=(0.5,0.1,0.4,0.4))
    #plt.tight_layout()
    plt.savefig("TestPlot.png")
    plt.show()
    
