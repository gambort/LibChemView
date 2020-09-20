from mplChemView.mplChemView import MoleculeDrawer
import numpy as np

if __name__=="__main__":
    import matplotlib.pyplot as plt
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,4))

    
    # Add LiH
    M1 = MoleculeDrawer()
    # Li at the origin
    M1.AddAtom('Li',[0,0,0])
    # H along the x axiis
    M1.AddAtom('H' ,[1.4,0,0])
    # Add the bonds
    M1.AddBonds()
    # And the bond label
    M1.LabelBond(0,1,"Li-H")

    # Add a mock LiH energy plot
    R = np.linspace(1.0,6)
    x = 2.*(R/1.4-1.)
    E = -5.*(1.+x)*np.exp(-x)
    X0 = 1.4
    Y0 = -0.5
    ax1.plot(R,E)

    # ME is e1,e2,e3 as a molecule
    # Set up a molecule drawer
    ME = MoleculeDrawer()
    # Add a C atom at the origin
    ME.AddAtom('C',[0,0,0])
    # Add an H along e1
    ME.AddAtom('H',[0.75,0,0])
    # Add an O along e2
    ME.AddAtom('O',[0,1.4,0])
    # Add an V along e3
    ME.AddAtom('V',[0,0,2.0])
    # Automatically add bonds using covalent radii
    ME.AddBonds()
    # Add some labels between molecules (indexed in order)
    ME.LabelBond(0,2,"C-O")
    # Add some labels between molecules (indexed in order)
    ME.LabelBond(0,3,"C-V")

    # Set up a molecule drawer
    M2 = MoleculeDrawer()
    # Read in an xyz file
    M2.ReadXYZFile("TestMolecule.xyz")
    # Embed the previous molecule offset by (2,2,2)
    M2.AddMolecule(ME, xyz=(2,2,2))

    
    # Show the two plots
    M1.DrawToAxes(ax1,
                  XY0=[X0,Y0], # Ensure molecule origin (Li) is here
                  ViewAng=(0,0,0), # Keep along the plane
                  ViewAxes='xz', # View the xz plane
    )
    M2.DrawToAxes(ax2,
                  ViewAng=(30,45,0), # Rotate to this angle
                  ViewAxes='yz', # View the yz plane
    )

    # Add M2 as an inset to ax1
    a3=M2.DrawToAxes(ax1, # Parent axis
                     # Show 50% along and 10% up with 40% width/height
                     Inset=(0.3,0.1,0.6,0.6),
                     # View along yz plane
                     ViewAxes='yz',
                     # Set to 50% transparency
                     alpha=0.5,
                     # Do not show bond labels
                     TextSize=0
    )

    # Don't use plt.tight_layout() because of inset
    # Save the file
    plt.savefig("TestPlot.png")
    # Show it
    plt.show()
    
