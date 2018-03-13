from molSimplify.Classes.mol3D import * #import the mol3D class
path = '/mol.xyz'
my_mol = mol3D() #Assign the mol3D class in local frame
my_mol.readfromxyz(path) # Look at specified XYZ File
metal_ind = my_mol.findMetal()[0] # Obtain index from metal list
bonded_atoms = my_mol.getBondedAtoms(metal_ind) #Check what's bonded to that metal
metal_coord = my_mol.getAtomCoords(metal_ind)
bondlengths = []
for i in bonded_atoms:
        connect_coord = my_mol.getAtomCoords(i)
        dist = distance(metal_coord,connect_coord)
        bondlengths.append(dist) #store in list
        print('The bond length between the metal and atom '+str(i)+' is '+str(dist)+' Angstrom.')
