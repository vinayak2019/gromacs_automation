from ocelot.task.readcif import ReadCif
import numpy as np
import sys

__author__ = 'Vinayak Bhat'


def make_gro(location,code_name):

    # read the cif file and select a clean structure
    rc = ReadCif.from_ciffile(location,"")
    rc.read()
    config = rc.results["configurations"][0]
    pstructure = config.pstructure

    # get the molecules
    mols = config.molconformers

    # get number of atoms in unit cell and in a molecule
    num_atoms = len(pstructure.sites)
    mol_atoms = len(mols[0].pmgmol)

    # lattice parameters for periodic box
    flat_lattice = list(np.concatenate(list(pstructure.lattice.as_dict()["matrix"])).flat)
    lattice_order = [0,4,8,1,2,3,5,6,7]
    flat_lattice = [flat_lattice[i] for i in lattice_order]

    # write the periodic and molecular gro files
    periodic_file = open("{}.gro".format(code_name),"w")
    molecule_file = open("{}_mol.gro".format(code_name), "w")

    # first two lines of the gro file
    atoms_counter = 1
    print(code_name,file=periodic_file)
    print(str(num_atoms), file=periodic_file)
    print(code_name,file=molecule_file)
    print(str(mol_atoms), file=molecule_file)

    # adding the atomic coordinates
    for mol_idx, mol in enumerate(mols):
        for site in mol.sites:
            element_name = site.specie.name
            coord_string = "{}{}{}".format(format(site.coords[0] /10,"8.3f"), format(site.coords[1]/10,"8.3f"), format(site.coords[2]/10,"8.3f"))
            res_name = str(mol_idx + 1) + code_name
            atom_number = str(atoms_counter)
            site_string = "{0:>10s}{1:>5s}{2:>5s}{3}".format(res_name, element_name, atom_number, coord_string)
            print(site_string,file=periodic_file)
            if mol_idx == 0:
                print(site_string, file=molecule_file)
            atoms_counter += 1

    # add the box dimensions
    lattice_string = ""
    for lattice in flat_lattice:
        lattice_string += format(lattice/10,'.5f') + " "
    print(lattice_string,file=periodic_file)
    print("5.0 5.0 5.0", file=molecule_file)

    # close the files
    periodic_file.close()
    molecule_file.close()

if __name__ == "__main__":
    location = sys.argv[1]
    code_name = sys.argv[2]
    make_gro(location,code_name)

