#!/usr/bin/env python3

#******************************************
#
#    SHARC Program Suite
#
#    SHARC-MN Extension
#
#    Copyright (c) 2023 University of Vienna
#    Copyright (c) 2025 University of Minnesota
#
#    This file is part of SHARC-MN.
#
#    SHARC is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SHARC is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    inside the SHARC manual.  If not, see <http://www.gnu.org/licenses/>.
#
#******************************************

# This script is a wrapper for sampling initial conditions for bimolecular processes. 
# This script will use state_selected.py script 
# by Yinan Shu, Jan. 19, 2025. 

import os
import subprocess
from optparse import OptionParser
import random 
import datetime
import math
import sys

# some constants
DEBUG = False
HARTREE_TO_EV = 27.211396132    # conversion factor from Hartree to eV
U_TO_AMU = 1./5.4857990943e-4            # conversion from g/mol to amu
ANG_TO_BOHR = 1./0.529177211    #1.889725989      # conversion from Angstrom to bohr

version = '1.0'
versiondate = datetime.date(2025, 1, 19)

NUMBERS = {'H':1, 'He':2,
'Li':3, 'Be':4, 'B':5, 'C':6,  'N':7,  'O':8, 'F':9, 'Ne':10,
'Na':11, 'Mg':12, 'Al':13, 'Si':14,  'P':15,  'S':16, 'Cl':17, 'Ar':18,
'K':19, 'Ca':20,
'Sc':21, 'Ti':22, 'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30,
'Ga':31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36,
'Rb':37, 'Sr':38,
'Y':39,  'Zr':40, 'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48,
'In':49, 'Sn':50, 'Sb':51, 'Te':52,  'I':53, 'Xe':54,
'Cs':55, 'Ba':56,
'La':57,
'Ce':58, 'Pr':59, 'Nd':60, 'Pm':61, 'Sm':62, 'Eu':63, 'Gd':64, 'Tb':65, 'Dy':66, 'Ho':67, 'Er':68, 'Tm':69, 'Yb':70, 'Lu':71,
         'Hf':72, 'Ta':73,  'W':74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80,
'Tl':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86,
'Fr':87, 'Ra':88,
'Ac':89,
'Th':90, 'Pa':91,  'U':92, 'Np':93, 'Pu':94, 'Am':95, 'Cm':96, 'Bk':97, 'Cf':98, 'Es':99,'Fm':100,'Md':101,'No':102,'Lr':103,
        'Rf':104,'Db':105,'Sg':106,'Bh':107,'Hs':108,'Mt':109,'Ds':110,'Rg':111,'Cn':112,
'Nh':113,'Fl':114,'Mc':115,'Lv':116,'Ts':117,'Og':118
}

# Atomic Weights of the most common isotopes
# From https://chemistry.sciences.ncsu.edu/msf/pdf/IsotopicMass_NaturalAbundance.pdf
MASSES = {'H' :   1.007825 * U_TO_AMU,
          'He':   4.002603 * U_TO_AMU,
          'Li':   7.016004 * U_TO_AMU,
          'Be':   9.012182 * U_TO_AMU,
          'B' :  11.009305 * U_TO_AMU,
          'C' :  12.000000 * U_TO_AMU,
          'N' :  14.003074 * U_TO_AMU,
          'O' :  15.994915 * U_TO_AMU,
          'F' :  18.998403 * U_TO_AMU,
          'Ne':  19.992440 * U_TO_AMU,
          'Na':  22.989770 * U_TO_AMU,
          'Mg':  23.985042 * U_TO_AMU,
          'Al':  26.981538 * U_TO_AMU,
          'Si':  27.976927 * U_TO_AMU,
          'P' :  30.973762 * U_TO_AMU,
          'S' :  31.972071 * U_TO_AMU,
          'Cl':  34.968853 * U_TO_AMU,
          'Ar':  39.962383 * U_TO_AMU,
          'K' :  38.963707 * U_TO_AMU,
          'Ca':  39.962591 * U_TO_AMU,
          'Sc':  44.955910 * U_TO_AMU,
          'Ti':  47.947947 * U_TO_AMU,
          'V' :  50.943964 * U_TO_AMU,
          'Cr':  51.940512 * U_TO_AMU,
          'Mn':  54.938050 * U_TO_AMU,
          'Fe':  55.934942 * U_TO_AMU,
          'Co':  58.933200 * U_TO_AMU,
          'Ni':  57.935348 * U_TO_AMU,
          'Cu':  62.929601 * U_TO_AMU,
          'Zn':  63.929147 * U_TO_AMU,
          'Ga':  68.925581 * U_TO_AMU,
          'Ge':  73.921178 * U_TO_AMU,
          'As':  74.921596 * U_TO_AMU,
          'Se':  79.916522 * U_TO_AMU,
          'Br':  78.918338 * U_TO_AMU,
          'Kr':  83.911507 * U_TO_AMU,
          'Rb':  84.911789 * U_TO_AMU,
          'Sr':  87.905614 * U_TO_AMU,
          'Y' :  88.905848 * U_TO_AMU,
          'Zr':  89.904704 * U_TO_AMU,
          'Nb':  92.906378 * U_TO_AMU,
          'Mo':  97.905408 * U_TO_AMU,
          'Tc':  98.907216 * U_TO_AMU,
          'Ru': 101.904350 * U_TO_AMU,
          'Rh': 102.905504 * U_TO_AMU,
          'Pd': 105.903483 * U_TO_AMU,
          'Ag': 106.905093 * U_TO_AMU,
          'Cd': 113.903358 * U_TO_AMU,
          'In': 114.903878 * U_TO_AMU,
          'Sn': 119.902197 * U_TO_AMU,
          'Sb': 120.903818 * U_TO_AMU,
          'Te': 129.906223 * U_TO_AMU,
          'I' : 126.904468 * U_TO_AMU,
          'Xe': 131.904154 * U_TO_AMU,
          'Cs': 132.905447 * U_TO_AMU,
          'Ba': 137.905241 * U_TO_AMU,
          'La': 138.906348 * U_TO_AMU,
          'Ce': 139.905435 * U_TO_AMU,
          'Pr': 140.907648 * U_TO_AMU,
          'Nd': 141.907719 * U_TO_AMU,
          'Pm': 144.912744 * U_TO_AMU,
          'Sm': 151.919729 * U_TO_AMU,
          'Eu': 152.921227 * U_TO_AMU,
          'Gd': 157.924101 * U_TO_AMU,
          'Tb': 158.925343 * U_TO_AMU,
          'Dy': 163.929171 * U_TO_AMU,
          'Ho': 164.930319 * U_TO_AMU,
          'Er': 165.930290 * U_TO_AMU,
          'Tm': 168.934211 * U_TO_AMU,
          'Yb': 173.938858 * U_TO_AMU,
          'Lu': 174.940768 * U_TO_AMU,
          'Hf': 179.946549 * U_TO_AMU,
          'Ta': 180.947996 * U_TO_AMU,
          'W' : 183.950933 * U_TO_AMU,
          'Re': 186.955751 * U_TO_AMU,
          'Os': 191.961479 * U_TO_AMU,
          'Ir': 192.962924 * U_TO_AMU,
          'Pt': 194.964774 * U_TO_AMU,
          'Au': 196.966552 * U_TO_AMU,
          'Hg': 201.970626 * U_TO_AMU,
          'Tl': 204.974412 * U_TO_AMU,
          'Pb': 207.976636 * U_TO_AMU,
          'Bi': 208.980383 * U_TO_AMU,
          'Po': 208.982416 * U_TO_AMU,
          'At': 209.987131 * U_TO_AMU,
          'Rn': 222.017570 * U_TO_AMU,
          'Fr': 223.019731 * U_TO_AMU,
          'Ra': 226.025403 * U_TO_AMU,
          'Ac': 227.027747 * U_TO_AMU,
          'Th': 232.038050 * U_TO_AMU,
          'Pa': 231.035879 * U_TO_AMU,
          'U' : 238.050783 * U_TO_AMU,
          'Np': 237.048167 * U_TO_AMU,
          'Pu': 244.064198 * U_TO_AMU,
          'Am': 243.061373 * U_TO_AMU,
          'Cm': 247.070347 * U_TO_AMU,
          'Bk': 247.070299 * U_TO_AMU,
          'Cf': 251.079580 * U_TO_AMU,
          'Es': 252.082972 * U_TO_AMU,
          'Fm': 257.095099 * U_TO_AMU,
          'Md': 258.098425 * U_TO_AMU,
          'No': 259.101024 * U_TO_AMU,
          'Lr': 262.109692 * U_TO_AMU,
          'Rf': 267. * U_TO_AMU,
          'Db': 268. * U_TO_AMU,
          'Sg': 269. * U_TO_AMU,
          'Bh': 270. * U_TO_AMU,
          'Hs': 270. * U_TO_AMU,
          'Mt': 278. * U_TO_AMU,
          'Ds': 281. * U_TO_AMU,
          'Rg': 282. * U_TO_AMU,
          'Cn': 285. * U_TO_AMU,
          'Nh': 286. * U_TO_AMU,
          'Fl': 289. * U_TO_AMU,
          'Mc': 290. * U_TO_AMU,
          'Lv': 293. * U_TO_AMU,
          'Ts': 294. * U_TO_AMU,
          'Og': 294. * U_TO_AMU
}

# Isotopes used for the masses
ISOTOPES={'H' : 'H-1' ,
          'He': 'He-4',
          'Li': 'Li-7',
          'Be': 'Be-9*',
          'B' : 'B_11' ,
          'C' : 'C-12' ,
          'N' : 'N-14' ,
          'O' : 'O-16' ,
          'F' : 'F-19*' ,
          'Ne': 'Ne-20',
          'Na': 'Na-23*',
          'Mg': 'Mg-24',
          'Al': 'Al-27*',
          'Si': 'Si-28',
          'P' : 'P-31*' ,
          'S' : 'S-32' ,
          'Cl': 'Cl-35',
          'Ar': 'Ar-40',
          'K' : 'K-39' ,
          'Ca': 'Ca-40',
          'Sc': 'Sc-45*',
          'Ti': 'Ti-48',
          'V' : 'V-51' ,
          'Cr': 'Cr-52',
          'Mn': 'Mn-55*',
          'Fe': 'Fe-56',
          'Co': 'Co-59*',
          'Ni': 'Ni-58',
          'Cu': 'Cu-63',
          'Zn': 'Zn-64',
          'Ga': 'Ga-69',
          'Ge': 'Ge-74',
          'As': 'As-75*',
          'Se': 'Se-80',
          'Br': 'Br-79',
          'Kr': 'Kr-84',
          'Rb': 'Rb-85',
          'Sr': 'Sr-88',
          'Y' : 'Y-89*' ,
          'Zr': 'Zr-90',
          'Nb': 'Nb-93*',
          'Mo': 'Mo-98',
          'Tc': 'Tc-98',
          'Ru': 'Ru-102',
          'Rh': 'Rh-103*',
          'Pd': 'Pd-106',
          'Ag': 'Ag-107',
          'Cd': 'Cd-114',
          'In': 'In-115',
          'Sn': 'Sn-120',
          'Sb': 'Sb-121',
          'Te': 'Te-130',
          'I' : 'I-127*' ,
          'Xe': 'Xe-132',
          'Cs': 'Cs-133*',
          'Ba': 'Ba-138',
          'La': 'La-139',
          'Ce': 'Ce-140',
          'Pr': 'Pr-141*',
          'Nd': 'Nd-142',
          'Pm': 'Pm-145',
          'Sm': 'Sm-152',
          'Eu': 'Eu-153',
          'Gd': 'Gd-158',
          'Tb': 'Tb-159*',
          'Dy': 'Dy-164',
          'Ho': 'Ho-165*',
          'Er': 'Er-166',
          'Tm': 'Tm-169*',
          'Yb': 'Yb-174',
          'Lu': 'Lu-175',
          'Hf': 'Hf-180',
          'Ta': 'Ta-181',
          'W' : 'W-184' ,
          'Re': 'Re-187',
          'Os': 'Os-192',
          'Ir': 'Ir-193',
          'Pt': 'Pt-195',
          'Au': 'Au-197*',
          'Hg': 'Hg-202',
          'Tl': 'Tl-205',
          'Pb': 'Pb-208',
          'Bi': 'Bi-209*',
          'Po': 'Po-209',
          'At': 'At-210',
          'Rn': 'Rn-222',
          'Fr': 'Fr-223',
          'Ra': 'Ra-226',
          'Ac': 'Ac-227',
          'Th': 'Th-232*',
          'Pa': 'Pa-231*',
          'U' : 'U-238' ,
          'Np': 'Np-237',
          'Pu': 'Pu-244',
          'Am': 'Am-243',
          'Cm': 'Cm-247',
          'Bk': 'Bk-247',
          'Cf': 'Cf-251',
          'Es': 'Es-252',
          'Fm': 'Fm-257',
          'Md': 'Md-258',
          'No': 'No-259',
          'Lr': 'Lr-262',
              'Rf': 'Rf-267',
              'Db': 'Db-268',
              'Sg': 'Sg-269',
              'Bh': 'Bh-270',
              'Hs': 'Hs-270',
              'Mt': 'Mt-278',
              'Ds': 'Ds-281',
              'Rg': 'Rg-282',
              'Cn': 'Cn-285',
              'Nh': 'Nh-286',
              'Fl': 'Fl-289',
              'Mc': 'Mc-290',
              'Lv': 'Lv-293',
              'Ts': 'Ts-294',
              'Og': 'Og-294'
}


def write_xyz_coordinates(combined_data, output_xyz_file):
    """
    Writes the x, y, z coordinates for each initial condition to a separate file in the specified format.

    Args:
        combined_data (dict): Combined data structure.
        output_xyz_file (str): Path to the output XYZ file.
    """
    with open(output_xyz_file, "w") as f:
        for condition in combined_data["conditions"]:
            # Number of atoms
            num_atoms = len(condition["atoms"])
            f.write(f"{num_atoms}\n")
            f.write(f"{condition['index']}\n")
            for atom in condition["atoms"]:
                f.write(f"{atom['element']} {atom['x']/ANG_TO_BOHR:.6f} {atom['y']/ANG_TO_BOHR:.6f} {atom['z']/ANG_TO_BOHR:.6f}\n")
    print(f"XYZ coordinates written to {output_xyz_file}.")


def write_combined_data(combined_data, output_file):
    """
    Writes the combined data to a file in the same format as the initconds.

    Args:
        combined_data (dict): Combined data structure.
        output_file (str): Path to the output file.
    """
    with open(output_file, "w") as f:
        # Write the header
        header = combined_data["header"]
        f.write(f"Ninit     {header['Ninit']}\n")
        f.write(f"Natom     {header['Natom']}\n")
        f.write(f"Repr      {header['Repr']}\n")
        f.write(f"Temp      {header['Temp']:.10f}\n")
        f.write(f"Eref      {header['Eref']:.10f}\n")
        f.write(f"Eharm     {header['Eharm']:.10f}\n")
        f.write("\n")

        # Write the equilibrium geometry
        f.write("Equilibrium\n")
        for atom in combined_data["equilibrium"]:
            f.write(f"{atom['element']}   {atom['atomic_number']:.1f}   {atom['x']:.8f}   {atom['y']:.8f}   {atom['z']:.8f}   "
                    f"{atom['mass']:.8f}   {atom['vx']:.8f}   {atom['vy']:.8f}   {atom['vz']:.8f}\n")
        f.write("\n")

        # Write initial conditions
        for condition in combined_data["conditions"]:
            f.write(f"Index     {condition['index']}\n")
            f.write("Atoms\n")
            for atom in condition["atoms"]:
                f.write(f"{atom['element']}   {atom['atomic_number']:.1f}   {atom['x']:.8f}   {atom['y']:.8f}   {atom['z']:.8f}   "
                        f"{atom['mass']:.8f}   {atom['vx']:.8f}   {atom['vy']:.8f}   {atom['vz']:.8f}\n")
            f.write("States\n")
            for state, values in condition["states"].items():
                f.write(f"{state:<12} {values['value_au']:.12f} a.u.    {values['value_ev']:.12f} eV\n")
            f.write("\n")
    print(f"Combined data written to {output_file}.")


def combine_molecule_data(molecule1_data, molecule2_data, E_col):
    """
    Combines molecule1_data and molecule2_data into a single data structure for merged molecules.

    Args:
        molecule1_data (dict): Data structure for molecule 1.
        molecule2_data (dict): Data structure for molecule 2.

    Returns:
        dict: Combined data structure.
    """
    # Validate that both molecules have the same number of initial conditions
    if len(molecule1_data["conditions"]) != len(molecule2_data["conditions"]):
        raise ValueError("Molecule 1 and Molecule 2 must have the same number of initial conditions.")

    # Combine headers
    combined_data = {
        "header": {
            "Ninit": molecule1_data["header"]["Ninit"],
            "Natom": molecule1_data["header"]["Natom"] + molecule2_data["header"]["Natom"],
            "Repr": None,
            "Temp": 0.0,  # Can be updated based on requirements
            "Eref": 0.0,  # Can be updated based on requirements
            "Eharm": 0.0,  # Can be updated based on requirements
        },
        "equilibrium": molecule1_data["equilibrium"] + molecule2_data["equilibrium"],
        "conditions": []
    }

    # Combine conditions
    for cond1, cond2 in zip(molecule1_data["conditions"], molecule2_data["conditions"]):
        combined_condition = {
            "index": cond1["index"],
            "atoms": cond1["atoms"] + cond2["atoms"],  # Combine atoms
            "states": {
                "Ekin": {
                    "value_au": cond1["states"]["Ekin"]["value_au"] + cond2["states"]["Ekin"]["value_au"] + E_col,
                    "value_ev": cond1["states"]["Ekin"]["value_ev"] + cond2["states"]["Ekin"]["value_ev"] + E_col*27.211396132,
                },
                "Epot": {
                    "value_au": cond1["states"]["Epot"]["value_au"] + cond2["states"]["Epot"]["value_au"],
                    "value_ev": cond1["states"]["Epot"]["value_ev"] + cond2["states"]["Epot"]["value_ev"],
                },
                "Etot": {
                    "value_au": cond1["states"]["Etot"]["value_au"] + cond2["states"]["Etot"]["value_au"] + E_col,
                    "value_ev": cond1["states"]["Etot"]["value_ev"] + cond2["states"]["Etot"]["value_ev"] + E_col*27.211396132,
                },
            }
        }
        combined_data["conditions"].append(combined_condition)

    print("Molecule data combined successfully.")
    return combined_data


def add_collision_velocity(molecule2_data, E_col):
    """
    Adds a uniform velocity to each atom of molecule 2 in the -z direction based on the collision energy.

    Args:
        molecule2_data (dict): Data structure for molecule 2.
        E_col (float): Collision energy in the same units as the masses in molecule2_data.
    """
    # Compute the total mass of molecule 2
    total_mass = sum(atom["mass"] for atom in molecule2_data["equilibrium"])

    if total_mass == 0:
        raise ValueError("Total mass of molecule 2 is zero. Check molecule2_data for missing atom masses.")

    # Compute the velocity
    velocity = math.sqrt(2 * E_col / total_mass)

    print(f"Computed uniform velocity for molecule 2 atoms: {velocity:.6f} bohr/(atomic_time)")

    # Add velocity to each atom in each initial condition
    for condition in molecule2_data["conditions"]:
        for atom in condition["atoms"]:
            atom["vz"] -= velocity  # Add the velocity in the -z direction

    print("Collision velocity added to molecule 2 atoms.")


def generate_atom_data(n, atom, position, z_x_values=None):
    """
    Generates initial conditions for a single atom.

    Args:
        n (int): Number of initial conditions.
        atom (str): Symbol of the atom (e.g., 'H' for hydrogen).
        position (int): 1 for placing atom at origin, 2 for placing atom at z_x_values.
        z_x_values (list): List of tuples [(x, y), ...] specifying atom positions for each initial condition if position=2.

    Returns:
        dict: A data structure similar to molecule1_data/molecule2_data.
    """
    if atom not in MASSES:
        raise ValueError(f"Atom symbol '{atom}' is not in the MASSES dictionary.")

    if position == 2 and (z_x_values is None or len(z_x_values) != n):
        raise ValueError("z_x_values must be provided and match the number of initial conditions when position=2.")

    atom_mass = MASSES[atom]/U_TO_AMU
    atom_number = NUMBERS[atom]

    # Generate the atom data structure
    atom_data = {
        "header": {
            "Ninit": n,
            "Natom": 1,
            "Repr": None,
            "Temp": 0.0,
            "Eref": 0.0,
            "Eharm": 0.0,
        },
        "equilibrium": [
            {
                "element": atom,
                "atomic_number": float(atom_number),
                "x": 0.0 if position == 1 else z_x_values[0][1],
                "y": 0.0,
                "z": 0.0 if position == 1 else z_x_values[0][0],
                "mass": atom_mass,
                "vx": 0.0,
                "vy": 0.0,
                "vz": 0.0,
            }
        ],
        "conditions": []
    }

    # Add initial conditions
    for i in range(n):
        atom_data["conditions"].append({
            "index": i + 1,
            "atoms": [
                {
                    "element": atom,
                    "atomic_number": float(atom_number),
                    "x": 0.0 if position == 1 else z_x_values[i][1],
                    "y": 0.0,
                    "z": 0.0 if position == 1 else z_x_values[i][0],
                    "mass": atom_mass,
                    "vx": 0.0,
                    "vy": 0.0,
                    "vz": 0.0,
                }
            ],
            "states": {
                "Ekin": {"value_au": 0.0, "value_ev": 0.0},
                "Epot": {"value_au": 0.0, "value_ev": 0.0},
                "Etot": {"value_au": 0.0, "value_ev": 0.0},
            },
        })

    return atom_data


def sample_z_x(n, bmin, bmax, strata, separation):
    """
    Samples z and x values for each initial condition.
    Args:
        n (int): Number of initial conditions.
        bmin (float): Minimum value of x.
        bmax (float): Maximum value of x.
        strata (int): Number of strata for x sampling.
        separation (float): Value of a, fixed for all initial conditions.
    Returns:
        list: A list of tuples (z,x) for each initial condition.
    """
    z = separation
    z_x_values = []

    for istrata in range(strata):
        indexmin = istrata * int(n / strata)
        indexmax = (istrata + 1) * int(n / strata)
        if istrata==strata-1:
            indexmax = n
        for i in range(indexmin,indexmax):
            bmin_stratum = bmin + istrata * (bmax - bmin) / strata
            bmax_stratum = bmin + (istrata + 1) * (bmax - bmin) / strata
            x = random.uniform(bmin_stratum, bmax_stratum)
            z_x_values.append((z, x))

    return z_x_values


def move_to_origin(conditions):
    """
    Moves the center of mass of the molecule to the origin (0, 0, 0).
    Args:
        conditions (list): List of initial conditions, where each condition contains atom data.
    """
    for condition in conditions:
        current_com = compute_center_of_mass([condition])[0]
        shift_x = current_com["com_x"]
        shift_y = current_com["com_y"]
        shift_z = current_com["com_z"]

        for atom in condition["atoms"]:
            atom["x"] -= shift_x
            atom["y"] -= shift_y
            atom["z"] -= shift_z


def compute_center_of_mass(conditions):
    """
    Computes the center of mass for each initial condition in the conditions list.
    Args:
        conditions (list): List of initial conditions, where each condition contains atom data.
    Returns:
        list: List of dictionaries with COM for each condition.
    """
    com_results = []
    for condition in conditions:
        total_mass = 0.0
        com_x, com_y, com_z = 0.0, 0.0, 0.0
        for atom in condition["atoms"]:
            mass = atom["mass"]
            total_mass += mass
            com_x += mass * atom["x"]
            com_y += mass * atom["y"]
            com_z += mass * atom["z"]
        
        # Normalize by total mass
        com = {
            "index": condition["index"],
            "com_x": com_x / total_mass,
            "com_y": com_y / total_mass,
            "com_z": com_z / total_mass
        }
        com_results.append(com)
    
    return com_results


def move_to_z_x(conditions, z_x_values):
    """
    Moves the center of mass of the molecule to the specified (x, 0, z) for each initial condition.
    Args:
        conditions (list): List of initial conditions, where each condition contains atom data.
        z_x_values (list): List of (z, x) tuples for each initial condition.
    """
    for condition, (z, x) in zip(conditions, z_x_values):
        current_com = compute_center_of_mass([condition])[0]
        shift_x = current_com["com_x"] - x
        shift_y = current_com["com_y"]
        shift_z = current_com["com_z"] - z

        for atom in condition["atoms"]:
            atom["x"] -= shift_x
            atom["y"] -= shift_y
            atom["z"] -= shift_z


def parse_initconds_file(filename):
    """
    Parses Initial Conditions file and extracts the relevant information.
    Args:
        filename (str): Path to the initconds file (e.g., file from --o1 or --o2).
    Returns:
        dict: A dictionary containing parsed information with keys for global data and each index.
    """
    data = {"header": {}, "equilibrium": [], "conditions": []}
    with open(filename, "r") as f:
        lines = f.readlines()

    # Parse the header
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("Ninit") or line.startswith("Natom") or line.startswith("Repr") or \
           line.startswith("Temp") or line.startswith("Eref") or line.startswith("Eharm"):
            key, value = line.split()
            data["header"][key] = float(value) if key != "Repr" else value
        elif line.startswith("Equilibrium"):
            i += 1
            while lines[i].strip():
                eq_line = lines[i].strip().split()
                data["equilibrium"].append({
                    "element": eq_line[0],
                    "atomic_number": float(eq_line[1]),
                    "x": float(eq_line[2]),
                    "y": float(eq_line[3]),
                    "z": float(eq_line[4]),
                    "mass": float(eq_line[5]),
                    "vx": float(eq_line[6]),
                    "vy": float(eq_line[7]),
                    "vz": float(eq_line[8]),
                })
                i += 1
        elif line.startswith("Index"):
            condition = {
                "index": int(line.split()[1]),
                "atoms": [],
                "states": {}
            }
            i += 1
            while i < len(lines) and lines[i].strip() != "States":
                atom_line = lines[i].strip().split()
                if len(atom_line) == 9:
                    condition["atoms"].append({
                        "element": atom_line[0],
                        "atomic_number": float(atom_line[1]),
                        "x": float(atom_line[2]),
                        "y": float(atom_line[3]),
                        "z": float(atom_line[4]),
                        "mass": float(atom_line[5]),
                        "vx": float(atom_line[6]),
                        "vy": float(atom_line[7]),
                        "vz": float(atom_line[8]),
                    })
                i += 1
            i += 1  # Skip "States"
            while i < len(lines) and lines[i].strip():
                state_line = lines[i].strip().split()
                condition["states"][state_line[0]] = {
                    "value_au": float(state_line[1]),
                    "value_ev": float(state_line[3])
                }
                i += 1
            data["conditions"].append(condition)
        i += 1


    return data


def read_initconds_files1(o):
    """
    Reads the initial condition files specified in --o1 and --o2.
    Args:
        o (str): File name for molecule 1 or 2.
    Returns:
        tuple: Parsed data for molecule 1 or 2..
    """
    if not os.path.isfile(o):
        raise FileNotFoundError(f"{o} not found.")

    print(f"Reading {o}...")
    molecule_data = parse_initconds_file(o)

    return molecule_data


def read_initconds_files(o1, o2):
    """
    Reads the initial condition files specified in --o1 and --o2.
    Args:
        o1 (str): File name for molecule 1.
        o2 (str): File name for molecule 2.
    Returns:
        tuple: Parsed data for molecule 1 and molecule 2.
    """
    if not os.path.isfile(o1):
        raise FileNotFoundError(f"{o1} not found.")
    if not os.path.isfile(o2):
        raise FileNotFoundError(f"{o2} not found.")

    print(f"Reading {o1}...")
    molecule1_data = parse_initconds_file(o1)
    print(f"Reading {o2}...")
    molecule2_data = parse_initconds_file(o2)

    return molecule1_data, molecule2_data


def delete_files(files):
    """
    Deletes the specified files.

    Args:
        files (list): List of file paths to delete.
    """
    for file in files:
        try:
            os.remove(file)
            print(f"Deleted file: {file}")
        except FileNotFoundError:
            print(f"File not found, could not delete: {file}")
        except Exception as e:
            print(f"Error deleting file {file}: {e}")


def check_state_selected_path():
    """
    Check if state_selected.py is in the same folder as bimolecular_collision.py.
    If not found, raise an error.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    state_selected_path = os.path.join(script_dir, "state_selected.py")
    if not os.path.isfile(state_selected_path):
        raise FileNotFoundError(f"Error: state_selected.py not found in the same folder as bimolecular_collision.py ({script_dir}).")
    return state_selected_path


def construct_args(options, molecule_suffix, input_file, mapping):
    """
    Construct the command-line arguments for state_selected.py for a specific molecule.
    - options: Parsed options from bimolecular_collision.py
    - molecule_suffix: '1' or '2', indicating the molecule
    - input_file: The input file (e.g., a.molden for molecule 1, b.molden for molecule 2)
    - mapping: Dictionary mapping molecule-specific options to state_specific.py generic options
    """
    args = ["python", os.path.join(os.path.dirname(__file__), "state_selected.py")]

    # Add the shared option
    if options.n is not None:
        args.append("-n")
        args.append(str(options.n))

    # Add molecule-specific options
    options_dict = vars(options)
    molecule_options = [key for key in options_dict.keys() if key.endswith(molecule_suffix)]
    for key in molecule_options:
        value = options_dict[key]
        base_key = key[:-1]  # Remove the suffix (e.g., "vibselect1" -> "vibselect")
        if base_key in mapping:
            generic_option = mapping[base_key]
            if isinstance(value, bool):
                if value:  # Add flags for boolean options
                    args.append(generic_option)
            elif value is not None:  # Add key-value pairs for other options
                if isinstance(value, list):  # Handle nargs=1 options
                    args.append(generic_option)
                    args.append(str(value[0]))
                else:
                    args.append(generic_option)
                    args.append(str(value))

    args.append(input_file)

    return args


def main():
  '''Main routine'''

  usage='''
bimolecular_collision.py [options] filename1.molden filename2.molden

Notice that if one of the systems is an atom, one only needs a single molden file.
For example, 

bimolecular_collision.py --system '1+2' filename2.molden
bimolecular_collision.py --system '2+1' filename1.molden

where "--system '1+2'" is for atom + diatom; "--system '2+1'" is for diatom + atom 

NOTICE:
    Our convention is using second molecule to collide towards the first molecule.

This script reads two MOLDEN files containing frequencies and normal modes [1]
and generates a state selected distribution of geometries and velocities for each molecule

Part of the code is adopted from wigner.py

Author: Yinan Shu
'''

  description=''

  parser = OptionParser(usage=usage, description=description)
  parser.add_option('-n', dest='n', type=int, nargs=1, default=3, help="Number of geometries to be generated (integer, default=3)")
  parser.add_option('--system', dest='system', type=str, nargs=1, default='3+3', help="type of the two systems, default is '3+3' which means its polyatomic molecules collides with polyatomic molecules, and '3' stands for polyatomic molecules with 3 or more atoms; options can be '1+2', '1+3', '2+1', '2+2', '2+3', '3+2', and '3+3', where 1 and 2 stand for atom and diatom respectively")
  parser.add_option('--atom', dest='atom', type=str, nargs=1, default='H', help="the atom involved in 1+2, 1+3, 2+1, and 3+1 systems")


  parser.add_option('--bmin', dest='bmin', type=float, nargs=1, default=0.0, help="minimal value of impact parameter in unit of angstrom")
  parser.add_option('--bmax', dest='bmax', type=float, nargs=1, default=5.0, help="maximum value of impact parameter in unit of angstrom")
  parser.add_option('--strata', dest='strata', type=int, nargs=1, default=1, help="number of strata between bmin and bmax, default value is 1, if a value larger than 1 is given, for example, --strata 3, then integer(n/3) initial conditions is given for each strata for the first 2 strata, and last strata has (n-2n/3) geometries. In addition, for ith strata, the bmin is bmin+(i-1)*(bmax-bmin), and bmax is bmin+i*(bmax-bmin)")
  parser.add_option('--separation', dest='separation', type=float, nargs=1, default=10.0, help="initial separation of the center of mass of two molecules, in unit of angstrom")
  parser.add_option('--relative_trans', dest='relative_trans',  type=float, nargs=1, default=1.0, help="initial relative translational energy between the center of mass of two molecules, in unit of eV")


  parser.add_option('-o', dest='o', type=str, nargs=1, default='initconds', help="Output filename (string, default=""initconds"")")
  parser.add_option('-x', dest='X', action='store_true',help="Generate a xyz file with the sampled geometries in addition to the initconds file")

  #options for molecule 1
  parser.add_option('--vibselect1', dest='vibselect1', type=int, nargs=1, default=6, help="Method of selection of vibrational mode energy (integer, default=1) for MOLECULE 1. 1 The user provides vibrational quantum numbers by the keyword vibstate1=(n1,n2,...,n3N-6) for a local minimum and vibstate1=(n1,n2,...,n3N-7) for a saddle point. 2 The program assigns vibrational quantum numbers at random, selected out of a Boltzmann distribution at a user-specified temperature that is specified by the keyword -t. 3 The program  generates an initial velocity from a Maxwell thermal distribution at a given temperature by -t keyword. This is an option for canonical ensemble, not an option for state selected ensemble. 4 The amount of energy in each mode is the same for each mode and is E1 by settting keyword --vibene1 E1. The unit for E1 is eV. 5 The amount of energy in mode m is Em, which can be different for each mode. Set --vibene1 E1, E2, ..., E3N-6 or --vibene1 E1, E2, ..., E3N-7. The units for Em are eV. 6 Like vibselect1=4 except that Em is calculated by the program as min[0.5hv, input E1]. 7 Like --vibselect1 5 except that Em is calculated by the program as  min[0.5hv, input Em]. 8 The user provides vibrational quantum numbers by keyword viblist1=(m1,n1;m2,n2;...,m3N-6,n3N-6), which only specifies the modes with non-zero vibrational quantum numbers, the vibrational quantum number of rest modes not specified in viblist1 are set to 0")
  parser.add_option('--vibdist1', dest='vibdist1', type=int, nargs=1, default=0, help="vibdist1 determines the type of phase space distribution for MOLECULE 1. 0 Default, classical or quasiclassical distribution. Uniform distribution. This distribution is quasiclassical if vibselect1 = 1 or 2, and it is classical if vibselect1>=4. 1 ground-state harmonic oscillator distribution. 2 wigner distribution.")
  parser.add_option('--vibstate1', dest='vibstate1', nargs=1, default="0", help="vibstate1 is a list of level of vibrational state for each mode, separated by comma, required by vibselect1=1. Example: --vibstate1 0,0,0,1,5")
  parser.add_option('--viblist1', dest='viblist1', nargs=1, default="0", help="viblist1 is a list of modes whose vibrational quantum numbers are non-zero, each pair (index of modes, vibrational quantum number, which are separated by comma) is separated question mark. Notice viblist1 is only used when set vibselect1 to 8, the modes that are not provided in viblist1 will have zero vibrational quantum number. Also notice that if a non-integer vibrational quantum number is provided, it will convert to the lower integer. Example: --viblist1 1,1?5,3")
  parser.add_option('--vibene1', dest='vibene1', nargs=1, default="0.0", help="vibene1 is a list of energies for each mode, separated by comma, required by vibselect1=4,5,6,7. Example: --vibene1 1.2,3.1,2.3")
  parser.add_option('--method1', dest='method1', type=int, nargs=1, default=0, help="method1 determins the level of energy approximation for MOLECULE 1. 0 use harmonic oscillator. 1 use directly computed potential energy (requires a lot calculations)")
  parser.add_option('--template1', dest='template1', type=str, nargs=1, default='MOLPRO.template1', help="Template filename")
  parser.add_option('--m1', dest='m1', action='store_true',help="Enter non-default atom masses for MOLECULE 1")
  parser.add_option('--s1', dest='s1', type=float, nargs=1, default=1.0, help="Scaling factor for the energies (float, default=1.0) for MOLECULE 1")
  parser.add_option('--t1', dest='t1', type=float, nargs=1, default=0., help="Temperature (float, default=0.0) for MOLECULE 1")
  parser.add_option('--T1', dest='T1', action='store_true', help="Discard high vibrational states in the temperature sampling for MOLECULE 1")
  parser.add_option('--L1', dest='L1', type=float, nargs=1, default=10.0, help="Discard frequencies below this value in cm-1 (float, default=10.) for MOLECULE 1")
  parser.add_option('--r1', dest='r1', type=int, nargs=1, default=16661, help="Seed for the random number generator (integer, default=16661) for MOLECULE 1")
  parser.add_option('--f1', dest='f1', type=int, nargs=1, default='0', help="Define the type of read normal modes for MOLECULE 1. 0 for automatic assignement, 1 for gaussian-type normal modes (Gaussian, Turbomole, Q-Chem, ADF, Orca), 2 for cartesian normal modes (Molcas, Molpro), 3 for Columbus-type (Columbus), or 4 for mass-weighted. (integer, default=0)")
  parser.add_option('--o1', dest='o1', type=str, nargs=1, default='initconds1', help="Output filename for MOLECULE 1 (string, default=""initconds"")")
  parser.add_option('--x1', dest='x1', action='store_true',help="Generate a xyz file with the sampled geometries in addition to the initconds file for MOLECULE 1")
  parser.add_option('--keep_trans_rot1', dest='KTR1', action='store_true',help="Keep translational and rotational components for MOLECULE 1")
  parser.add_option('--use_eq_geom1',    dest='UEG1', action='store_true',help="For all samples, use the equilibrium geometry (only sampled velocities are used, therefore, the mode energies are not correct) for MOLECULE 1")
  parser.add_option('--use_zero_veloc1', dest='UZV1', action='store_true',help="For all samples, set velocities to zero (only sampled geometries are used, therefore, the the mode energies are not correct) for MOLECULE 1")

  #options for molecule 2
  parser.add_option('--vibselect2', dest='vibselect2', type=int, nargs=1, default=6, help="Method of selection of vibrational mode energy (integer, default=1) for MOLECULE 2. 1 The user provides vibrational quantum numbers by the keyword vibstate2=(n1,n2,...,n3N-6) for a local minimum and vibstate2=(n1,n2,...,n3N-7) for a saddle point. 2 The program assigns vibrational quantum numbers at random, selected out of a Boltzmann distribution at a user-specified temperature that is specified by the keyword -t. 3 The program  generates an initial velocity from a Maxwell thermal distribution at a given temperature by -t keyword. This is an option for canonical ensemble, not an option for state selected ensemble. 4 The amount of energy in each mode is the same for each mode and is E1 by settting keyword --vibene2 E1. The unit for E1 is eV. 5 The amount of energy in mode m is Em, which can be different for each mode. Set --vibene2 E1, E2, ..., E3N-6 or --vibene2 E1, E2, ..., E3N-7. The units for Em are eV. 6 Like vibselect2=4 except that Em is calculated by the program as min[0.5hv, input E1]. 7 Like --vibselect2 5 except that Em is calculated by the program as  min[0.5hv, input Em]. 8 The user provides vibrational quantum numbers by keyword viblist2=(m1,n1;m2,n2;...,m3N-6,n3N-6), which only specifies the modes with non-zero vibrational quantum numbers, the vibrational quantum number of rest modes not specified in viblist2 are set to 0")
  parser.add_option('--vibdist2', dest='vibdist2', type=int, nargs=1, default=0, help="vibdist2 determines the type of phase space distribution for MOLECULE 2. 0 Default, classical or quasiclassical distribution. Uniform distribution. This distribution is quasiclassical if vibselect2 = 1 or 2, and it is classical if vibselect2>=4. 1 ground-state harmonic oscillator distribution. 2 wigner distribution.")
  parser.add_option('--vibstate2', dest='vibstate2', nargs=1, default="0", help="vibstate2 is a list of level of vibrational state for each mode, separated by comma, required by vibselect2=1. Example: --vibstate2 0,0,0,1,5")
  parser.add_option('--viblist2', dest='viblist2', nargs=1, default="0", help="viblist2 is a list of modes whose vibrational quantum numbers are non-zero, each pair (index of modes, vibrational quantum number, which are separated by comma) is separated question mark. Notice viblist2 is only used when set vibselect2 to 8, the modes that are not provided in viblist2 will have zero vibrational quantum number. Also notice that if a non-integer vibrational quantum number is provided, it will convert to the lower integer. Example: --viblist2 1,1?5,3")
  parser.add_option('--vibene2', dest='vibene2', nargs=1, default="0.0", help="vibene2 is a list of energies for each mode, separated by comma, required by vibselect2=4,5,6,7. Example: --vibene2 1.2,3.1,2.3")
  parser.add_option('--method2', dest='method2', type=int, nargs=1, default=0, help="method2 determins the level of energy approximation for MOLECULE 2. 0 use harmonic oscillator. 1 use directly computed potential energy (requires a lot calculations)")
  parser.add_option('--template2', dest='template2', type=str, nargs=1, default='MOLPRO.template2', help="Template filename")
  parser.add_option('--m2', dest='m2', action='store_true',help="Enter non-default atom masses for MOLECULE 2")
  parser.add_option('--s2', dest='s2', type=float, nargs=1, default=1.0, help="Scaling factor for the energies (float, default=1.0) for MOLECULE 2")
  parser.add_option('--t2', dest='t2', type=float, nargs=1, default=0., help="Temperature (float, default=0.0) for MOLECULE 2")
  parser.add_option('--T2', dest='T2', action='store_true', help="Discard high vibrational states in the temperature sampling for MOLECULE 2")
  parser.add_option('--L2', dest='L2', type=float, nargs=1, default=10.0, help="Discard frequencies below this value in cm-1 (float, default=10.) for MOLECULE 2")
  parser.add_option('--r2', dest='r2', type=int, nargs=1, default=16661, help="Seed for the random number generator (integer, default=16661) for MOLECULE 2")
  parser.add_option('--f2', dest='f2', type=int, nargs=1, default='0', help="Define the type of read normal modes for MOLECULE 2. 0 for automatic assignement, 1 for gaussian-type normal modes (Gaussian, Turbomole, Q-Chem, ADF, Orca), 2 for cartesian normal modes (Molcas, Molpro), 3 for Columbus-type (Columbus), or 4 for mass-weighted. (integer, default=0)")
  parser.add_option('--o2', dest='o2', type=str, nargs=1, default='initconds2', help="Output filename for MOLECULE 1 (string, default=""initconds"")")
  parser.add_option('--x2', dest='x2', action='store_true',help="Generate a xyz file with the sampled geometries in addition to the initconds file for MOLECULE 2")
  parser.add_option('--keep_trans_rot2', dest='KTR2', action='store_true',help="Keep translational and rotational components for MOLECULE 2")
  parser.add_option('--use_eq_geom2',    dest='UEG2', action='store_true',help="For all samples, use the equilibrium geometry (only sampled velocities are used, therefore, the mode energies are not correct) for MOLECULE 2")
  parser.add_option('--use_zero_veloc2', dest='UZV2', action='store_true',help="For all samples, set velocities to zero (only sampled geometries are used, therefore, the the mode energies are not correct) for MOLECULE 2")

  (options, args) = parser.parse_args()

  mapping = {
    "vibselect": "--vibselect",
    "vibdist": "--vibdist",
    "vibstate": "--vibstate",
    "viblist": "--viblist",
    "vibene": "--vibene",
    "method": "--method",
    "template": "--template",
    "m": "-m",
    "s": "-s",
    "t": "-t",
    "T": "-T",
    "L": "-L",
    "r": "-r",
    "f": "-f",
    "o": "-o",
    "x": "-x"
  }

  # check if state_selected.py exists
  state_selected_path = check_state_selected_path()

  system=options.system

  #==============
  # Sampling coordinates
  #==============
  # used for sampling coordinates 
  n=options.n
  bmin=options.bmin*ANG_TO_BOHR
  bmax=options.bmax*ANG_TO_BOHR
  strata=options.strata
  separation=options.separation*ANG_TO_BOHR
  #compute initial separation 
  z_x_values=sample_z_x(n, bmin, bmax, strata, separation)

  # first check how many inputs are given 
  if system=='2+3' or system=='3+2' or system=='3+3': 
      print("====================================================")
      print("Initial Condition Sampling for MOLECULE + MOLECULE")
      print("====================================================")
      if len(args) != 2:
          parser.error("You must provide exactly two input files (two molden files for two molecules)")
      input_file1, input_file2 = args
      print("perform initial condition sampling for MOLECULE 1 using state_selected.py")
      args1 = construct_args(options, "1", input_file1, mapping)
      print("commend:", args1)
      subprocess.run(args1)
      print("perform initial condition sampling for MOLECULE 2 using state_selected.py")
      args2 = construct_args(options, "2", input_file2, mapping)
      print("commend:", args2)
      subprocess.run(args2)
      o1=options.o1
      o2=options.o2
      molecule1_data, molecule2_data = read_initconds_files(o1, o2)
      #move molecule 1 to origin 
      print("center of mass of MOLECULE 1 is moved to (0, 0, 0)")
      move_to_origin(molecule1_data["conditions"])
      #move molecule 2 to (x,y,0)
      print("center of mass of MOLECULE 2 is moved to (impact_parameter, 0, initial_separation)")
      move_to_z_x(molecule2_data["conditions"], z_x_values)
      if not DEBUG:
          delete_files([o1, o2])
      if options.x1 and not DEBUG:
          delete_files([options.o1+'.xyz'])
      if options.x2 and not DEBUG:
          delete_files([options.o2+'.xyz'])
  elif system=='1+2' or system=='1+3': 
      print("====================================================")
      print("Initial Condition Sampling for ATOM + MOLECULE")
      print("====================================================")
      atom=options.atom
      print("Involved atom is:", atom)
      if len(args) !=1:
          arser.error("You must provide exactly one input file (one molden file for the second molecule)")
      input_file2 = args[0]
      # generate atom data (molecule 1)
      print("place atom at (0, 0, 0)")
      molecule1_data = generate_atom_data(n, atom, 1, z_x_values)
      print("perform initial condition sampling for MOLECULE 2 using state_selected.py")
      args2 = construct_args(options, "2", input_file2, mapping)
      print("commend:", args2)
      subprocess.run(args2)
      o2=options.o2    
      molecule2_data = read_initconds_files1(o2)
      # move molecule 2 to (x,y,0)
      print("center of mass of MOLECULE 2 is moved to (impact_parameter, 0, initial_separation)")
      move_to_z_x(molecule2_data["conditions"], z_x_values)
      if not DEBUG:
          delete_files([o2])
      if options.x2 and not DEBUG:
          delete_files([options.o2+'.xyz'])

  elif system=='2+1' or system=='3+1':
      print("====================================================")
      print("Initial Condition Sampling for MOLECULE + ATOM")
      print("====================================================")
      atom=options.atom
      print("Involved atom is:", atom)
      if len(args) !=1:
          arser.error("You must provide exactly one input file (one molden file for the first molecule)") 
      input_file1 = args[0]
      print("perform initial condition sampling for MOLECULE 1 using state_selected.py")
      args1 = construct_args(options, "1", input_file1, mapping)
      print("commend:", args1)
      subprocess.run(args1)
      o1=options.o1
      molecule1_data = read_initconds_files1(o1)
      # move molecule 1 to origin 
      print("center of mass of MOLECULE 1 is moved to (0, 0, 0)")
      move_to_origin(molecule1_data["conditions"])
      # generate atom data (molecule 2)
      print("place atom at (impact_parameter, 0, initial_separation)")
      molecule2_data = generate_atom_data(n, atom, 2, z_x_values)
      if not DEBUG:
          delete_files([o1])
      if options.x1 and not DEBUG:
          delete_files([options.o1+'.xyz'])

  if not DEBUG:
      delete_files(['KEYSTROKES.state_selected'])

  #==============
  # Sampling collision velocity
  #==============
  # add relative translation energy to molecule2
  E_col=options.relative_trans/HARTREE_TO_EV
  add_collision_velocity(molecule2_data, E_col)


  #==============
  # Put to initial condition together
  #==============
  combined_data=combine_molecule_data(molecule1_data, molecule2_data, E_col)
  outfile=options.o
  write_combined_data(combined_data, outfile)

  if options.X:
      write_xyz_coordinates(combined_data,options.o+'.xyz')

  # save the shell command
  command='python '+' '.join(sys.argv)
  f=open('KEYSTROKES.bimolecular_collision','w')
  f.write(command)
  f.close()

# ======================================================================================================================

if __name__ == '__main__':
    main()

