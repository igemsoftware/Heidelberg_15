#!/usr/bin/env python

from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
import mpmath as math
import numpy as np
import pexpect
import os
from  mpmath import mp as math
import sys
import xml.etree.ElementTree as etree
from sys import stdout
#from jug import *
#from jug import mapreduce as mr
#from jug.task import *
#from jug.compound import *
import shutil
import os.path as path
import time
from multiprocessing import *
import argparse

forcefield_name = "leaprc.ff12SB"

def gen_range(lower, upper, step):
    start = lower
    while start < upper:
        yield start
        start += step

def sum_a(nested_array):
    output = list(nested_array[0])
    n_a = list(nested_array)
    for i in range(1, len(n_a)):
        for j in range(len(n_a[0])):
            output[j] += n_a[i][j]
    return output
        
def arange(lower, upper, step):
    array = [elem for elem in gen_range(lower, upper, step)]
    return array

def times(number, array):
    from  mpmath import mp as math
    math.prec = 200
    output = [math.fmul(math.mpmathify(elem),math.mpmathify(number)) for elem in array]
    return output

def times_a(array1, array2):
    from  mpmath import mp as math
    math.prec = 200
    a1 = list(array1)
    a2 = list(array2)
    output = [math.fmul(math.mpmathify(a1[i]),math.mpmathify(a2[i])) for i in range(len(a1))]
    return output

def power(array,a):
    output = [elem**a for elem in array]
    return output

def power2(array_nested, a):
    output = [power(elem,a) for elem in array_nested]
    return output

def exp(array):
    from  mpmath import mp as math
    output = [math.exp(math.mpmathify(elem)) for elem in array]
    return output

def log(array):
    from  mpmath import mp as math
    output = [math.log(math.mpmathify(elem)) for elem in array]
    return output

def absolute(vector):
    return math.sqrt(sum([elem**2 for elem in vector]))

def vector_angle(vec1, vec2):
    abs_vec1 = absolute(vec1)
    abs_vec2 = absolute(vec2)
    prod_12 = sum([alem*blem for alem,blem in zip(vec1, vec2)])
    return math.acos(prod_12/(abs_vec1*abs_vec2))


def cross(vec1, vec2):
    return [vec1[1]*vec2[2]-vec1[2]*vec2[1],vec1[2]*vec2[0]-vec1[0]*vec2[2],vec1[0]*vec2[1]-vec1[1]*vec2[0]]

def constrainPO3(topology,system):
    C3s = []
    O3s = []
    Ps = []
    C5s = []
    for elem in topology.atoms():
        if elem.name == "C3'":
            C3s.append(elem.index)
        elif elem.name == "O3'":
            O3s.append(elem.index)
        elif elem.name == "P":
            Ps.append(elem.index)
        elif elem.name == "C5'":
            C5s.append(elem.index)
        else:
            pass
    force = mm.CustomTorsionForce('-k0*(theta^2)')
    force.setForceGroup(2)
    force.addPerTorsionParameter('k0')

    for i in range(len(Ps)):
        force.addTorsion(C3s[i],O3s[i],Ps[i],C5s[i+1],[20])
    system.addForce(force)
    print(C3s, O3s, Ps, C5s)
    
def constrainPO5(topology,system):
    C5s = []
    O5s = []
    Ps = []
    O3s = []
    for elem in topology.atoms():
        if elem.name == "O3'":
            O3s.append(elem.index)
        elif elem.name == "O5'":
            O5s.append(elem.index)
        elif elem.name == "P":
            Ps.append(elem.index)
        elif elem.name == "C5'":
            C5s.append(elem.index)
        else:
            pass
    force = mm.CustomTorsionForce('-k0*(theta^2)')
    force.setForceGroup(2)
    force.addPerTorsionParameter('k0')
    for i in range(len(Ps)):
        force.addTorsion(C5s[i+1],O5s[i+1],Ps[i],O3s[i],[20])
    system.addForce(force)
    print(C5s, O5s, Ps, O3s)
    
def constrainC5O5(topology,system):
    C5s = []
    O5s = []
    Ps = []
    H5s = []
    for elem in topology.atoms():
        if elem.name == "C5'":
            C5s.append(elem.index)
        elif elem.name == "O5'":
            O5s.append(elem.index)
        elif elem.name == "P":
            Ps.append(elem.index)
        elif elem.name == "C4'":
            H5s.append(elem.index)
        else:
            pass
    force = mm.CustomTorsionForce('-k0*(theta^2)')
    force.setForceGroup(2)
    force.addPerTorsionParameter('k0')
    for i in range(len(Ps)):
        force.addTorsion(H5s[i+1],C5s[i+1],O5s[i+1],Ps[i],[20])
    system.addForce(force)
    print(C5s, O5s, Ps, H5s)

def get_aptamer(ligand_range, positions):
    return positions[ligand_range[1]-1:]

def cross(vec1, vec2):
    return [vec1[1]*vec2[2]-vec1[2]*vec2[1],vec1[2]*vec2[0]-vec1[0]*vec2[2],vec1[0]*vec2[1]-vec1[1]*vec2[0]]

def get_ligand(topology):
    ligand_indices = []
    for a in topology.atoms():
        if a.residue.name not in ["DGN","DAN","DTN","DCN","DG","DA","DT","DC","DG5","DA5","DT5","DC5","DG3","DA3","DT3","DC3"]:
            ligand_indices.append(a.index)
    return ligand_indices

def get_ligand_range(topology):
    return [get_ligand(topology)[0],len(get_ligand(topology))]

def get_offset(positions_old, positions):
    vec_a = (positions[len(positions_old)-1]-positions[len(positions_old)-2])
    vec_b = (positions_old[len(positions_old)-1]-positions_old[-2])
    alpha = math.acos(sum([alem.value_in_unit(unit.angstroms)*blem.value_in_unit(unit.angstroms) for alem, blem in zip(vec_a,vec_b)])/(np.linalg.norm(vec_a.value_in_unit(unit.angstroms))*np.linalg.norm(vec_b.value_in_unit(unit.angstroms))))
    alpha_t = 0.
    #alpha_t = math.pi-113.3*math.pi/360.
    d_alpha = alpha_t-alpha
    axis = cross(vec_a.value_in_unit(unit.angstroms),vec_b.value_in_unit(unit.angstroms))/np.linalg.norm(cross(vec_a.value_in_unit(unit.angstroms),vec_b.value_in_unit(unit.angstroms)))
    offset = positions_old[-1]-positions[len(positions_old)-1]
    return d_alpha, axis, offset, vec_a, vec_b
    
def position_aptamer(positions_old, positions):
    alpha, axis, offset, vec_a, vec_b = get_offset(positions_old, positions)
    ps = positions
    #print(alpha.__class__.__name__)
    phi_2 = (alpha/2).real
    #print(alpha/2)
    x, y, z = axis
    s = np.math.sin(phi_2)
    c = np.math.cos(phi_2)
    rot = np.array([[2*(np.power(x,2)-1)*np.power(s,2)+1, 2*x*y*np.power(s,2)-2*z*c*s, 2*x*z*np.power(s,2)+2*y*c*s],
                    [2*x*y*np.power(s,2)+2*z*c*s, 2*(np.power(y,2)-1)*np.power(s,2)+1, 2*z*y*np.power(s,2)-2*x*c*s],
                    [2*x*z*np.power(s,2)-2*y*c*s, 2*z*y*np.power(s,2)+2*x*c*s, 2*(np.power(z,2)-1)*np.power(s,2)+1]])

    for j in range(len(positions_old)-2,len(positions)):
        roted = np.dot(np.array(positions[j].value_in_unit(unit.angstrom)),rot)
        ps[j] = mm.Vec3(roted[0],roted[1],roted[2])*unit.angstrom
    drift = positions_old[-1] - ps[len(positions_old)-1]
    for j in range(len(positions_old)-2,len(positions)):
        ps[j] += drift+vec_b.value_in_unit(unit.angstroms)/np.linalg.norm(vec_b.value_in_unit(unit.angstroms))*.6*unit.angstroms
    positions_new = positions_old[:-1]+ps[len(positions_old)-1:]
    return positions_new

def get_offset_five(positions_old, positions, ligand_length):
    vec_a = (positions[ligand_length+(len(positions)-len(positions_old))+1]-positions[ligand_length+(len(positions)-len(positions_old))])
    vec_b = (positions_old[ligand_length+1]-positions_old[ligand_length])
    alpha = math.acos(sum([alem.value_in_unit(unit.angstroms)*blem.value_in_unit(unit.angstroms) for alem, blem in zip(vec_a,vec_b)])/(np.linalg.norm(vec_a.value_in_unit(unit.angstroms))*np.linalg.norm(vec_b.value_in_unit(unit.angstroms))))
    alpha_t = 0.
    #alpha_t = math.pi-113.3*math.pi/360.
    d_alpha = alpha_t-alpha
    axis = cross(vec_a.value_in_unit(unit.angstroms),vec_b.value_in_unit(unit.angstroms))/np.linalg.norm(cross(vec_a.value_in_unit(unit.angstroms),vec_b.value_in_unit(unit.angstroms)))
    offset = positions_old[-1]-positions[len(positions_old)-1]
    return d_alpha, axis, offset, vec_a, vec_b

def position_aptamer_five(positions_old, positions, ligand_length):
    alpha, axis, offset, vec_a, vec_b = get_offset_five(positions_old, positions)
    ps = positions
    phi_2 = alpha/2
    x, y, z = axis
    s = np.math.sin(phi_2)
    c = np.math.cos(phi_2)
    rot = np.array([[2*(np.power(x,2)-1)*np.power(s,2)+1, 2*x*y*np.power(s,2)-2*z*c*s, 2*x*z*np.power(s,2)+2*y*c*s],
                    [2*x*y*np.power(s,2)+2*z*c*s, 2*(np.power(y,2)-1)*np.power(s,2)+1, 2*z*y*np.power(s,2)-2*x*c*s],
                    [2*x*z*np.power(s,2)-2*y*c*s, 2*z*y*np.power(s,2)+2*x*c*s, 2*(np.power(z,2)-1)*np.power(s,2)+1]])

    for j in range(ligand_length,ligand_length+(len(positions)-len(positions_old))):
        roted = np.dot(np.array(positions[j].value_in_unit(unit.angstrom)),rot)
        ps[j] = mm.Vec3(roted[0],roted[1],roted[2])*unit.angstrom
    drift = positions_old[-1] - ps[len(positions_old)-1]
    for j in range(len(positions_old)-2,len(positions)):
        ps[j] += drift+vec_b.value_in_unit(unit.angstroms)/np.linalg.norm(vec_b.value_in_unit(unit.angstroms))*.6*unit.angstroms
    positions_new = positions_old[:-1]+ps[len(positions_old)-1:]
    return positions_new
    

def fix_ligand(ligand):
    for prt in ligand:
        system.setParticleMass(prt,1e50)
        
def stratify(pos,energies,size,phi_size):
    energies = [math.mpmathify(elem) for elem in energies]
    countx1 = 0
    countx2 = 0
    meanx1 = 0
    meanx2 = 0
    county1 = 0
    county2 = 0
    meany1 = 0
    meany2 = 0
    countz1 = 0
    countz2 = 0
    meanz1 = 0
    meanz2 = 0
    counti1 = 0
    counti2 = 0
    meani1 = 0
    meani2 = 0
    countj1 = 0
    countj2 = 0
    meanj1 = 0
    meanj2 = 0
    countk1 = 0
    countk2 = 0
    meank1 = 0
    meank2 = 0
    countphi1 = 0
    countphi2 = 0
    meanphi1 = 0
    meanphi2 = 0
    for l in range(len(energies)):
        if pos[l][0] < 0.:
            countx1 += 1
            meanx1 += energies[l]
        else:
            countx2 += 1
            meanx2 += energies[l]
        if pos[l][1] < 0:
            county1 += 1
            meany1 += energies[l]
        else:
            county2 += 1
            meany2 += energies[l]
        if pos[l][2] < 0:
            countz1 += 1
            meanz1 += energies[l]
        else:
            countz2 += 1
            meanz2 = energies[l]
        if pos[l][3] < 0:
            counti1 += 1
            meani1 += energies[l]
        else:
            counti2 += 1
            meani2 += energies[l]
        if pos[l][4] < 0:
            countj1 += 1
            meanj1 += energies[l]
        else:
            countj2 += 1
            meanj2 += energies[l]
        if pos[l][5] < 0:
            countk1 += 1
            meank1 += energies[l]
        else:
            countk2 += 1
            meank2 += energies[l]
        if pos[l][6] < 0:
            countphi1 += 1
            meanphi1 += energies[l]
        else:
            countphi2 += 1
            meanphi2 += energies[l]
    meanx1 /= countx1+1e-20
    print(meanx1)
    meanx2 /= countx2+1e-20
    meany1 /= county1+1e-20
    meany2 /= county2+1e-20
    meanz1 /= countz1+1e-20
    meanz2 /= countz2+1e-20
    meani1 /= counti1+1e-20
    meani2 /= counti2+1e-20
    meanj1 /= countj1+1e-20
    meanj2 /= countj2+1e-20
    meank1 /= countk1+1e-20
    meank2 /= countk2+1e-20
    meanphi1 /= countphi1+1e-20
    meanphi2 /= countphi2+1e-20
    
    varx1 = np.math.sqrt(np.array([np.power(energies[l]-meanx1,2) for l in range(len(energies)) if pos[l][0] < 0.]).sum()/(countx1-1+1e-20))
    print(varx1)
    vary1 = np.math.sqrt(np.array([np.power(energies[l]-meany1,2) for l in range(len(energies)) if pos[l][1] < 0.]).sum()/(county1-1+1e-20))
    varz1 = np.math.sqrt(np.array([np.power(energies[l]-meanz1,2) for l in range(len(energies)) if pos[l][2] < 0.]).sum()/(countz1-1+1e-20))
    vari1 = np.math.sqrt(np.array([np.power(energies[l]-meani1,2) for l in range(len(energies)) if pos[l][3] < 0.]).sum()/(counti1-1+1e-20))
    varj1 = np.math.sqrt(np.array([np.power(energies[l]-meanj1,2) for l in range(len(energies)) if pos[l][4] < 0.]).sum()/(countj1-1+1e-20))
    vark1 = np.math.sqrt(np.array([np.power(energies[l]-meank1,2) for l in range(len(energies)) if pos[l][5] < 0.]).sum()/(countk1-1+1e-20))
    varphi1 = np.math.sqrt(np.array([np.power(energies[l]-meanphi1,2) for l in range(len(energies)) if pos[l][6] < 0.]).sum()/(countphi1-1+1e-20))
    varx2 = np.math.sqrt(np.array([np.power(energies[l]-meanx2,2) for l in range(len(energies)) if pos[l][0] >= 0.]).sum()/(countx2-1+1e-20))
    vary2 = np.math.sqrt(np.array([np.power(energies[l]-meany2,2) for l in range(len(energies)) if pos[l][1] >= 0.]).sum()/(county2-1+1e-20))
    varz2 = np.math.sqrt(np.array([np.power(energies[l]-meanz2,2) for l in range(len(energies)) if pos[l][2] >= 0.]).sum()/(countz2-1+1e-20))
    vari2 = np.math.sqrt(np.array([np.power(energies[l]-meani2,2) for l in range(len(energies)) if pos[l][3] >= 0.]).sum()/(counti2-1+1e-20))
    varj2 = np.math.sqrt(np.array([np.power(energies[l]-meanj2,2) for l in range(len(energies)) if pos[l][4] >= 0.]).sum()/(countj2-1+1e-20))
    vark2 = np.math.sqrt(np.array([np.power(energies[l]-meank2,2) for l in range(len(energies)) if pos[l][5] >= 0.]).sum()/(countk2-1+1e-20))
    varphi2 = np.math.sqrt(np.array([np.power(energies[l]-meanphi2,2) for l in range(len(energies)) if pos[l][6] >= 0.]).sum()/(countphi2-1+1e-20))
    
    x_var = [varx1,varx2]
    y_var = [vary1,vary2]
    z_var = [varz1,varz2]
    i_var = [vari1,vari2]
    j_var = [varj1,varj2]
    k_var = [vark1,vark2]
    phi_var = [varphi1,varphi2]
            
    return [x_var, y_var, z_var, i_var, j_var, k_var, phi_var]
            
def var_to_ratio(var):
    var = np.array(var)
    res = var/var.sum()
    return res

def uniform_strat(x_var,y_var,z_var,i_var,j_var,k_var,phi_var,size,phi_size):
    x = var_to_ratio(x_var)
    y = var_to_ratio(y_var)
    z = var_to_ratio(z_var)
    i = var_to_ratio(i_var)
    j = var_to_ratio(j_var)
    k = var_to_ratio(k_var)
    phi = var_to_ratio(phi_var)
    size = np.array(size)
    res = []
    oldlem = size
    
    stepx = (size[1]-size[0])/len(x)
    stepy = (size[1]-size[0])/len(y)
    stepz = (size[1]-size[0])/len(z)
    stepi = (phi_size[1]-phi_size[0])/len(i)
    stepj = (phi_size[1]-phi_size[0])/len(j)
    stepk = (phi_size[1]-phi_size[0])/len(k)
    stepphi = (phi_size[1]-phi_size[0])/len(phi)
    
    distx = np.random.choice(np.array([np.random.uniform(size[0]+(l-1)*stepx,size[0]+l*stepx) for l in range(1,len(x)+1)]),p=np.append(x[:-1],[1-sum(x[:-1])]))
    disty = np.random.choice(np.array([np.random.uniform(size[0]+(l-1)*stepy,size[0]+l*stepy) for l in range(1,len(y)+1)]),p=np.append(y[:-1],[1-sum(y[:-1])]))
    distz = np.random.choice(np.array([np.random.uniform(size[0]+(l-1)*stepz,size[0]+l*stepz) for l in range(1,len(z)+1)]),p=np.append(z[:-1],[1-sum(z[:-1])]))
    disti = np.random.choice(np.array([np.random.uniform(phi_size[0]+(l-1)*stepi,phi_size[0]+l*stepi) for l in range(1,len(i)+1)]),p=np.append(i[:-1],[1-sum(i[:-1])]))
    distj = np.random.choice(np.array([np.random.uniform(phi_size[0]+(l-1)*stepj,phi_size[0]+l*stepj) for l in range(1,len(j)+1)]),p=np.append(j[:-1],[1-sum(j[:-1])]))
    distk = np.random.choice(np.array([np.random.uniform(phi_size[0]+(l-1)*stepk,phi_size[0]+l*stepk) for l in range(1,len(k)+1)]),p=np.append(k[:-1],[1-sum(k[:-1])]))
    distphi = np.random.choice(np.array([np.random.uniform(phi_size[0]+(l-1)*stepphi,phi_size[0]+l*stepphi) for l in range(1,len(phi)+1)]),p=np.append(phi[:-1],[1-sum(phi[:-1])]))
    
    return distx, disty, distz, disti, distj, distk, distphi

def energy_strata(xyz,energies):
    energies = [math.mpmathify(elem) for elem in energies]
    strata = [[energies[l] for l in range(len(energies)) if (c*xyz[l][0] <= 0) and (b*xyz[l][1] <= 0) and (a*xyz[l][2] <= 0) and (d*xyz[l][3] <= 0) and (e*xyz[l][4] <= 0) and (f*xyz[l][5] <= 0) and (g*xyz[l][6] <= 0)] for a in [-1,1] for b in [-1,1] for c in [-1,1] for d in [-1,1] for e in [-1,1] for f in [-1,1] for g in [-1,1]]
    return strata

class Aptamer:

    global _FORMAT
    global _HYBRID

    def __init__(self, forcefield_name, ligand_mol2_path):
        self.process = pexpect.spawn('tleap -f'+forcefield_name)
        self.process.sendline('source leaprc.gaff')
        self.process.sendline("set default PBradii mbondi2")
        self.process.sendline("ligand = load"+_FORMAT+" "+ligand_mol2_path)
        if _FORMAT != "":
		self.process.sendline("loadamberparams "+_HYBRID)
	self.geometry = []
        self.position = [0, 0, 0]
        self.orientation = [0, 0, 0]
        #global lig_energy
        #self.lig_energy = lig_energy
        
    def atom_position(self, identifier, residue_ID, atom_ID):
        self.process.sendline("desc "+identifier+"."+str(residue_ID)+"."+str(atom_ID))
        self.process.expect("Atom position:")
        self.process.expect("Atom velocity")
        vector = eval("["+self.process.before.strip()+"]")
        return vector
    
    def bond_COM(self, identifier, residue_ID, atom1_ID, atom2_ID):
        vec1 = self.atom_position(identifier, residue_ID, atom1_ID)
        vec2 = self.atom_position(identifier, residue_ID, atom2_ID)
        return [(alem+blem)/2 for alem, blem in zip(vec1,vec2)]
        
    
    def command(self, command_text):
        #self.process.expect(">")
        self.process.sendline(command_text)
    
    def sequence(self, identifier, string_of_residues):
        inputstring = "{"+string_of_residues+"}"
        self.command(identifier+" = sequence "+inputstring)
        #self.command("union = combine { "+identifier+" ligand }")
        
    def unify(self, identifier):
        self.command("union = combine { ligand "+identifier+" }")
        
        
    def seq_first(self, identifier, string_of_residues):
        inputstring = "{"+string_of_residues+"}"
        self.command(identifier+" = sequence "+inputstring)
        
    def ligand(self, identifier, path_to_pdb):
        self.command(identifier+" = load"+_FORMAT+" "+path_to_pdb)
    
    def save_all(self, identifier, identifier1, identifier2, path):
        self.command(identifier+" = combine {"+identifier1+" "+identifier2+"}")
        self.command("saveamberparm "+identifier+" "+path+identifier+".prmtop "+path+identifier+".inpcrd")
    
    def torsionCN(self, identifier, residue, angle):
        self.command("impose "+identifier+" {"+str(residue)+"""} {{ "C2'" "C1'" "N9" "C4" """+str(angle)+"""}}""")
        self.command("impose "+identifier+" {"+str(residue)+"""} {{ "C2'" "C1'" "N1" "C6" """+str(angle)+"""}}""")
        self.geometry.append([residue, 1, angle])
        
    def torsionCN_first(self, identifier, residue, angle):
        self.command("impose "+identifier+" {"+str(residue)+"""} {{ "C2'" "C1'" "N9" "C4" """+str(angle)+"""}}""")
        self.command("impose "+identifier+" {"+str(residue)+"""} {{ "C2'" "C1'" "N1" "C6" """+str(angle)+"""}}""")
        self.geometry.append([residue, 1, angle])
        self.command("union = combine { "+identifier+" ligand }")
        
    def torsionPO5(self, identifier, residue, angle):
        self.command("impose "+identifier+" {"+str(residue)+"""} {{ "OP1" "P" "O5'" "C5'" """+str(angle)+"""}}""")
        self.geometry.append([residue, 2, angle])

    def torsionC5O5(self, identifier, residue, angle):
        self.command("impose "+identifier+" {"+str(residue)+"""} {{ "H5'" "C5'" "O5'" "P" """+str(angle)+"""}}""")
        self.geometry.append([residue, 3, angle])
        
    def torsionPO3(self, identifier, residue, angle):
        self.command("impose "+identifier+" {"+str(residue)+"""} {{ "C3'" "O3'" "P" "O5'" """+str(angle)+"""}}""")
        self.geometry.append([residue, 4, angle])
        
    def build_geometry(self, identifier):
        for elem in self.geometry:
            (residue, bond, angle) = elem
            bond_dict = {1 : torsionCN, 2 : torsionPO5, 3 : torsionC5O5, 4 : torsionPO3}
            apply(bond_dict[bond],[identifier, residue, angle])
        
    def translate(self, identifier, vector, step=1):
        self.command("translate "+identifier+"."+str(step+1)+" {%s %s %s}"%tuple(vector))
        #self.command("translate "+identifier+"."+str(step+1)+" {-%s -%s -%s}"%tuple(vector))
        for i in range(0,3):
            self.position[i] += vector[i]
            
    def rotate(self, identifier, angles, step=1):
        ang_num = 2*math.pi/360
        rotateZ = [math.cos(angles[0]*ang_num),-math.sin(angles[0]*ang_num),0,math.sin(angles[0]*ang_num),math.cos(angles[0]*ang_num),0,0,0,1]
        rotateY = [math.cos(angles[1]*ang_num),0,-math.sin(angles[1]*ang_num),0,1,0,math.sin(angles[1]*ang_num),0,math.cos(angles[1]*ang_num)]
        rotateX = [1,0,0,0,math.cos(angles[2]*ang_num),-math.sin(angles[2]*ang_num),0,math.sin(angles[2]*ang_num),math.cos(angles[2]*ang_num)]
        antirotateZ = [math.cos(-angles[0]*ang_num),-math.sin(-angles[0]*ang_num),0,math.sin(-angles[0]*ang_num),math.cos(-angles[0]*ang_num),0,0,0,1]
        antirotateY = [math.cos(-angles[1]*ang_num),0,-math.sin(-angles[1]*ang_num),0,1,0,math.sin(-angles[1]*ang_num),0,math.cos(-angles[1]*ang_num)]
        antirotateX = [1,0,0,0,math.cos(-angles[2]*ang_num),-math.sin(-angles[2]*ang_num),0,math.sin(-angles[2]*ang_num),math.cos(-angles[2]*ang_num)]
        #self.command("transform "+identifier+""+"{ {%s %s %s} {%s %s %s} {%s %s %s} }"%tuple(rotateZ))
        self.command("transform "+identifier+"."+str(step+1)+"{ {%s %s %s} {%s %s %s} {%s %s %s} }"%tuple(antirotateZ))
        #self.command("transform "+identifier+""+"{ {%s %s %s} {%s %s %s} {%s %s %s} }"%tuple(rotateY))
        self.command("transform "+identifier+"."+str(step+1)+"{ {%s %s %s} {%s %s %s} {%s %s %s} }"%tuple(antirotateY))
        #self.command("transform "+identifier+""+"{ {%s %s %s} {%s %s %s} {%s %s %s} }"%tuple(rotateX))
        self.command("transform "+identifier+"."+str(step+1)+"{ {%s %s %s} {%s %s %s} {%s %s %s} }"%tuple(antirotateX))
        for i in range(0,3):
            self.orientation[i] += angles[i]
            
    def get_offset(self,identifier,residue_idx_1,residue_idx_2):
        self.process.sendline("desc "+identifier+"."+str(residue_idx_1)+".1")
        self.process.expect("Atom position:")
	print self.process.before
	print self.process.after
        self.process.expect("Atom velocity")
        vector1 = eval("["+self.process.before.strip()+"]")
        #print vector1
	time.sleep(1)
        self.process.sendline("desc "+identifier+"."+str(residue_idx_2)+".1")
        print self.process.before
	print self.process.after
	self.process.expect("Atom position:")
	print self.process.before
	print self.process.after
        self.process.expect("Atom velocity")
        vector2 = eval("["+self.process.before.strip()+"]")
        #print vector2
        vector = [vector2[i] - vector1[i] for i in range(len(vector2))]
        return vector
    
    def reshift(self,identifier,to_residue,step):
        vector = self.get_offset(identifier,1,to_residue)
        self.translate(identifier,vector,step)
        
    def respin(self,identifier,step):
        self.rotate(identifier,[0,0,180],step=step)
        
    def rotate_axis(self, identifier, angle, axis, step):
        rotateM = rotate_axis(angle, axis)
        #self.command("transform "+identifier+""+"{ {%s %s %s} {%s %s %s} {%s %s %s} }"%tuple(rotateZ))
        self.command("transform "+identifier+"."+str(step+1)+"{ {%s %s %s} {%s %s %s} {%s %s %s} }"%tuple(rotateM))
            
    def relative_rescale(self, identifier, to_residue, step):
        self.process.sendline("desc union.1.1")
        self.process.expect("Atom position:")
        self.process.expect("Atom velocity")
        vector1 = eval("["+self.process.before.strip()+"]")
        if (step-2)%4 != 0:
            for i in range(step):
                self.respin("union",i-1)
        self.process.sendline("desc union.1.1")
        self.process.expect("Atom position:")
        self.process.expect("Atom velocity")
        vector2 = eval("["+self.process.before.strip()+"]")
        vector = [vector2[i]-vector1[i] for i in range(len(vector2))]
        self.translate("union", vector, step=step)
        self.reshift("union", to_residue, step)


def ligand_box(padding):
        from  mpmath import mp as math
        math.prec = 200
        global forcefield_name
        process = pexpect.spawn('tleap -f'+forcefield_name)
        process.sendline('source leaprc.gaff')
        process.sendline("set default PBradii mbondi2")
        process.sendline("ligand = load"+_FORMAT+" "+_INFILE)
        process.sendline("saveamberparm ligand ligand.prmtop ligand.inpcrd")
        time.sleep(0.5)
        lig_crd = app.AmberInpcrdFile("ligand.inpcrd")
        positions = lig_crd.positions.value_in_unit(unit.angstroms)
        longest_distance = max(power(sum_a(power2(list(positions),2)),0.5))+padding
        box_x = max([abs(positions[i][0]) for i in range(len(positions))]) + padding
        box_y = max([abs(positions[i][1]) for i in range(len(positions))]) + padding
        box_z = max([abs(positions[i][2]) for i in range(len(positions))]) + padding
        return box_x, box_y, box_z, longest_distance
    

def ligand_energy():
        #self.command("saveamberparm ligand ligand.prmtop ligand.inpcrd")
        #time.sleep(1)
        from  mpmath import mp as math
        math.prec = 200
        lig_top = app.AmberPrmtopFile("ligand.prmtop")
        lig_crd = app.AmberInpcrdFile("ligand.inpcrd")
        lig_system = lig_top.createSystem(nonbondedMethod=NoCutoff, nonbondedCutoff=10*nanometer, constraints=HAngles, implicitSolvent=OBC1)
        lig_integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        lig_simulation = Simulation(lig_top.topology, lig_system, lig_integrator)
        lig_simulation.context.setPositions(lig_crd.positions)
        lig_state = lig_simulation.context.getState(getEnergy = True)
        lig_energy = lig_state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        return lig_energy
    
    
def get_PO3(positions_old, positions):
    pos = positions
    vec_a = (positions[len(positions_old)-1]-positions[len(positions_old)-2])
    x, y, z = vec_a.value_in_unit(unit.angstroms)
    x, y, z = np.array([x,y,z])/np.linalg.norm(np.array([x,y,z]))
    shift_forward = mm.Vec3(0,0,0)*unit.angstroms-positions[len(positions_old)-1]
    phi_2 = np.random.uniform(-np.math.pi/2,np.math.pi/2)
    s = np.math.sin(phi_2)
    c = np.math.cos(phi_2)
    rot = np.array([[2*(np.power(x,2)-1)*np.power(s,2)+1, 2*x*y*np.power(s,2)-2*z*c*s, 2*x*z*np.power(s,2)+2*y*c*s],
                    [2*x*y*np.power(s,2)+2*z*c*s, 2*(np.power(y,2)-1)*np.power(s,2)+1, 2*z*y*np.power(s,2)-2*x*c*s],
                    [2*x*z*np.power(s,2)-2*y*c*s, 2*z*y*np.power(s,2)+2*x*c*s, 2*(np.power(z,2)-1)*np.power(s,2)+1]])
    
    for j in range(len(positions_old)-1,len(positions)):
        pos[j] += shift_forward

    for j in range(len(positions_old)-1,len(positions)):
        #pos[j] += drift
        roted = np.dot(np.array(pos[j].value_in_unit(unit.angstrom)),rot)
        pos[j] = mm.Vec3(roted[0],roted[1],roted[2])*unit.angstrom
        pos[j] -= shift_forward
    
    positions_new = pos
    return positions_new

def get_PO5(positions_old, positions):
    pos = positions
    vec_a = (positions[len(positions_old)+2]-positions[len(positions_old)-1])
    x, y, z = vec_a.value_in_unit(unit.angstroms)
    x, y, z = np.array([x,y,z])/np.linalg.norm(np.array([x,y,z]))
    shift_forward = mm.Vec3(0,0,0)*unit.angstroms-positions[len(positions_old)+2]
    phi_2 = np.random.uniform(-np.math.pi/2,np.math.pi/2)
    s = np.math.sin(phi_2)
    c = np.math.cos(phi_2)
    rot = np.array([[2*(np.power(x,2)-1)*np.power(s,2)+1, 2*x*y*np.power(s,2)-2*z*c*s, 2*x*z*np.power(s,2)+2*y*c*s],
                    [2*x*y*np.power(s,2)+2*z*c*s, 2*(np.power(y,2)-1)*np.power(s,2)+1, 2*z*y*np.power(s,2)-2*x*c*s],
                    [2*x*z*np.power(s,2)-2*y*c*s, 2*z*y*np.power(s,2)+2*x*c*s, 2*(np.power(z,2)-1)*np.power(s,2)+1]])
    
    for j in range(len(positions_old)+2,len(positions)):
        pos[j] += shift_forward

    for j in range(len(positions_old)+2,len(positions)):
        #pos[j] += drift
        roted = np.dot(np.array(pos[j].value_in_unit(unit.angstrom)),rot)
        pos[j] = mm.Vec3(roted[0],roted[1],roted[2])*unit.angstrom
        pos[j] -= shift_forward
    
    positions_new = pos
    return positions_new

def get_C5O5(positions_old, positions):
    pos = positions
    vec_a = (positions[len(positions_old)+3]-positions[len(positions_old)+2])
    x, y, z = vec_a.value_in_unit(unit.angstroms)
    x, y, z = np.array([x,y,z])/np.linalg.norm(np.array([x,y,z]))
    shift_forward = mm.Vec3(0,0,0)*unit.angstroms-positions[len(positions_old)+3]
    phi_2 = np.random.uniform(-np.math.pi/2,np.math.pi/2)
    s = np.math.sin(phi_2)
    c = np.math.cos(phi_2)
    rot = np.array([[2*(np.power(x,2)-1)*np.power(s,2)+1, 2*x*y*np.power(s,2)-2*z*c*s, 2*x*z*np.power(s,2)+2*y*c*s],
                    [2*x*y*np.power(s,2)+2*z*c*s, 2*(np.power(y,2)-1)*np.power(s,2)+1, 2*z*y*np.power(s,2)-2*x*c*s],
                    [2*x*z*np.power(s,2)-2*y*c*s, 2*z*y*np.power(s,2)+2*x*c*s, 2*(np.power(z,2)-1)*np.power(s,2)+1]])
    
    for j in range(len(positions_old)+3,len(positions)):
        pos[j] += shift_forward

    for j in range(len(positions_old)+3,len(positions)):
        #pos[j] += drift
        roted = np.dot(np.array(pos[j].value_in_unit(unit.angstrom)),rot)
        pos[j] = mm.Vec3(roted[0],roted[1],roted[2])*unit.angstrom
        pos[j] -= shift_forward
    
    positions_new = pos
    return positions_new

def get_C5(positions_old, positions):
    pos = positions
    vec_a = (positions[len(positions_old)+6]-positions[len(positions_old)+3])
    x, y, z = vec_a.value_in_unit(unit.angstroms)
    x, y, z = np.array([x,y,z])/np.linalg.norm(np.array([x,y,z]))
    shift_forward = mm.Vec3(0,0,0)*unit.angstroms-positions[len(positions_old)+6]
    phi_2 = np.random.uniform(-np.math.pi/2,np.math.pi/2)
    s = np.math.sin(phi_2)
    c = np.math.cos(phi_2)
    rot = np.array([[2*(np.power(x,2)-1)*np.power(s,2)+1, 2*x*y*np.power(s,2)-2*z*c*s, 2*x*z*np.power(s,2)+2*y*c*s],
                    [2*x*y*np.power(s,2)+2*z*c*s, 2*(np.power(y,2)-1)*np.power(s,2)+1, 2*z*y*np.power(s,2)-2*x*c*s],
                    [2*x*z*np.power(s,2)-2*y*c*s, 2*z*y*np.power(s,2)+2*x*c*s, 2*(np.power(z,2)-1)*np.power(s,2)+1]])
    
    for j in range(len(positions_old)+6,len(positions)):
        pos[j] += shift_forward

    for j in range(len(positions_old)+6,len(positions)):
        #pos[j] += drift
        roted = np.dot(np.array(pos[j].value_in_unit(unit.angstrom)),rot)
        pos[j] = mm.Vec3(roted[0],roted[1],roted[2])*unit.angstrom
        pos[j] -= shift_forward
    
    positions_new = pos
    return positions_new

def get_base(positions_old, positions):
    pos = positions
    vec_a = (positions[len(positions_old)+11]-positions[len(positions_old)+9])
    x, y, z = vec_a.value_in_unit(unit.angstroms)
    x, y, z = np.array([x,y,z])/np.linalg.norm(np.array([x,y,z]))
    shift_forward = mm.Vec3(0,0,0)*unit.angstroms-positions[len(positions_old)+11]
    phi_2 = np.random.uniform(-np.math.pi/2,np.math.pi/2)
    s = np.math.sin(phi_2)
    c = np.math.cos(phi_2)
    rot = np.array([[2*(np.power(x,2)-1)*np.power(s,2)+1, 2*x*y*np.power(s,2)-2*z*c*s, 2*x*z*np.power(s,2)+2*y*c*s],
                    [2*x*y*np.power(s,2)+2*z*c*s, 2*(np.power(y,2)-1)*np.power(s,2)+1, 2*z*y*np.power(s,2)-2*x*c*s],
                    [2*x*z*np.power(s,2)-2*y*c*s, 2*z*y*np.power(s,2)+2*x*c*s, 2*(np.power(z,2)-1)*np.power(s,2)+1]])
    end = 0
    #print(len(positions)-len(positions_old))
    if len(positions)-len(positions_old) == 30:
        end = len(positions_old)+24
    elif len(positions)-len(positions_old) == 32:
        end = len(positions_old)+26
    elif len(positions)-len(positions_old) == 33:
        end = len(positions_old)+27
        
    for j in range(len(positions_old)+11,end-1):
        pos[j] += shift_forward

    for j in range(len(positions_old)+11,end-1):
        #pos[j] += drift
        roted = np.dot(np.array(pos[j].value_in_unit(unit.angstrom)),rot)
        pos[j] = mm.Vec3(roted[0],roted[1],roted[2])*unit.angstrom
        pos[j] -= shift_forward
    
    positions_new = pos
    return positions_new


    
#
#def join(partials):
#    energies = partials
#    return energies
    

def initial_sample(topology,coordinates,Nsteps,index,box=50,rang=[0,60]):
    print("Index is: ",index)
    aptamer_top = topology
    aptamer_crd = coordinates
    en = []
    xyz = []
    positions = []
    free_E_old = 1e50
    cnt = 0
    centre = np.math.ceil((rang[1]-rang[0])/2)
    #print("a")
    system = aptamer_top.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, implicitSolvent=app.OBC1)
    integrator = mm.LangevinIntegrator(300.*unit.kelvin, 1./unit.picosecond, 0.002*unit.picoseconds)
    simulation = app.Simulation(aptamer_top.topology, system, integrator)
    #print("b")
    for i in range(Nsteps):

        pos = aptamer_crd.positions
        pos0 = aptamer_crd.positions[int(centre)]
        shift = mm.Vec3(np.random.uniform(-box,box),np.random.uniform(-box,box),np.random.uniform(-box,box))*unit.angstrom
        x = np.random.uniform(-1,1)
        y = np.random.uniform(-1,1)
        z = np.random.uniform(-1,1)
        x, y, z = np.array([x,y,z])*1/(np.linalg.norm(np.array([x,y,z])))

        phi_2 = np.random.uniform(-np.math.pi,np.math.pi)
        x, y, z = np.array([x,y,z])*1/(np.linalg.norm(np.array([x,y,z])))
        xyz.append([shift[0].value_in_unit(unit.angstroms), shift[1].value_in_unit(unit.angstroms), shift[2].value_in_unit(unit.angstroms), x, y, z, phi_2])

        #phi_2 = 2*np.random.uniform(-np.math.pi,np.math.pi)
        s = np.math.sin(phi_2)
        c = np.math.cos(phi_2)
        rot = np.array([[2*(np.power(x,2)-1)*np.power(s,2)+1, 2*x*y*np.power(s,2)-2*z*c*s, 2*x*z*np.power(s,2)+2*y*c*s],
                        [2*x*y*np.power(s,2)+2*z*c*s, 2*(np.power(y,2)-1)*np.power(s,2)+1, 2*z*y*np.power(s,2)-2*x*c*s],
                        [2*x*z*np.power(s,2)-2*y*c*s, 2*z*y*np.power(s,2)+2*x*c*s, 2*(np.power(z,2)-1)*np.power(s,2)+1]])
        #print(rot)
        drift = get_aptamer(get_ligand_range(aptamer_top.topology), aptamer_crd.positions)[10]
        for j in range(get_ligand_range(aptamer_top.topology)[1],len(pos)):
            pos[j] -= drift
        for j in range(0,get_ligand_range(aptamer_top.topology)[1]):
            pos[j] -= pos0

        #pos[126] += drift
        for j in range(get_ligand_range(aptamer_top.topology)[1],len(pos)):
            #pos[j] += drift
            roted = np.dot(np.array(pos[j].value_in_unit(unit.angstrom)),rot)
            pos[j] = mm.Vec3(roted[0],roted[1],roted[2])*unit.angstrom
            pos[j] += shift 
        simulation.context.setPositions(pos)
        #print("minimizing ...")
        #simulation.minimizeEnergy()
        state = simulation.context.getState(getPositions=True,getEnergy=True,groups=1)
        free_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        #print(free_E)
        en.append(free_E)
        #monte_pos.append(shift)
        if free_E < free_E_old:
            free_E_old = free_E
            cnt += 1
            positions = pos
            #print(free_E)
	    fil = open("montetest%s.pdb"%i,"w")
	    app.PDBFile.writeModel(aptamer_top.topology,pos,file=fil,modelIndex=i)
	    fil.close()
	    del fil

    return en, pos, xyz, free_E, index



def stratified_sample(topology, coordinates, variances, Nsteps, index, box=50, rang=[0,400]):
    aptamer_top = topology
    aptamer_crd = coordinates
    en = []
    xyz = []
    cnt = 0
    for i in range(Nsteps):
        pos = aptamer_crd.positions
        pos0 = aptamer_crd.positions[(rang[1]-rang[0])/2]
        shift = [0,0,0]
        shift[0], shift[1], shift[2], x, y, z, phi_2 = uniform_strat(variances[0],variances[1],variances[2],variances[3],variances[4],variances[5],variances[6],[-box,box],[-np.math.pi,np.math.pi])
        x, y, z = np.array([x,y,z])*1/(np.linalg.norm(np.array([x,y,z])))
        xyz.append([shift[0], shift[1], shift[2], x, y, z, phi_2])
        shift = mm.Vec3(*shift)*unit.angstroms
        s = np.math.sin(phi_2)
        c = np.math.cos(phi_2)
        rot = np.array([[2*(np.power(x,2)-1)*np.power(s,2)+1, 2*x*y*np.power(s,2)-2*z*c*s, 2*x*z*np.power(s,2)+2*y*c*s],
                        [2*x*y*np.power(s,2)+2*z*c*s, 2*(np.power(y,2)-1)*np.power(s,2)+1, 2*z*y*np.power(s,2)-2*x*c*s],
                        [2*x*z*np.power(s,2)-2*y*c*s, 2*z*y*np.power(s,2)+2*x*c*s, 2*(np.power(z,2)-1)*np.power(s,2)+1]])
        #print(rot)
        drift = get_aptamer(get_ligand_range(aptamer_top.topology), aptamer_crd.positions)[10]
        for j in range(get_ligand_range(aptamer_top.topology)[1],len(pos)):
            pos[j] -= drift
        for j in range(0,get_ligand_range(aptamer_top.topology)[1]):
            pos[j] -= pos0
        #pos[126] += drift
        for j in range(get_ligand_range(aptamer_top.topology)[1],len(pos)):
            #pos[j] += drift
            roted = np.dot(np.array(pos[j].value_in_unit(unit.angstrom)),rot)
            pos[j] = mm.Vec3(roted[0],roted[1],roted[2])*unit.angstrom
            pos[j] += shift 
        simulation.context.setPositions(pos)
        #print("minimizing ...")
        #simulation.minimizeEnergy()
        state = simulation.context.getState(getPositions=True,getEnergy=True,groups=1)
        #print("getting positions ...")
        #fil = open("montetest%s.pdb"%i,"w")
        #app.PDBFile.writeModel(aptamer_top.topology,pos,file=fil,modelIndex=i)
        #fil.close()
        #del fil
        #simulation.step(1)
        free_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        en.append(free_E)
        #monte_pos.append(shift)

        if free_E < free_E_old:
            free_E_old = free_E
            cnt += 1
            positions = pos
            
        return en, pos, xyz, free_E, index


def mcmc_sample(topology, coordinates, old_coordinates, index, stepsize=200, Nsteps=5000, ligand_heavyness=1e50, aptamer_heaviness=1e3):
    aptamer_top = topology
    aptamer_crd = coordinates
    pos = old_coordinates
    #print(len(pos))
    #print(len(pos)-len(coordinates.positions))
    #*unit.angstroms
    len_heavy = len(pos)-2-len(get_ligand(aptamer_top.topology))
    system = aptamer_top.createSystem(nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=1.2*unit.nanometers, constraints=app.HBonds, implicitSolvent=app.OBC1)
    #print(index,index,index,index)
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
    simulation = app.Simulation(aptamer_top.topology, system, integrator)
    en = []
    xyz = []
    positions = []
    free_E_old = 1e20
    simulation.context.setPositions(get_C5(pos, get_C5O5(pos, get_base(pos, get_PO5(pos,get_PO3(pos, position_aptamer(pos, aptamer_crd.positions)))))))

    state = simulation.context.getState(getPositions=True,getEnergy=True,groups=1)

    for i in range(Nsteps):
        simulation.context.setPositions(get_base(pos, get_C5(pos, get_C5O5(pos, get_PO5(pos,get_PO3(pos, position_aptamer(pos, aptamer_crd.positions)))))))
        state = simulation.context.getState(getPositions=True,getEnergy=True,groups=1)
        free_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        #print(free_E)
        fil = open("montestep%s.pdb"%i,"w")
	app.PDBFile.writeModel(aptamer_top.topology,state.getPositions(),file=fil,modelIndex=i)
	fil.close()
	del fil

        if free_E < free_E_old:
            positions = state.getPositions()
        en.append(free_E)
        #print(free_E)

    return en, positions, free_E



def mcmc_sample_five(topology, coordinates, old_coordinates, index, stepsize=200, Nsteps=20, ligand_heavyness=1e50, aptamer_heaviness=1e3):
    aptamer_top = topology
    aptamer_crd = coordinates
    pos = np.array(old_coordinates)
    heavies = [len(get_ligand(aptamer_top.topology))+(len(positions)-len(positions_old))-1,len(aptamer_crd.coordinates)]
    system = aptamer_top.createSystem(nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=1.2*unit.nanometers, constraints=app.HBonds, implicitSolvent=app.OBC1)

    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
    simulation = app.Simulation(aptamer_top.topology, system, integrator)
    en = []
    xyz = []
    positions = []

    simulation.context.setPositions(position_aptamer_five(pos, aptamer_crd.positions))
    simulation.minimizeEnergy(maxIterations=5000)
    state = simulation.context.getState(getPositions=True,getEnergy=True,groups=1)
    free_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    en.append(free_E)
    posit = state.getPositions()
    integrator = mm.LangevinIntegrator(5000*unit.kelvin, 1.0/unit.picoseconds, 1.0*unit.femtoseconds)
    system = aptamer_top.createSystem(nonbondedMethod=app.NoCutoff, constraints=app.HBonds, implicitSolvent=app.OBC1)
    constrainPO3(aptamer_top.topology, system)
    constrainPO5(aptamer_top.topology, system)
    constrainC5O5(aptamer_top.topology, system)
    for prt in get_ligand(aptamer_top.topology):
        system.setParticleMass(prt,1e50)
    for prt in range(heavies[0],heavies[1]):
        system.setParticleMass(prt,1e3)
    simulation = app.Simulation(aptamer_top.topology, system, integrator)

    simulation.context.setPositions(posit)

    for i in range(Nsteps):
        simulation.step(stepsize)
        state = simulation.context.getState(getPositions=True,getEnergy=True,groups=1)
        free_E = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        if free_E < free_E_old:
            positions = append(state.getPositions())
        en.append(free_E)

    return en, positions, free_E
    
        

def initial(Ntide):
    global _NINIT
    global _INFILE
    global _BETA
    beta = _BETA
    print("Constructing Ligand/Aptamer complex ...")
    
    internal = Aptamer("leaprc.ff12SB", _INFILE)
    internal.sequence(Ntide,Ntide)
    internal.unify(Ntide)
    internal.command("saveamberparm union %s.prmtop %s.inpcrd"%(Ntide,Ntide))
    
    time.sleep(1)
    
    print("Aptamer/Ligand complex constructed.")
    print("Loading Aptamer/Ligand complex ...")
    
    aptamer_top = app.AmberPrmtopFile("%s.prmtop"%Ntide)
    aptamer_crd = app.AmberInpcrdFile("%s.inpcrd"%Ntide)
    
    ligand_range = get_ligand_range(aptamer_top.topology)
    #sample_box_task = ligand_box(0.1)
    #sample_box = bvalue(sample_box_task)[3]
    sample_box = 5
    #print(sample_box)
    volume = (2*sample_box)**3*(2*math.pi)**3
    print("Sampling parameter space ...")
    print(Ntide)
    en_pos_xyz = [initial_sample(aptamer_top, aptamer_crd, _NINIT, i, box=sample_box, rang=ligand_range) for i in range(100)]
    #en_pos_xyz = bvalue(en_pos_xyz_tasks)
    
    print("done.")
    print("Harvesting results ...")
    en = []
    xyz = []
    positions_s = []
    positions = []
    for elem in en_pos_xyz:
        en += elem[0]
        positions_s.append([elem[3], elem[1]])
        xyz += elem[2]
    positions = min(positions_s)[1]
    #print(positions)
    #print("Stratifying and bootstrapping ...")
    #
    #variances = stratify(xyz, en, [-sample_box,sample_box], [-np.math.pi,np.math.pi])
    #
    #en_pos_xyz_strat = [bvalue(stratified_sample(variances, aptamer_top, aptamer_crd, 10000, i, box=sample_box, rang=ligand_range)) for i in range(10)]
    #
    #for elem in en_pos_xyz_strat:
    #    en += elem[0]
    #    positions_s += (elem[3], elem[1])
    #    xyz += elem[2]
        
    Z = volume*sum([math.exp(-beta*elem) for elem in en])/len(en)
    #print(Z)
    P = [math.exp(-beta*elem)/Z for elem in en]
    #print(P)
    S = volume*sum([-elem*math.log(elem*volume) for elem in P])/len(P)
    print("Ntide: %s entropy: %s"%(Ntide,S))
    
    return positions, Ntide, S


def evaluate(positions, Ntides, entropies, threshold=0.1):
    res_Ntide_positions = []
    res_positions = []
    print("Chosen ntides with entropies:")
    for pos, alem, blem in zip(positions, Ntides, entropies):
        if blem <= min(entropies)+threshold:
            res_Ntide_positions.append([pos,alem])
            #res_positions.append(pos)
            print("%s : %s"%(alem,blem))
    return res_Ntide_positions


def step(array):
    global _INFILE
    global _BETA
    global _NSTEP
    beta = _BETA
    old_positions, Ntides, is_3prime = array
    internal = Aptamer("leaprc.ff12SB",_INFILE)
    identifier = Ntides.replace(" ","")
    internal.sequence(identifier,Ntides.strip())
    internal.unify(identifier)
    internal.command("saveamberparm union %s.prmtop %s.inpcrd"%(identifier,identifier))
    time.sleep(2)
    
    #print("WhereamI?")
    print("Identifier: "+Ntides)

    volume = (2*math.pi)**5
    aptamer_top = app.AmberPrmtopFile("%s.prmtop"%identifier)
    aptamer_crd = app.AmberInpcrdFile("%s.inpcrd"%identifier)
    
#    print("loaded")
    
#    if is_3prime == 1:
    en_pos = [mcmc_sample(aptamer_top, aptamer_crd, old_positions, index, Nsteps=_NSTEP) for index in range(10)]
#    else:
#        en_pos_task = [mcmc_sample_five(aptamer_top, aptamer_crd, old_positions, index, Nsteps=200) for index in range(20)]
#    barrier()
#    en_pos = value(en_pos_task)
    en = []
    positions = []
    positions_s = []
    for elem in en_pos:
        en += elem[0]
        #print(elem[2], elem[1])
        positions_s.append([elem[2], elem[1]])
    
    positions = min(positions_s)[1]
    
    fil = open("best_structure%s.pdb"%Ntides,"w")
    app.PDBFile.writeModel(aptamer_top.topology,positions,file=fil)
    fil.close()
    del fil
    
    Z = volume*math.fsum([math.exp(-beta*elem) for elem in en])/len(en)
    P = [math.exp(-beta*elem)/Z for elem in en]
    S = volume*math.fsum([-elem*math.log(elem*volume) for elem in P])/len(P)
    
    print("%s : %s"%(Ntides,S))
    
    return positions, Ntides, S


def choice():
    pass

#print("Initializing with alphabet ...")
    
#alphabet = ["DGN","DAN","DTN","DCN"]
    
#print(alphabet)
#print("Choosing from candidates ...")
    
#pos_Nt_S = []
#for elem in alphabet:
#    pos_Nt_S += [initial(elem)]
    #pos_Nt_S = bvalue(pos_Nt_S_task)
    
#print("Nucleotides with their respective entropies are: ")
#print(pos_Nt_S[1],pos_Nt_S[2])
    
#positions = [elem[0] for elem in pos_Nt_S]
#Ntides = [elem[1] for elem in pos_Nt_S]
#entropies = [elem[2] for elem in pos_Nt_S]
#pos_Nt = evaluate(positions, Ntides, entropies, threshold=0.1)

def loop():
    global _THRESHOLD
    global _NMER
    #print("Initializing with alphabet ...")
    
    alphabet = ["DGN","DAN","DTN","DCN"]
    
    print(alphabet)
    print("Choosing from candidates ...")
    
    pos_Nt_S_task = []
    pool = Pool(4)
    pos_Nt_S_task = pool.map(initial,alphabet)

    pos_Nt_S = pos_Nt_S_task
    
    #print("Nucleotides with their respective entropies are: ")
    #print(pos_Nt_S[1],pos_Nt_S[2])
    
    positions = [elem[0] for elem in pos_Nt_S]
    #print(pos_Nt_S)
    
    Ntides = [elem[1] for elem in pos_Nt_S]
    entropies = [elem[2] for elem in pos_Nt_S]

    pos_Nt = evaluate(positions, Ntides, entropies, threshold=0.5)
    
    #barrier()
    positions = []
    Ntides = []
    for elem in pos_Nt:
        positions.append(elem[0])
        Ntides.append(elem[1])
    positions = positions
    print([len(elem) for elem in positions])
    #sleep(20)
    
    print("Chosen nucleotides: ")
    print(Ntides)
    
    for i in range(_NMER):
        
        if i == 0:
            print("Initializing 2nd step ...")
        elif i == 1:
            print("Initializing 3rd step ...")
        else:
            print("Initializing %sth step ..."%(i+2))
        
        pos_Nt_S = []
        pos_Nt_S_tasks = []
        pool = Pool(6)
        #print(zip(positions, Ntides))
        pos_Nt_S_tasks = pool.map(step,[[alem,blem.replace("3","").replace("N","5").strip()+ntide,1] for alem, blem in zip(positions,Ntides) for ntide in [" DG3"," DA3"," DT3"," DC3"]])            
                
        pos_Nt_S = pos_Nt_S_tasks
        positions = [elem[0] for elem in pos_Nt_S]
        positions = positions
        #print(positions)
        Ntides = [elem[1] for elem in pos_Nt_S]
        entropies = [elem[2] for elem in pos_Nt_S]
        
        #print("Sequences with entropies are: ")
        #print(pos_Nt_S[1], pos_Nt_S[2])
        
        pos_Nt_task = evaluate(positions, Ntides, entropies, threshold=0.005)
        #barrier()
        pos_Nt = pos_Nt_task
        Ntides = []
	positions = []
        for elem in pos_Nt:
            Ntides.append(elem[1])
        for elem in pos_Nt:
	    positions.append(elem[0])
        print("Chosen sequences are: ")
        for elem in pos_Nt:
            print(elem[1])
        
    return pos_Nt


def result(pos_Nt):
    """Write a lot of files"""
    count = 0
    for elem in pos_Nt:
        internal = Aptamer("leaprc.ff12SB",_INFILE)
        identifier = elem[1].replace(" ","")
        internal.sequence(identifier, elem[1])
        internal.unify(identifier)
        internal.command("saveamberparm union Aptamer%s.prmtop Aptamer%s.inpcrd"%(count,count))
        time.sleep(2)
        #aptamer_top = app.AmberPrmtopFile("casGC.prmtop")
        #aptamer_crd = app.AmberInpcrdFile("casGC.inpcrd")
    print("Run successful! Have fun with your Aptamers")
    return 1

parser = argparse.ArgumentParser(description='MAWS - Make Aptamers Without SELEX', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('path', metavar='PATH', help='Path to the calculation directory.')
parser.add_argument('infile', metavar='INFILE', help='3D structure of the ligand.')
parser.add_argument('-b', '--beta', type=float, default=0.01, help='lagrange multiplier β.')
parser.add_argument('-i', '--ninit', type=int, default=200, help='number of initial steps as multiple of 100.')
parser.add_argument('-s', '--nstep', type=int, default = 200, help='number of samples in every step after the first.')
parser.add_argument('-l', '--nmer', type=int, default = 15, help='The final length of the aptamer.')
parser.add_argument('-t', '--threshold', type=int, default = 0.01, help='Threshold for candidate selection.')
parser.add_argument('-f', '--format', type=str, default = "pdb", help='input file format, may be one of pdb or mol2.')
parser.add_argument('-h', '--hybrid', type=str, default = "", help='parameter modifying file for hybrid calculations')


args = parser.parse_args()

_INFILE = args.infile
_BETA = args.beta
_NINIT = args.ninit
_NSTEP = args.nstep
_NMER = args.nmer
_THRESHOLD = args.threshold
_HYBRID = args.hybrid
if args.format in ['pdb','mol2']:
    _FORMAT = args.format
else:
    raise ValueError("INVALID FORMAT %s"%args.format)
    

positions_and_Ntides = loop()
result(positions_and_Ntides)

print("Run successful!")
print("Please come again!")


