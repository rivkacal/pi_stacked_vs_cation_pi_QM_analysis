'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Details:
This Module contatins all planes creation and distances/angles calculations, including rotation matrices and assignment according to specific residue atom-level details!
Also pre-analysis (geometrical parameters calculations) and post-analysis categorization of interactions are provided. 
Some functions include inside flags. If you modify the code make sure that they do not exclude the case of your interest.
# ! look for all # ! warnings before running these functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
import numpy as np
import math
import LoadStruct as lst
from math import sqrt
from collections import namedtuple
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Bio.PDB import Superimposer
from Bio.PDB.Residue import DisorderedResidue, Residue
from Bio.PDB.vectors import Vector
from termcolor import cprint
from Bio.PDB import PDBIO, PDBList
from Bio.Seq import Seq
from Bio.PDB.Polypeptide import Polypeptide
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain

Point = namedtuple("Point", ["x", "y"])

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SIMPLE FORMATTING/DATA STRUCTURE FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
def flatten_points(point):
    '''This is neccessary working with biopython returned atoms positions'''
    if point.shape == (1,3):
        return point.flatten()
    return point

def validate_atoms_existence(residue, atoms_list) -> bool:
    ''' checks if all atom types in input list are present in the current residue structure'''
    for atom in atoms_list:
        is_current_atom_present = residue.has_id(atom)
        if not is_current_atom_present:
            return False
    return True


'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ANGLES/PLANES CALCULATIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

def calc_res_center(shape, pt1, pt2, pt3, pt4, pt5, pt6) -> float:
    ''' This function recieves current residue shape, and atom coordinates coords up to 6
        points to calculate the centroid according to the polygon shape'''
    
    points = []
    centroid = np.zeros((1,3))

    if shape == 'hexagon':
        points = [pt1, pt2, pt3, pt4, pt5, pt6]
                    
    if shape == 'pentagon':
        points = [pt1, pt2, pt3, pt4, pt5]
    
    for idx,point in enumerate(points):
        points[idx] = flatten_points(point)
        
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    z = [p[2] for p in points]
    centroid = (sum(x) / len(points), sum(y) / len(points), sum(z) / len(points))
    return centroid

def vec_plane_angle(vec1, a: float, b: float, c: float) :
    ''' This function gets one plane coefficents and another vector and returns the angle between the plane normal and the given vector '''
    v1 = vec1
    v2 = np.array([a, b, c])
    if (np.linalg.norm(v1)*np.linalg.norm(v2)) < 0.000001:
        return None
    cosalpha = np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
    angle = np.arccos(np.clip(cosalpha, -1, 1))
    angle_deg = angle * 180 / np.pi
    if angle_deg > 90:
        angle_deg = 180 - angle_deg # return the smallest angle as plane normal could point difeferently
    return angle_deg

def two_pt_dist(point1, point2):
    ''' This function returns the Euclidean distance between two input points '''
    # calculating Euclidean distance
    # using linalg.norm()
    dist = np.linalg.norm(point1 - point2)
    
    return dist
                
def two_plane_angle(a1: float,b1: float,c1: float,a2: float,b2: float,c2: float) -> float:
    ''' This function gets two different planes coefficents and returns the angle between them '''
    # say v1 is normal to plane 1 such that: v1=(a1,b1,c1)
    # v2 is normal to the second plane: v2=(a2,b2,c2)
    # cos(alpha) = v1*v2/(|v1||v2|)
    v1 = np.array([a1, b1, c1])
    v2 = np.array([a2, b2, c2])
    if (np.linalg.norm(v1)*np.linalg.norm(v2)) < 0.000001:
        return None
    cosalpha = np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
    #angle = np.arccos(np.clip(cosalpha, -1, 1))
    return cosalpha

#from: https://stackoverflow.com/questions/53591350/plane-fit-of-3d-points-with-singular-value-decomposition
def plane_svd_fit(array_of_coords, plot_flag):
    '''recieves np array of 3D points (N*3) and returns best fitted a,b,c,d for ax+by+cz+d=0'''
    xyz = array_of_coords
    if xyz.shape == (0,) or xyz.shape == (3,) or xyz is None:
        return 0,0,0,0
    ''' best plane fit'''
    #1.calculate centroid of points and make points relative to it
    centroid         = xyz.mean(axis = 0)
    xyzT             = np.transpose(xyz)
    xyzR             = xyz - centroid                         #points relative to centroid
    #xyzRT            = np.transpose(xyzR)                       

    #2. calculate the singular value decomposition of the xyzT matrix and get the normal as the last column of u matrix
    u, sigma, v       = np.linalg.svd(xyzR)
    normal            = v[2]                                 
    normal            = normal / np.linalg.norm(normal)       #we want normal vectors normalized to unity

    '''matplotlib display'''
    #prepare normal vector for display
    forGraphs = list()

    centroid = np.array(centroid)
    centroid = centroid.flatten()
    forGraphs.append(np.array([centroid[0],centroid[1],centroid[2],normal[0],normal[1], normal[2]]))

    #get d coefficient to plane for display
    d = normal[0] * centroid[0] + normal[1] * centroid[1] + normal[2] * centroid[2]
        
    if plot_flag:
        # create x,y for display
        minPlane = int(math.floor(min(min(xyzT[0]), min(xyzT[1]), min(xyzT[2]))))
        maxPlane = int(math.ceil(max(max(xyzT[0]), max(xyzT[1]), max(xyzT[2]))))
        xx, yy = np.meshgrid(range(minPlane,maxPlane), range(minPlane,maxPlane))

        # calculate corresponding z for display
        z = (-normal[0] * xx - normal[1] * yy + d) * 1. /normal[2]

        #matplotlib display code
        forGraphs = np.asarray(forGraphs)
        X, Y, Z, U, V, W = zip(*forGraphs)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        #ax.plot_surface(xx, yy, z, alpha=0.2)
        ax.scatter(xyzT[0],xyzT[1],xyzT[2])
        ax.quiver(X, Y, Z, U, V, W)
        ax.set_xlim([min(xyzT[0])- 0.1, max(xyzT[0]) + 0.1])
        ax.set_ylim([min(xyzT[1])- 0.1, max(xyzT[1]) + 0.1])
        ax.set_zlim([min(xyzT[2])- 0.1, max(xyzT[2]) + 0.1])
        plt.xlabel('x')
        plt.ylabel('y')
        #plt.zlabel('z')
        ax.scatter(xyz[0,0], xyz[0,1],  xyz[0,2], c='red')
        #if the CG is the first in list of ring only-atms then it will be colored red, otherwise - ignore this color (treat as blue)
        ax.scatter(0,0,0, c='green') #origin for convinency
        plt.show()
    return normal[0], normal[1], normal[2], -d
    
def get_rotation_matrix(axis, theta):
    """
    Find the rotation matrix associated with counterclockwise rotation
    about the given axis by theta radians.
    Credit: http://stackoverflow.com/users/190597/unutbu

    Args:
        axis (list): rotation axis of the form [x, y, z]
        theta (float): rotational angle in radians

    Returns:
        array. Rotation matrix.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def rot_plane_to_XY(a,b,c):
    '''To create a reference residue first need to rotate to the XY plane (Then in the XY plane rotate about Z to align X direction)'''
    common_sqrt = np.sqrt(a**2 + b**2 + c**2)
    rot_mat = np.array([[b**2/(a**2+b**2)+(1-b**2/(a**2+b**2))*c/common_sqrt, -a*b*(1-c/common_sqrt), -a/common_sqrt],
                        [-a*b*(1-c/common_sqrt),a**2/(a**2+b**2)+(1-a**2/(a**2+b**2))*c/common_sqrt, -b/common_sqrt],
                        [-a/common_sqrt,b/common_sqrt, c/common_sqrt]])
    
    return rot_mat

def get_angle(p1: Point, p2: Point) -> float:
    """Get the angle of this line with the horizontal axis."""
    dx = p2.x - p1.x
    dy = p2.y - p1.y
    theta = math.atan2(dy, dx)
    angle = math.degrees(theta)  # angle is in (-180, 180]
    if angle < 0:
        angle = angle + 360
    return angle

def angle_between_3_points(a, b, c):
    '''Three points angle, where the angle is at the b position! '''
    ab = np.array(b) - np.array(a)
    bc = np.array(b) - np.array(c)
    dot_product = np.dot(ab, bc)
    magnitude_ab = np.linalg.norm(ab)
    magnitude_bc = np.linalg.norm(bc)
    angle_rad = np.arccos(dot_product / (magnitude_ab * magnitude_bc))
    angle_deg = angle_rad * 180 / np.pi
    return angle_deg

def get_superimpose_mat(refrence_res, res_to_rotran):
    'input: refrence residue and a given residue to rotate and translate on, return matrix'
    #use superimposer and vector/matrix operations!

    sup = Superimposer()
    # Specify the atom lists
    # 'fixed' and 'moving' are lists of Atom objects
    # The moving atoms will be put on the fixed atoms

    fixed = refrence_res.get_list()
    moving = res_to_rotran.get_list()
    #print(moving)
    fixed_list = []
    for item in fixed:
        if not item.name.startswith('H') and not item.name.startswith('CA') and not item.name.startswith('CB'):
            fixed_list.append(item)
    moving_list = []
    for item in moving:
        if not item.name.startswith('H') and not item.name.startswith('CA') and not item.name.startswith('CB'):
            moving_list.append(item)

    sup.set_atoms(fixed_list, moving_list)
    return sup.rotran

def apply_rotation_translation(rotran, residue):
    '''input rotation_matrix+translation_vector appended arrays'''
    residue_transformed = residue.copy()

    rotation = rotran[0]
    translation = rotran[1]

    for atom in residue_transformed:
        coordi2 = np.dot(atom.get_coord(), rotation) + translation

        atom.set_coord(coordi2)

    return residue_transformed


'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IDENTIFYING ATOMS FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''


def assign_side_chain_atoms(residue) -> list:
    ''' returns a list of all side chain atoms for the given residue, only if all the side chain atoms appear in structure. Else returns an empty list'''
    atoms_list = []    
    
    if residue.resname == "PRO":
        atoms_list = ["CA","CB", "CG", "CD", "N"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []            
                    
    if residue.resname == "GLY":
        atoms_list = []
        
    if residue.resname == "ALA":
        atoms_list = ["CB"]
        if not residue.has_id("CB"):
            atoms_list = []
            
    if residue.resname == "ARG":
        atoms_list = ["CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []            
              
    if residue.resname == "ASN":
        atoms_list = ["CB", "CG", "ND2"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []            
              
    if residue.resname == "ASP":
        atoms_list = ["CB", "CG", "OD1", "OD2"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []   
            
    if residue.resname == "CYS":
        atoms_list = ["CB", "SG"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = [] 
            
    if residue.resname == "GLN":
        atoms_list = ["CB", "CG", "CD", "NE2", "OE1"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = [] 
            
    if residue.resname == "GLU":
        atoms_list = ["CB","CG","CD","OE1","OE2"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []  
            
    if residue.resname == "HIS":
        atoms_list = ["CB","CG","CD2","ND1","NE2","CE1"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = [] 
            
    if residue.resname == "ILE":
        atoms_list = ["CB","CG1","CD1","CG2"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = [] 

    if residue.resname == "LEU":
        atoms_list = ["CB","CG","CD1","CD2"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = [] 

    if residue.resname == "LYS":
        atoms_list = ["CB","CG","CD","CE","NZ"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = [] 
            
    if residue.resname == "MET":
        atoms_list = ["CB","CG","SD","CE"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []
            
    if residue.resname == "PHE":
        atoms_list = ["CB","CG","CD1","CE1","CZ","CE2","CD2"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = [] 

    if residue.resname == "SER":
        atoms_list = ["CB","OG"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []     
    
    if residue.resname == "THR":
        atoms_list = ["CB","OG1","CG2"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []      

    if residue.resname == "TRP":
        atoms_list = ["CB","CG","CD1","NE1","CE2","CZ2","CH2","CZ3","CE3","CD2"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []  
            
    if residue.resname == "TYR":
        atoms_list = ["CB","CG","CD1","CE1","CZ","CE2","CD2"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []  
            
    if residue.resname == "VAL":
        atoms_list = ["CB","CG2","CG1"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []  
            
    return atoms_list

def assign_ring_atoms(residue) -> list:
    ''' returns a list of all rings' atoms for the given residue, only if the residue has a ring (excluding Pro! only aromatic). Else returns an empty list
    Note that for TRP only the 6 membered ring is considered'''
    atoms_list = []    
                           
    if residue.resname == "HIS":
        atoms_list = ["CG","CD2","ND1","NE2","CE1"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []             
            
    if residue.resname == "PHE":
        atoms_list = ["CG","CD1","CE1","CZ","CE2","CD2"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []  
    # TODO implement a flag for 5 atoms ring of TRP
    if residue.resname == "TRP":
        atoms_list = ["CD2","CE2","CZ2","CH2","CZ3","CE3"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []  
            
    if residue.resname == "TYR":
        atoms_list = ["CG","CD1","CE1","CZ","CE2","CD2"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []  

    return atoms_list

def is_aromatic_pair(residue_i,residue_j, thresold: float) -> bool:
    ''' check if any given residue in a pair has at least one ring (non-Hydrogen)atom within threshold from an atom in its pair sidechain's ring'''
    res_i_ring = assign_ring_atoms(residue_i)
    res_j_ring = assign_ring_atoms(residue_j)
    
    if res_i_ring != [] and res_j_ring != []:
        for atom_i in res_i_ring:
            for atom_j in res_j_ring:
                dist = two_pt_dist(residue_i[atom_i].get_coord(),residue_j[atom_j].get_coord())
                if dist <= thresold:
                    return True
    return False

def is_neigbor_within_thresh(residue_i,residue_j, thresold: float) -> bool:
    '''check if any of the sidechain's atom is within this threshold from its pair's sidechain (does not have to be aromatic)'''
    
    res_i_sidechain = assign_side_chain_atoms(residue_i)
    res_j_sidechain = assign_side_chain_atoms(residue_j)
    
    if res_i_sidechain != [] and res_j_sidechain != []:
        for atom_i in res_i_sidechain:
            for atom_j in res_j_sidechain:
                dist = two_pt_dist(residue_i[atom_i].get_coord(),residue_j[atom_j].get_coord())
                if dist <= thresold:
                    return True
    return False

def assign_cationpi_atoms(residue) -> list:
    '''returns a list of all given side chain atoms for the given residue that can participate as positive charge center, only if these side chain atoms appear in structure. 
    Else returns an empty list.
    Note Histidine was not discussed in the literature before as cation and the center atom of positive charge is not well defined'''
    atoms_list = []

    if residue.resname == "ARG":
        atoms_list = ["CZ"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []
    if residue.resname == "LYS":
        atoms_list = ["NZ"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []
    if residue.resname == "HIS":
        # treat both cases of Nitrogens seperately as LYS
        atoms_list = ["ND1","NE2"]
        if not validate_atoms_existence(residue, atoms_list):
            atoms_list = []

    return atoms_list

def assign_res_center(residue):
    ''' This function returns the position of the aromatic/imidazole ring center according to the residue input'''
    #initialize in case there is residue with not all below atoms observed in X-ray:
    centroid = (-10000, -10000, -10000)
    
    if residue.resname == 'TRP':
        shape = 'hexagon'
        if residue.has_id("CD2") and residue.has_id("CE3") and residue.has_id("CZ2") and residue.has_id("CZ3") and residue.has_id("CH2") and residue.has_id("CE2"):
            centroid = calc_res_center(shape, residue['CD2'].get_coord(), residue['CE3'].get_coord(), residue['CZ2'].get_coord(), residue['CZ3'].get_coord(), residue['CH2'].get_coord(), residue['CE2'].get_coord())
                                           
    if residue.resname == 'TYR' or residue.resname == 'PHE':
        shape = 'hexagon'
        if residue.has_id("CG") and residue.has_id("CD1") and residue.has_id("CE1") and residue.has_id("CZ") and residue.has_id("CE2") and residue.has_id("CD2"):
            centroid = calc_res_center(shape, residue['CG'].get_coord(), residue['CD1'].get_coord(), residue['CE1'].get_coord(), residue['CZ'].get_coord(), residue['CE2'].get_coord(), residue['CD2'].get_coord()) 
                                           
    if residue.resname == 'HIS':
        shape = 'pentagon'
        if residue.has_id("CG") and residue.has_id("ND1") and residue.has_id("CE1") and residue.has_id("NE2") and residue.has_id("CD2"):
            centroid = calc_res_center(shape, residue['CG'].get_coord(), residue['ND1'].get_coord(), residue['CE1'].get_coord(), residue['NE2'].get_coord(), residue['CD2'].get_coord(), None) 
        
    if np.sum(centroid) < -29000:
        centroid = None
    
    return centroid


'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CALCULATING PARMETERS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
#################################################################### CATION PI #################################################################
def cation_pi_dist(residue_i, residue_j, thresold: float):
    ''' return the cation-pi parameters defined in the paper: D,theta1,theta2. 
    first residue is aromatic acceptor (residue_i) is either HSE,PHE,TYR,TRP
    second is cation (residue_j) is either HSP,LYS,ARG
    threhold was not yet implemented (to screen some strctures), still and argument to avoid errors of outer functions calling this.
    !! This function currently excludes HIS-HIS pairs (HSP-HSE) when the inside flag: is_HIS_HIS_pair == False
    Please find the angles according to:
    J Mol Biol. 2021 Aug 20; 433(17): 167035.   doi: 10.1016/j.jmb.2021.167035
    or use the following link:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8338773/#:~:text=A%20range%20of%20distance%20thresholds,the%20molecules%20involved%20%5B16%5D'''
    # ! important to check this flag before running
    is_HIS_HIS_pair = False

    res_i_centroid = assign_res_center(residue_i)
    # the postion of the 'interesting' atom or the center of positive charge will be stored here:
    res_j_atom = assign_cationpi_atoms(residue_j)

    theta1 = None
    theta2 = None
       
    if res_j_atom == [] or res_i_centroid is None:        
        if residue_i.resname=='HIS' and residue_j.resname=='HIS':
            distance = [10000,10000]
            return distance, theta1, theta2
        # not two of them are HIS so we have 1 distance and not 2 values (array) in return
        distance = 10000
        return distance, theta1, theta2

    # if the cation is not HSP (HIS) check separately ARG,LYS. The distance and theta1 are similarly calculated
    if not residue_j.resname == 'HIS':
        # D parameter's magnitude
        distance = two_pt_dist(res_i_centroid, residue_j[res_j_atom[0]].get_coord())
        # to calculate aromatic plane first choose all heavy atoms of rings
        mat = []
        ring_atoms = assign_ring_atoms(residue_i)
        for atom in ring_atoms:
            mat.append(residue_i[atom].get_coord().flatten())
        mat = np.array(mat)
        # find the best plane fitting all 5/6 points (more than 3...) using svd algorithm found in stackoverflow: 
        a,b,c,d = plane_svd_fit(mat, False)
        plane_normal = np.array([a,b,c])
        # in case the above normal was assigned in the 'wrong' direction yielding angle > 90 later, calculate the opposite one (50% statistically):
        plane_normal2 = np.array([-a,-b,-c])
        plane_normal_vec = Vector(plane_normal)
        plane_normal_vec2 = Vector(plane_normal2)
        d_vec = residue_j[res_j_atom[0]].get_coord().flatten() - res_i_centroid
        d_vec = Vector(d_vec)
        angle_x = d_vec.angle(plane_normal_vec)
        theta1 = math.degrees(angle_x)
        if not theta1 <= 90:
            angle_y = d_vec.angle(plane_normal_vec2)
            theta1 = math.degrees(angle_y)
            if not theta1 <= 90:
                raise ValueError('angle Theta1 does not have a physical value [0,90]')

        # only for ARG we have additional angle theta2 that tells about the orientation of itsplane from the dirction of d_vec
        if residue_j.resname == 'ARG':
            arg_plane_atoms = ['NH1','NH2','NE','CZ']
            mat2 = []
            for atom in arg_plane_atoms:
                if residue_j.has_id(atom):
                    mat2.append(residue_j[atom].get_coord().flatten())
            if mat2 != []:
                mat2 = np.array(mat2)
            arg_a,arg_b,arg_c,arg_d = plane_svd_fit(mat2,False)
            if (arg_a==0 and arg_b==0 and arg_c==0 and arg_d==0):
                return distance, theta1, theta2
            arg_normal = np.array([arg_a,arg_b,arg_c])
            arg_normal2 = np.array([-arg_a,-arg_b,-arg_c])
            arg_normal = Vector(arg_normal)
            arg_normal2 = Vector(arg_normal2)
            angle2_x = d_vec.angle(arg_normal)
            theta2 = math.degrees(angle2_x)
            if not theta2 <= 90:
                angle2_y = d_vec.angle(arg_normal2)
                theta2 = math.degrees(angle2_y)
                if not theta2 <= 90:
                    raise ValueError('Theta2 does not have a physical value [0,90]')
                        

    if residue_j.resname == 'HIS':
        # planes calculations, aromatic acceptor:
        mat = []
        ring_atoms = assign_ring_atoms(residue_i)
        for atom in ring_atoms:
            mat.append(residue_i[atom].get_coord().flatten())
        mat = np.array(mat)
        a,b,c,d = plane_svd_fit(mat, False)

        plane_normal = np.array([a,b,c])
        plane_normal2 = np.array([-a,-b,-c])
        plane_normal_vec = Vector(plane_normal)
        plane_normal_vec2 = Vector(plane_normal2)

        # donor HIS
        his_plane_aytoms = ["ND1","CE1","NE2","CD2","CG"]
        mat2 = []
        for atom in his_plane_aytoms:
            if residue_j.has_id(atom):
                mat2.append(residue_j[atom].get_coord().flatten())
        if mat2 != []:
            mat2 = np.array(mat2)
        his_a,his_b,his_c,his_d = plane_svd_fit(mat2,False)
        his_normal = np.array([his_a,his_b,his_c])
        his_normal2 = np.array([-his_a,-his_b,-his_c])
        his_normal = Vector(his_normal)
        his_normal2 = Vector(his_normal2)

        distance = []
        theta1 = []
        theta2 = []
        for i in range(0,len(res_j_atom)):
            distance.append(two_pt_dist(res_i_centroid, residue_j[res_j_atom[i]].get_coord()))    
            d_vec = residue_j[res_j_atom[i]].get_coord().flatten() - res_i_centroid
            d_vec = Vector(d_vec)
            angle_x = d_vec.angle(plane_normal_vec)
            theta1_pt = math.degrees(angle_x)
            if not theta1_pt <= 90:
                angle_y = d_vec.angle(plane_normal_vec2)
                theta1_pt = math.degrees(angle_y)
                if not theta1_pt <= 90:
                    raise ValueError('angle Theta1 does not have a physical value [0,90]')
            theta1.append(theta1_pt)

            angle2_x = d_vec.angle(his_normal)
            theta2_pt = math.degrees(angle2_x)
            if not theta2_pt <= 90:
                angle2_y = d_vec.angle(his_normal2)
                theta2_pt = math.degrees(angle2_y)
                if not theta2_pt <= 90:
                    raise ValueError('Theta2 does not have a physical value [0,90]')
            theta2.append(theta2_pt)
    
    return distance, theta1, theta2

# ! call the pre-analysis function (where you known first residue is aromatic, second is either LYS,ARG)
# ! No histidine for the following find_cation_pi_Hb: 
def find_cation_pi_Hb(residue_i, residue_i_idx, residue_j, residue_j_idx, shortest_thresh, pdb_file, chain_id, list_pairs):
    '''This function gets two residues, checks if there exists different residue atoms  within the threshold, and updates the list '''
 #first check basic distance pair defenitions and if true:
    if is_H_bond(residue_i, residue_j):
        is_H_b = True
    else:
        is_H_b = False

    current_distance, theta1, theta2 = cation_pi_dist(residue_i, residue_j, shortest_thresh)
    if residue_j.resname == 'HIS':
        current_distance = np.array(current_distance)
        if not theta1 is None and not theta2 is None:
            closest_atom_idx = np.argmin(current_distance)
            current_distance = np.min(current_distance)
            theta1 = theta1[closest_atom_idx]
            theta2 = theta2[closest_atom_idx]

    if current_distance <= shortest_thresh:
        #calculate planes angle and centroids distance:
        #found a pair--> assign centers (D)

            list_pairs.append({"pdb_file": pdb_file[-8:],
                                    "chain_id": chain_id,
                                    "res_i_id": residue_i_idx,
                                    "res_j_id": residue_j_idx,
                                    "distance": current_distance,
                                        "Theta1":  theta1,
                                        "Theta2": theta2,
                                        "H-bond": is_H_b})
            
            # cprint(f"A cation-pi interaction found within {current_distance} [A] for {pdb_file} chain {chain_id} {residue_i.resname} {residue_i_idx} {residue_j.resname} {residue_j_idx} with angles {theta1} {theta2}",'red')

    return list_pairs


# ! after calculations: analyzing QM files where the first residue is not necessarily aromatic and you cannot tell the order. 
# ! use find_cation_pi_uknown_idx if cation residue is NOT HIS

def find_cation_pi_uknown_idx(residue_i, residue_i_idx, residue_j, residue_j_idx, shortest_thresh, pdb_file, chain_id, list_pairs):
    '''This function gets two residues, checks if there exists different residue atoms  within the threshold, and updates the list '''
 #first check basic distance pair defenitions and if true:
    cations = ['ARG','LYS'] #not his yet
    aromatics = ['PHE', 'TYR', 'TRP', 'HIS']
    if residue_i.resname in aromatics and residue_j.resname in cations:
        current_distance, theta1, theta2 = cation_pi_dist(residue_i, residue_j, shortest_thresh)
    elif residue_j.resname in aromatics and residue_i.resname in cations:
        current_distance, theta1, theta2 = cation_pi_dist(residue_j, residue_i, shortest_thresh)
        #fliped residues above
    else:
        raise ValueError('missing either a cation or an aromatic residue')
    

    if current_distance <= shortest_thresh:
        #calculate planes angle and centroids distance:
        #found a pair--> assign centers (D)

            list_pairs.append({"pdb_file": pdb_file[-8:],
                                       "chain_id": chain_id,
                                       "res_i_id": residue_i_idx,
                                       "res_j_id": residue_j_idx,
                                       "distance": current_distance,
                                        "Theta1":  theta1,
                                        "Theta2": theta2})
            # cprint(f"A cation-pi interaction found within {current_distance} [A] for {pdb_file} chain {chain_id} {residue_i.resname} {residue_i_idx} {residue_j.resname} {residue_j_idx} with angles {theta1} {theta2}",'red')

            return list_pairs
    return list_pairs


# ! if there is HIS as cation use the following:
def count_H_atoms(residue):
    cnt = 0
    for atom in residue.get_list():
        #print(atom.get_name())
        if 'H' in atom.get_name():
            cnt = cnt + 1
    print(cnt)
    return cnt

def find_cation_pi_uknown_idx_HIS_cation(residue_i, residue_i_idx, residue_j, residue_j_idx, shortest_thresh, pdb_file, chain_id, list_pairs):
    '''This function gets two residues, checks if there exists different residue atoms  within the threshold, and updates the list '''
 #first check basic distance pair defenitions and if true:
    cations = ['ARG','LYS','HIS'] #not his yet
    aromatics = ['PHE', 'TYR', 'TRP', 'HIS']
    if not (residue_i.resname == 'HIS' and residue_j.resname == 'HIS'):
        if residue_i.resname == 'HIS' and residue_j.resname in aromatics:
            current_distance, theta1, theta2 = cation_pi_dist(residue_j, residue_i, shortest_thresh)
        elif residue_j.resname == 'HIS' and residue_i.resname in aromatics:
            current_distance, theta1, theta2 = cation_pi_dist(residue_i, residue_j, shortest_thresh)
        elif residue_i.resname in cations and residue_j.resname in aromatics:
            current_distance, theta1, theta2 = cation_pi_dist(residue_i, residue_j, shortest_thresh)
        elif residue_j.resname in cations and residue_i.resname in aromatics:
            current_distance, theta1, theta2 = cation_pi_dist(residue_j, residue_i, shortest_thresh)
            #fliped residues above
        else:
            raise ValueError('missing either a cation or an aromatic residue')
        
        current_distance = np.array(current_distance)
        closest_atom = np.argmin(current_distance)
        theta1 = np.array(theta1)
        theta2 = np.array(theta2)
        theta1 = theta1[closest_atom]
        theta2 = theta2[closest_atom]
        current_distance = np.min(current_distance)
    # ! carefull next:
    '''if you know H positions use the following where you know which histidine is positive'''
    # if residue_i.resname == 'HIS' and residue_j.resname == 'HIS':
    #     cnt_i = count_H_atoms(residue_i)
    #     cnopt_j = count_H_atoms(residue_j)
    #     res_i_copy = residue_i.copy()
    #     res_j_copy = residue_j.copy()
    #     if cnt_i > cnt_j: 
    #         # then residue i is HSP, j is HSE/D
    #         current_distance, theta1, theta2 = cation_pi_dist(res_j_copy, res_i_copy, shortest_thresh)
    #     elif cnt_j > cnt_i:
    #         current_distance, theta1, theta2 = cation_pi_dist(res_i_copy, res_j_copy, shortest_thresh)
    #     else:
    #         raise ValueError("Both histidines have same amount of H-->impossible to calculate cation-pi")
    '''if you do not know then search for the closest interaction regardless of charge'''

    if residue_i.resname == 'HIS' and residue_j.resname == 'HIS':
        res_i_copy = residue_i.copy()
        res_j_copy = residue_j.copy()
        current_distance_i, theta1_i, theta2_i = cation_pi_dist(res_j_copy, res_i_copy, shortest_thresh)
        current_distance_j, theta1_j, theta2_j = cation_pi_dist(res_i_copy, res_j_copy, shortest_thresh)

        current_distance_i = np.array(current_distance_i)
        current_distance_j = np.array(current_distance_j)

        current_distance = min(min(current_distance_i),min(current_distance_j))
        closest_atom_i_idx = np.argmin(current_distance_i)
        closest_atom_j_idx = np.argmin(current_distance_j)
        if  current_distance_i[closest_atom_i_idx] <  current_distance_j[closest_atom_j_idx]:
            # then say i is HSP:
            theta1 = theta1_i[closest_atom_i_idx]
            theta2 = theta2_i[closest_atom_i_idx]
        else:
            theta1 = theta1_j[closest_atom_j_idx]
            theta2 = theta2_j[closest_atom_j_idx]

        
    if current_distance <= shortest_thresh:
        #calculate planes angle and centroids distance:
        #found a pair--> assign centers (D)

            list_pairs.append({"pdb_file": pdb_file[-8:],
                                       "chain_id": chain_id,
                                       "res_i_id": residue_i_idx,
                                       "res_j_id": residue_j_idx,
                                       "distance": current_distance,
                                        "Theta1":  theta1,
                                        "Theta2": theta2})
            # cprint(f"A cation-pi interaction found within {current_distance} [A] for {pdb_file} chain {chain_id} {residue_i.resname} {residue_i_idx} {residue_j.resname} {residue_j_idx} with angles {theta1} {theta2}",'red')

            return list_pairs
    return list_pairs







############################################# PI-PI ########################################################################################################################

# first create a reference residue: ring plane is in XY plane such that CG (for non TRP) is along +X.
def create_ref_residue(residue_inp, mat, centroid):
    ''' gets input residue with its centroid and list of atoms coordnates to impose in XY plane with specific O-CG to x+ orientation'''
    residue = residue_inp.copy()
    
    #translate matrix such that the origin is the center of residue:
    mat = np.array(mat-centroid)
    #extract ring's plane
    a,b,c,d = plane_svd_fit(mat, False)
    #rotation matrix to land on XY:
    xy_rot = rot_plane_to_XY(a,b,-c)
    #apply the above rotation:
    res = np.matmul(mat,xy_rot) #coordinaes after rotation to xy plane and center at origin
    rotation1 = xy_rot
    translation = centroid
    # print('translation vec is:')
    # print(translation)
    list_ring_atoms = assign_ring_atoms(residue)
    for atom in residue.get_list():
        if not atom.get_name() in list_ring_atoms:
            residue.detach_child(atom.get_id())

    for atom in residue.get_list():
        coord1 = atom.get_coord() - translation
        coord2 = np.matmul(coord1,rotation1)
        #atom.transform(rotation, translation)
        atom.set_coord(coord2)
    
    #find the CG-O angle with x -axis  
    Point2 = Point(res[0,0], res[0,1])
    Point1 = Point(0,0)
    CG_X_angle = get_angle(Point1, Point2)
    #rotate about z axis with -angle
    theta = CG_X_angle * np.pi/ 180
    #rotation matrix to align CG-O bond on x axis direction finally...
    CG_O_rot = get_rotation_matrix([0,0,1], theta)
    final_coord = np.matmul(res,CG_O_rot)
    rotation2 = CG_O_rot
    # print('rotation ma is:')
    # print(rotation)
    for atom in residue.get_list():
        coordi1 = atom.get_coord()
        coordi2 = np.matmul(coordi1,rotation2)
        atom.set_coord(coordi2)
        #visualize
    #a1,b1,c1,d1 = plane_svd_fit(final, False)
    
    return residue, final_coord

# after creating a refernce once, can reuse it assigning the same input PDB single residue according to the name of the residue we iterate over now 
def assign_reference_res(residue_name, pdbs_path):
    '''This function intends to create the best refrence residues structures to superimpose all residues of the same name on!
    please note to have all the following files containing references in the same  working directory'''
    if residue_name == 'HIS' or residue_name == 'HSE' or residue_name == 'HSP':
        his_path = pdbs_path + 'HIS_ref_1VLA_chainB_minus3.pdb'
        ref_struct = lst.get_pdb_structure(his_path)
        chains_dict = lst.create_chain_dict(ref_struct)
    elif residue_name == 'PHE':
        #phe_path = pdbs_path + 'PHE_ref_6IX8_chainA_201.pdb'
        phe_path = pdbs_path +'PHE_ref_4IAB_C_59_conv.pdb'
        ref_struct = lst.get_pdb_structure(phe_path)
        chains_dict = lst.create_chain_dict(ref_struct)
    elif residue_name == 'TYR':
        tyr_path = pdbs_path + 'TYR_ref_4PF8_chainA_150.pdb'
        ref_struct = lst.get_pdb_structure(tyr_path)
        chains_dict = lst.create_chain_dict(ref_struct)
    elif residue_name == 'TRP':
        trp_path = pdbs_path + 'TRP_ref_1LJ8_chainA_264.pdb'
        ref_struct = lst.get_pdb_structure(trp_path)
        chains_dict = lst.create_chain_dict(ref_struct)
    else:
        raise ValueError('The following residue is not aromatic and cannot be operated on by the following functions. Please choose eithe: HIS,PHE,TYR,TRP')

    for chain_id, chain_obj in chains_dict.items():
        for residue in chain_obj.get_list():
            centroid = assign_res_center(residue)
            list_ring_atoms = assign_ring_atoms(residue)
            mat = []
            for atom in residue.get_list():
                if not atom.get_name() in list_ring_atoms:
                    residue.detach_child(atom.get_id())
            for atom in list_ring_atoms:
                if atom == 'CG' and not residue.has_id(atom):
                    raise ValueError("The given refrence residue does not have a CG atom object")
                if atom == 'CD2' and not residue.has_id(atom):
                    raise ValueError("The given refrence residue does not have a CD2 atom object")
                if residue.has_id(atom):
                    mat.append(residue[atom].get_coord().flatten())
            mat = np.array(mat)
            refrence_residue, refrence_coord = create_ref_residue(residue,mat,centroid)

    return refrence_residue, refrence_coord


def calculate_polar_Ttheta(transformed_res):
    """Get the angle of this line with the z axis."""
    #if residue is HIS or PHE    
    O2 = assign_res_center(transformed_res)
    #the disatnce from z axis (hirzontally shortest from (0,0,O2.z))
    dist_from_z = np.sqrt(O2[0]**2+O2[1]**2)
    O2_height = np.absolute(O2[2]) #just Z    
    theta = np.arctan(O2_height/dist_from_z)
    angle = math.degrees(theta)
    #print(angle)    
    return angle  

def calculate_polar_Tphi(transformed_res):
    """Get the angle of this residue's O2 with the x axis."""        
    O2 = assign_res_center(transformed_res)     
    #project the O2 on X,Y -> z=0
    projected_O2 = Point(O2[0], O2[1])
    Point1 = Point(0,0)    
    #need angle with x axis (1,0)        
    Tphi = get_angle(Point1, projected_O2)
    return Tphi


# ! Now call all these calculations wrapped in a single function
def assign_X_Y_phi_geometry_params_Hb(residue_i, residue_i_idx, residue_j, residue_j_idx, shortest_thresh, centroid_dist_thresh, pdb_file, chain_id, list_pairs, refrece_res_i, refrece_res_j):
    '''This function gets two residues, maximal distance between centroids, and a threshold for shortes interatom contact distance, and returns a dataframe with geometric parameters:
    P, angle between planes, Ttheta1, Ttheta2, Tphi1, Tphi2 and if the interaction is H-bonds. Also returns PDB,chains and residues id's '''
 #first check basic distance pair defenitions and if true:
    if is_H_bond(residue_i, residue_j):
        is_H_b = True
    else:
        is_H_b = False

    if residue_j.resname == 'HIS':
        temp_res = residue_i.copy()
        residue_i = residue_j.copy()
        residue_j = temp_res
    if is_aromatic_pair(residue_i, residue_j, shortest_thresh):
        #calculate planes angle and centroids distance:
        #found a pair--> assign centers (D)

        res_i_center = assign_res_center(residue_i)
        res_j_center = assign_res_center(residue_j)
        res_i_center = np.array(res_i_center)
        res_j_center = np.array(res_j_center)
        centers_dist = two_pt_dist(res_i_center,res_j_center)

        if centers_dist <= centroid_dist_thresh:
        # assign angle between planes (P)
            list_ring_i_atoms = assign_ring_atoms(residue_i)
            mat_i = []
            list_ring_j_atoms = assign_ring_atoms(residue_j)
            mat_j = []

            for atom in residue_i.get_list():
                if not atom.get_name() in list_ring_i_atoms:
                    residue_i.detach_child(atom.get_id())

            for atom in residue_j.get_list():
                if not atom.get_name() in list_ring_j_atoms:
                    residue_j.detach_child(atom.get_id())

            for atom in list_ring_i_atoms:
                if atom == 'CG' and not residue_i.has_id(atom):
                    return list_pairs
                if atom == 'CD2' and not residue_i.has_id(atom):
                    return list_pairs
                #all calculation below expect the CG atom to exist and be the first in the list!!!!!!!
                if residue_i.has_id(atom):
                    mat_i.append(residue_i[atom].get_coord().flatten())

            for atom in list_ring_j_atoms:
                if atom == 'CG' and not residue_j.has_id(atom):
                    return list_pairs
                if atom == 'CD2' and not residue_j.has_id(atom):
                    return list_pairs
                if residue_j.has_id(atom):
                    mat_j.append(residue_j[atom].get_coord().flatten())

            mat_i = np.array(mat_i)
            mat_j = np.array(mat_j)
            ai,bi,ci,di = plane_svd_fit(mat_i, False)
            aj,bj,cj,dj = plane_svd_fit(mat_j, False)

            if (ai==0 and bi==0 and ci==0 and di==0) or (aj==0 and bj==0 and cj==0 and dj==0):
                return list_pairs

            res_i_j_planes_cosa = two_plane_angle(ai,bi,ci,aj,bj,cj)
            planes_angle = np.arccos(np.clip(res_i_j_planes_cosa, -1, 1))*180/np.pi
            if planes_angle > 90.0:
                planes_angle = 180 - planes_angle


            res_i_rotran = get_superimpose_mat(refrece_res_i, residue_i)

            res_j_transformed = apply_rotation_translation(res_i_rotran, residue_j)

            #calculate Ttheta1:
            Ttheta1 = calculate_polar_Ttheta(res_j_transformed)
            Tphi1 = calculate_polar_Tphi(res_j_transformed)

            #same process only now residue_j is refrenced to xy plane:
            res_j_rotran = get_superimpose_mat(refrece_res_j, residue_j)
            res_i_transformed = apply_rotation_translation(res_j_rotran, residue_i)

            Ttheta2 = calculate_polar_Ttheta(res_i_transformed)
            Tphi2 = calculate_polar_Tphi(res_i_transformed)

            if not (0<=Ttheta1<=90 and 0<=Ttheta2<=90 and 0<=planes_angle<=90):
                raise ValueError(f"Angle has non-physical value, expected 0-90 deg. Got Ttheta1: {Ttheta1}, Ttheta2: {Ttheta2}, plains_angle: {planes_angle}")
            ## visualize
            # joined_mat = []
            # for row in mat_i:
            #     joined_mat.append(row)
            # for row in mat_j:
            #     joined_mat.append(row)

            # joined_mat = np.array(joined_mat)
            # a,b,c,d = plane_svd_fit(joined_mat, False)
            #now categorize to ineraction geometery starting with the plane angles P:
            cprint(f"pairwise interaction in file {pdb_file[-8:]} chain {chain_id} res {residue_i_idx} and res {residue_j_idx} ", "yellow")
            list_pairs.append({"pdb_file": pdb_file[-8:],
                                       "chain_id": chain_id,
                                       "res_i_id": residue_i_idx,
                                       "res_j_id": residue_j_idx,
                                       "P": planes_angle,
                                       "Ttheta1": Ttheta1,
                                       "Ttheta2": Ttheta2,
                                       "Centroids_dist": centers_dist,
                                       "Tphi1": Tphi1,
                                       "Tphi2": Tphi2,
                                       "H-bond": is_H_b
                                      })

            return list_pairs
        return list_pairs
    return list_pairs









'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CATEGORIZING INTERACTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
##################################################################### H-bonds #####################################################################
def check_bonded_N_H_O_H(residue):
    ''' residue should be donor residue'''
    atoms_list = residue.get_atoms()

    H_list = []
    N_O_list = []

    H_donors = [] # after distance calculation contain only possible H donors
    donors_elec_atoms = [] # their corresponding electronegative atoms
    for atom in atoms_list:
        if atom.get_name().startswith('H'):
            H_list.append(atom)
        elif atom.get_name().startswith('N') or atom.get_name().startswith('O'):
            N_O_list.append(atom)
    
    for h_atom in H_list:
        H_coord = h_atom.get_coord()
        for atom in N_O_list:
            elec_coord = atom.get_coord()
            if two_pt_dist(H_coord, elec_coord) < 1.2:
                H_donors.append(h_atom)
                donors_elec_atoms.append(atom)
            
    return H_donors, donors_elec_atoms

def check_H_b_accpetors_bound_atoms(residue):
    ''' if the residue is non-aromatic then return closest atom that is not H or same atom, if aromatic return None'''
    atoms_list = residue.get_atoms()
    N_O_list = []
    N_O_bound_atoms = []

    for atom in atoms_list:
        if atom.get_name().startswith('N') or atom.get_name().startswith('O'):
            N_O_list.append(atom)
    
    # check if residue is  aromatic:
    aromatic_res = ['HIS','TYR','TRP','PHE']
    if residue.resname in aromatic_res:
        return N_O_list, None

    atoms_list = list(residue.get_atoms())
    dist_mat = np.empty((len(N_O_list),len(atoms_list)))

    for n_o_idx,n_o_atom in enumerate(N_O_list):
        for atm_id, atom in enumerate(atoms_list):
            n_o_coord = n_o_atom.get_coord()
            atom_coord = atom.get_coord()
            dist_mat[n_o_idx][atm_id] = two_pt_dist(n_o_coord, atom_coord)

    # now sort the matrix and note that for any row (N,O atom line) the firxt 0 column includes the 
    # closest atom which is the same atom (N,O appear in both list), So the connected bound atom is in
    # column 1

    for line in dist_mat:
        idx = line.argsort()[1]
        if not atoms_list[idx].get_name().startswith('H'):
            N_O_bound_atoms.append(atoms_list[idx])
        else:
            idx = line.argsort()[2]
            N_O_bound_atoms.append(atoms_list[idx])
    return N_O_list, N_O_bound_atoms


def check_H_b_geometry(donor_res,acceptor_res) -> bool:
    ''' following David Baker's paper: doi:10.1016/S0022-2836(03)00021-4
    Args: 
    - donor_res : coorsponds to a H-b donor
    - acceptor_res : correpsponds to a H-b acceptor
    Returns:
    True for H-b defenition satisfied, else False'''
    
    donor_H, donor_elec_atom = check_bonded_N_H_O_H(donor_res)
    acceptor_elec_atom, acceptor_bound_atom = check_H_b_accpetors_bound_atoms(acceptor_res)

    # print(donor_H)
    # print(donor_elec_atom)
    # print(acceptor_elec_atom)
    # print(acceptor_bound_atom)

    if acceptor_bound_atom == None:
        #print('acceptor res is aromatic')
        list_ring_atoms = assign_ring_atoms(acceptor_res)
        mat = []
        
        acceptor_copy = acceptor_res.copy()
        for atom in acceptor_copy.get_list():
            if not atom.get_name() in list_ring_atoms:
                acceptor_copy.detach_child(atom.get_id())

        for atom in list_ring_atoms:
            if atom == 'CG' and not acceptor_copy.has_id(atom) or atom == 'CD2' and not acceptor_copy.has_id(atom):
                raise ValueError('residue is missing essential atoms')
            #all calculation below expect the CG atom to exist and be the first in the list!!!!!!!
            if acceptor_copy.has_id(atom):
                mat.append(acceptor_copy[atom].get_coord().flatten())
        mat = np.array(mat)        
        a,b,c,d = plane_svd_fit(mat, False)


    ''' ~~~~~~~~~~~~~~~~~ parameters calculations ~~~~~~~~~~~~~~~~~~~~'''
    # for any pair of H donor and any acceptor
    for h_idx,h_atom in enumerate(donor_H):
        for accept_idx,acceptor_atom in enumerate(acceptor_elec_atom):
            h_coord = h_atom.get_coord()
            acceptor_coord = acceptor_atom.get_coord()
            d_h_bound_coord = donor_elec_atom[h_idx].get_coord()
            if not acceptor_bound_atom == None:
                acceptor_bound_coord = acceptor_bound_atom[accept_idx].get_coord()
            # for dihedral cal, define the second plane
            else:
                hb_plane_points = [d_h_bound_coord, h_coord, acceptor_coord]
                hb_plane_points = np.array(hb_plane_points)
                d_h_a_plane_a, d_h_a_plane_b, d_h_a_plane_c, dx = plane_svd_fit(hb_plane_points, False)
            '''~~~~~ accorind to the reference paper'''
            x_dihedral = None # this will be overwritten if the acceptor is aromatic
            delta_HA = two_pt_dist(h_coord, acceptor_coord)
            theta = angle_between_3_points(d_h_bound_coord, h_coord, acceptor_coord)
            if not acceptor_bound_atom == None:
                # acceptor is not aromatic so psi is a 3 point vector
                psi = angle_between_3_points(acceptor_bound_coord, acceptor_coord, h_coord)
            else: # acceptor is aromatic so psi is angle between A-H vector and ring plane
                vec_normal_plane_ang = vec_plane_angle((h_coord-acceptor_coord), a, b, c)
                if not vec_normal_plane_ang == None:
                    psi = 90 + vec_normal_plane_ang
                else:
                    return False 
                # now calculate also dihedral x
                x_cos_alpha = two_plane_angle(a,b,c, d_h_a_plane_a, d_h_a_plane_b, d_h_a_plane_c)
                x_angle = np.arccos(np.clip(x_cos_alpha, -1, 1))
                x_dihedral = x_angle * 180 / np.pi

            if delta_HA <= 2.1 and delta_HA >= 1.4:
                if theta >= 110 and theta <= 180:
                    if psi >= 70 and psi <= 180:
                        return True
            
            elif delta_HA <= 2.8 and delta_HA > 2.1:
                if theta >= 80 and theta <= 180:
                    if psi >= 60 and psi <= 180:
                        return True            
    
    # if delta_HA <= 2.8:
    #     print('deltaHA:', delta_HA)
    #     print('theta:', theta)
    #     print('psi:', psi)
    #     print('X:', x_dihedral)

def is_H_bond(res_i,res_j):
    '''Treating here oly amino acid sidechains'''
    H_b_donors = ['HIS','ARG','ASN','GLN','LYS','SER','THR','TRP','TYR','HSP','HSE']
    H_b_acceptors = ['ASN','ASP','GLN','HIS','SER','THR','TYR','HSP','HSE']

    #check if the residues provided fit the acceptor and donors to begin with
    if not ((res_i.resname in H_b_donors and res_j.resname in H_b_acceptors) or (res_j.resname in H_b_donors and res_i.resname in H_b_acceptors)):
        return False
    if res_i.resname in H_b_donors and res_j.resname in H_b_acceptors:
        if check_H_b_geometry(res_i,res_j):
            return True

    if res_j.resname in H_b_donors and res_i.resname in H_b_acceptors:
        if check_H_b_geometry(res_j,res_i):
            return True
    else:
        return False
    
##################################################################### CH-pi bonds #####################################################################

def check_bonded_C_H(residue):
    ''' find CH donors lists for CH-pi interactions'''
    atoms_list = residue.get_atoms()

    H_list = []
    C_list = []    
    H_donors = [] # after distance calculation contain only possible H donors
    donors_elec_atoms = [] # their corresponding electronegative atoms

    for atom in atoms_list:
        if atom.get_name().startswith('H'):
            H_list.append(atom)
        elif atom.get_name().startswith('C'):
            C_list.append(atom)
    for h_atom in H_list:
        H_coord = h_atom.get_coord()
        for atom in C_list:
            elec_coord = atom.get_coord()
            if two_pt_dist(H_coord, elec_coord) < 1.2:
                H_donors.append(h_atom)
                donors_elec_atoms.append(atom)
            
    return H_donors, donors_elec_atoms

def check_CH_b_geometry(donor_res,acceptor_res) -> bool:
    ''' following Brandl-Weiss system: doi:https://doi.org/10.1006/jmbi.2000.4473
    Args: 
    - donor_res : coorsponds to a H-b donor
    - acceptor_res : correpsponds to a H-b acceptor
    Returns:
    True for H-b defenition satisfied, else False'''
    

    # print(donor_H)
    # print(donor_elec_atom)
    # print(acceptor_elec_atom)
    # print(acceptor_bound_atom)

    list_ring_atoms = assign_ring_atoms(acceptor_res)
    donor_copy = donor_res.copy()    
    acceptor_copy = acceptor_res.copy()

    if list_ring_atoms == []:
        return False
    for atom in list_ring_atoms:
        if not acceptor_copy.has_id(atom):
            return False
    for atom in acceptor_copy.get_list():
        #print(f'{atom}+1')
        if not atom.get_name() in list_ring_atoms:
            acceptor_copy.detach_child(atom.get_id())

    pdbs_path = '/home_d/rivka/Bioinformatics_HIS/python_codes/'
    reference_res, reference_coord = assign_reference_res(acceptor_copy.resname, pdbs_path)
    acceptor_rotran = get_superimpose_mat(reference_res, acceptor_copy)
    donor_transformed = apply_rotation_translation(acceptor_rotran, donor_copy)        
    donor_H, donor_elec_atom = check_bonded_C_H(donor_transformed)


    # now the acceptor is in xy plane with centroid as the origin and the donor_transformed is to calculate from
    ''' ~~~~~~~~~~~~~~~~~ parameters calculations ~~~~~~~~~~~~~~~~~~~~'''
    # for any pair of H donor and any acceptor
    if len(donor_elec_atom) < 1 :
        return False
    
    for c_idx,c_atom in enumerate(donor_elec_atom):
        c_coords = c_atom.get_coord()
        C_X_dist = two_pt_dist(np.array([0,0,0]), c_coords)
        if not C_X_dist < 4.5:
            continue
        h_coord = donor_H[c_idx].get_coord()
        Hp_X_dist = two_pt_dist(np.array([h_coord[0],h_coord[1],0]),np.array([0,0,0]))
        if not Hp_X_dist < 1.2:
            continue
        C_H_X_angle = angle_between_3_points(c_coords, h_coord, np.array([0,0,0]))
        if not C_H_X_angle > 120:
            continue            
        else: # it fits all parameters
            #print(C_X_dist,C_H_X_angle,Hp_X_dist)
            return True   
 

def is_CH_bond(res_i,res_j):
    # ! only aromatic so far
    '''Treating here oly amino acid sidechains'''
    H_b_donors = ['HIS','ARG','PHE','TYR','TRP','HSP','HSE','LYS']
    H_b_acceptors = ['HIS','PHE','TYR','TRP','HSE','HSP']

    #check if the residues provided fit the acceptor and donors to begin with

    if not ((res_i.resname in H_b_donors and res_j.resname in H_b_acceptors) or (res_j.resname in H_b_donors and res_i.resname in H_b_acceptors)):
        return False
    if res_i.resname in H_b_donors and res_j.resname in H_b_acceptors:
        if check_CH_b_geometry(res_i,res_j):
            return True

    if res_j.resname in H_b_donors and res_i.resname in H_b_acceptors:
        if check_CH_b_geometry(res_j,res_i):            
            return True
    else:
        return False
        #print('first is acceptor')  

##################################################################### stacked interactions #####################################################################


def check_stacked_geometry(residue_i,residue_j):
    '''Following https://doi.org/10.1021/acs.jpcb.5b08126'''
    first_ring = residue_i.copy()
    second_ring = residue_j.copy()

    # find first ring normal n1
    list_ring1_atoms = assign_ring_atoms(first_ring)
    mat1 = []

    for atom in list_ring1_atoms:
        if atom == 'CG' and not first_ring.has_id(atom) or atom == 'CD2' and not first_ring.has_id(atom):
            raise ValueError('residue is missing essential atoms')
        #all calculation below expect the CG atom to exist and be the first in the list!!!!!!!
        if first_ring.has_id(atom):
            mat1.append(first_ring[atom].get_coord().flatten())
    mat1 = np.array(mat1)        
    a1,b1,c1,d1 = plane_svd_fit(mat1, False)
    # we donot know to which direction the normal will point
    n1 = np.array([a1,b1,c1])
    n1m = np.array([-a1,-b1,-c1])


    # similarly for the second ring
    list_ring2_atoms = assign_ring_atoms(second_ring)
    mat2 = []

    for atom in list_ring2_atoms:
        if atom == 'CG' and not second_ring.has_id(atom) or atom == 'CD2' and not second_ring.has_id(atom):
            raise ValueError('residue is missing essential atoms')
        #all calculation below expect the CG atom to exist and be the first in the list!!!!!!!
        if second_ring.has_id(atom):
            mat2.append(second_ring[atom].get_coord().flatten())
    mat2 = np.array(mat2)        
    a2,b2,c2,d2 = plane_svd_fit(mat2, False)
    n2 = np.array([a2,b2,c2])
    n2m = np.array([-a2,-b2,-c2])
    centroid1 = assign_res_center(first_ring)
    centroid2 = assign_res_center(second_ring)
    centroid1 = np.array(centroid1)
    centroid2 = np.array(centroid2)
    R = two_pt_dist(centroid1,centroid2)
    if not R < 5:
        return False
    
    # calculate alpha brterrn R and n1:
    n1_vec = Vector(n1)
    n1m_vec = Vector(n1m)
    R_vec = centroid2 - centroid1
    R_vec = Vector(R_vec)
    angle_alpha = R_vec.angle(n1_vec)
    alpha = math.degrees(angle_alpha)
    if not alpha <= 90:
        angle_compl = R_vec.angle(n1m_vec)
        alpha = math.degrees(angle_compl)
        if not alpha <=90:
            raise ValueError('angle alpha does not have a physical value [0,90]')
    if not alpha < 45:
        return False
    
    # calculate beta: the angle between n2 and n1
    n2_vec = Vector(n2)
    angle_beta = n2_vec.angle(n1_vec)
    beta = math.degrees(angle_beta)
    if not beta <90:
        compl_beta = n2_vec.angle(n1m_vec)
        beta = math.degrees(compl_beta)
        if not beta <=90:
             raise ValueError('angle beta does not have a physical value [0,90]')
    if not beta < 45:
        return False   

    return True        

def is_pi_stacked(res_i,res_j):
    # ! only aromatic so far
    '''Treating here oly amino acid sidechains'''
    aromatics = ['HIS','PHE','TYR','TRP','HSE','HSP']

    #check if the residues provided fit the acceptor and donors to begin with

    if not (res_i.resname in aromatics and res_j.resname in aromatics):
        #print('the given pair cannot form side-chain H-bonds')
        return False    

    if check_stacked_geometry(res_i,res_j):
        return True

    if check_stacked_geometry(res_j,res_i):
        return True
    else:
        return False
    





            
'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~BACKGROUND SEARCH~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

def extract_one_coord(residue_i, residue_j,pdbs_path: str, atom_name: str,coords_list: list, is_CH_arr: list, is_Hb_arr: list, is_cationpi_arr: list ):
    ''' Extract the coordinates of input atom name after imposing its partner residue on a refrence.
    Params:
    - chain_dict: should contain a single chain from PDB structure with only two residues pairwise interacting)
    - pdbs_path: directory containing the reference residues
    - ref_res_name: user input for the residue to refer the second residue coordinates with respect to
    - atom_name: name of atom to identify within the pair residue of the reference one
    - coords_list: list containing x,y,z values of all previously analyzed non-reference residue coordinates

    Returns:
    - list_ring_atoms: the names of atoms (in order) of the reference file
    - reference_coord: the coordinates of the reference residue 
    - coords_list: updated with the current coordinates with respect to the input

    '''

    cations = ['ARG','LYS','HIS','HSP'] # not his yet
    aromatics = ['PHE', 'TYR', 'TRP', 'HIS', 'HSE', 'HSP']
    shortest_thresh = 6

    is_cation_pi_marker = False
    if is_H_bond(residue_i, residue_j):
        is_H_b = True
    else:
        is_H_b = False
    if is_CH_bond(residue_i, residue_j):
        is_CH_b = True
        #print('is_CHb',residue_i.resname,residue_j.resname)
    else:
        is_CH_b = False
                    
    is_pipi = False
    # if first residue is HIS then we consider HIS as a cation! (we donot know). Unless there is LYS or ARG there
    # so unless HIS is first, analyze HIS as the acceptor in this analysis
    if (not residue_i.resname in cations) and (residue_i.resname in aromatics and residue_j.resname in aromatics):
        is_pipi = True
    if is_pipi:
        print('is pipi')
        centroid1 = assign_res_center(residue_i)
        centroid2 = assign_res_center(residue_j)
        if centroid1 is None or centroid2 is None:
            return coords_list
        centroid_dist = two_pt_dist(np.array(centroid1),np.array(centroid2))
        if (not centroid_dist < 10):
            return coords_list, is_CH_arr, is_Hb_arr, is_cationpi_arr

        if is_aromatic_pair(residue_i, residue_j, 4.6):
        # here assume first residue is HIS and the second uknown
            reference_res, reference_coord = assign_reference_res(residue_j.resname, pdbs_path)
            ref_res_copy = residue_j.copy()
            partner_res = residue_i.copy()
            list_ring_atoms = assign_ring_atoms(ref_res_copy)

            for atom in ref_res_copy.get_list():
                if not atom.get_name() in list_ring_atoms:
                    ref_res_copy.detach_child(atom.get_id())
            if list_ring_atoms == []:
                return coords_list, is_CH_arr, is_Hb_arr, is_cationpi_arr
            for atom in list_ring_atoms:
                if not ref_res_copy.has_id(atom):
                    return coords_list, is_CH_arr, is_Hb_arr, is_cationpi_arr
                
            res_ref_rotran = get_superimpose_mat(reference_res, ref_res_copy)
            partner_res_transformed = apply_rotation_translation(res_ref_rotran, partner_res) # apply same matrix calculations on the cation
            if partner_res_transformed.has_id(atom_name):
                atom_coords = partner_res_transformed[atom_name].get_coord()
                coords_list.append(atom_coords)
                is_CH_arr.append(is_CH_b)
                is_Hb_arr.append(is_H_b)
                is_cationpi_arr.append(is_cation_pi_marker)
            return coords_list, is_CH_arr, is_Hb_arr, is_cationpi_arr
        
        return coords_list, is_CH_arr, is_Hb_arr, is_cationpi_arr
    
    # if is not pi-pi then check for cation-pi analysis
    #
    if residue_i.resname == 'HIS' and not residue_j.resname == 'HIS':
        current_distance, theta1, theta2 = cation_pi_dist(residue_j, residue_i, shortest_thresh)
        current_distance = np.array(current_distance)
        closest_atom_idx = np.argmin(current_distance)
        current_distance = np.min(current_distance)
    else:
        current_distance, theta1, theta2 = cation_pi_dist(residue_i, residue_j, shortest_thresh)
    

    if residue_i.resname == 'HIS' and residue_j.resname == 'HIS':
        res_i_copy = residue_i.copy()
        res_j_copy = residue_j.copy()
        current_distance_i, theta1_i, theta2_i = cation_pi_dist(res_j_copy, res_i_copy, shortest_thresh)
        if theta1_i==None or theta2_i == None:
            return coords_list, is_CH_arr, is_Hb_arr, is_cationpi_arr
        current_distance_j, theta1_j, theta2_j = cation_pi_dist(res_i_copy, res_j_copy, shortest_thresh)
        if theta1_j==None or theta2_j == None:
            return coords_list, is_CH_arr, is_Hb_arr, is_cationpi_arr
        current_distance_i = np.array(current_distance_i)
        current_distance_j = np.array(current_distance_j)
        if current_distance_i is None or current_distance_j is None:
            return coords_list, is_CH_arr, is_Hb_arr, is_cationpi_arr
        current_distance_i = np.array(current_distance_i)
        current_distance_j = np.array(current_distance_j)
        theta1_i = np.array(theta1_i)
        theta2_i = np.array(theta2_i)
        theta1_j = np.array(theta1_j)
        theta2_j = np.array(theta2_j)
     
        current_distance = min(min(current_distance_i),min(current_distance_j))
        closest_atom_i_idx = np.argmin(current_distance_i)
        closest_atom_j_idx = np.argmin(current_distance_j)

        if current_distance_i[closest_atom_i_idx] < current_distance_j[closest_atom_j_idx]:
            # then say i is HSP:
            theta1 = theta1_i[closest_atom_i_idx]
            theta2 = theta2_i[closest_atom_i_idx]

        elif current_distance_i[closest_atom_i_idx] > current_distance_j[closest_atom_j_idx] and not (current_distance_j[0]==current_distance_j[1]):
            theta2 = theta2_j[closest_atom_j_idx]
        else:
            return coords_list, is_CH_arr, is_Hb_arr, is_cationpi_arr
    if not current_distance <= shortest_thresh:
        return coords_list, is_CH_arr, is_Hb_arr, is_cationpi_arr
    
    is_cation_pi = False
    #check for LYS/ARG first (HIS-LYS)
    if residue_i.resname in aromatics and residue_j.resname in cations:
        reference_res, reference_coord = assign_reference_res(residue_i.resname, pdbs_path)
        ref_res_copy = residue_i.copy()
        ref_res = residue_i.copy() #non modified for later representation
        cation_res = residue_j.copy()
        is_cation_pi = True

    elif residue_j.resname in aromatics and residue_i.resname in cations:
        reference_res, reference_coord = assign_reference_res(residue_j.resname, pdbs_path)
        ref_res_copy = residue_j.copy()
        cation_res = residue_i.copy()
        ref_res = residue_j.copy()
        is_cation_pi = True

    if is_cation_pi:
        residue_i_idx = residue_i.get_id()[1]
        residue_j_idx = residue_j.get_id()[1]
        chain_id = 'A'
        title='RES'
        shortest_thresh = 10
        list_pairs = []
        list_pairs = find_cation_pi_uknown_idx(residue_i, residue_i_idx, residue_j, residue_j_idx, shortest_thresh, pdbs_path, chain_id, list_pairs)
        # ! carefull running HIS:
        #list_pairs = find_cation_pi_uknown_idx_HIS_cation(residue_i, residue_i_idx, residue_j, residue_j_idx, shortest_thresh, pdbs_path, chain_id, list_pairs)
        
        D_val = list_pairs[0]['distance']
        T1_val = list_pairs[0]['Theta1']
        Rx_vals = np.abs(D_val*np.sin(T1_val*np.pi/180))
        if Rx_vals <= 2.3:
            is_cation_pi_marker = True
        list_ring_atoms = assign_ring_atoms(ref_res_copy)

        for atom in ref_res_copy.get_list():
            if not atom.get_name() in list_ring_atoms:
                ref_res_copy.detach_child(atom.get_id())
        
        for atom in list_ring_atoms:
            # all calculation below expect the CG atom to exist and be the first in the list!!!!!!!
            if not ref_res_copy.has_id(atom):
                return coords_list, is_CH_arr, is_Hb_arr, is_cationpi_arr


        res_ref_rotran = get_superimpose_mat(reference_res, ref_res_copy)
        cation_res_transformed = apply_rotation_translation(res_ref_rotran, cation_res) # apply same matrix calculations on the cation
        if cation_res_transformed.has_id(atom_name):
            atom_coords = cation_res_transformed[atom_name].get_coord()
            coords_list.append(atom_coords)
            is_CH_arr.append(is_CH_b)
            is_Hb_arr.append(is_H_b)
            is_cationpi_arr.append(is_cation_pi_marker)
            return coords_list, is_CH_arr, is_Hb_arr, is_cationpi_arr

                    
    return coords_list, is_CH_arr, is_Hb_arr, is_cationpi_arr	



