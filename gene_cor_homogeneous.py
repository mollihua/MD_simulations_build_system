#! /usr/bin/env python
# PURPOSE: to generate a molecular system containing randomly positioned BGLC's around a protein.
# SYNTAX: rancor_atomdist.py prot.cor bglc1.cor prot_bglcs.cor
# Note: prot.cor is the coordinate file of protein solute
#       bglc1.cor is the coordinate of a BGLC unit
#       BGLC refers to beta-glucopyranose


import random
import math
import numpy as np
import sys, string, os
from scipy.spatial import cKDTree


# define the number of BGLC's
bglc_num = 20


# the radius of the smallest sphere that can accommodate a BGLC
bglc_r = 6  # unit: angstrom


# calculate the boxsize that allows positioning of the geometric centers of BGLC's
watboxsize = 40    # dimension of the cubic molecular system (unit: angstrom)
boxsize = watboxsize - 2*(bglc_r + 2)   # with a margin of 2 angstroms


# read in two pdb files
prot = open (sys.argv[1], 'r').readlines()[2:]
bglc_unit = open (sys.argv[2], 'r').readlines()[2:]
filename = sys.argv[3]
try:
    os.remove( filename )
except OSError:
    pass
bglc_sys = open( filename, 'a' )

atomnum_unit = len( bglc_unit )
bglc_sys.write( '* bglc_sys' + '\n' )
bglc_sys.write( str( bglc_num * atomnum_unit ) + '\n' )


# Calculate the geometric center
def get_center ( corfile ):
    xlist, ylist, zlist = [], [], []
    for line in corfile:
        items = string.split( line )
        xlist.append( float( items[4] ) )
        ylist.append( float( items[5] ) )
        zlist.append( float( items[6] ) )
    xarr = np.array( xlist )
    yarr = np.array( ylist )
    zarr = np.array( zlist )
    xave = np.mean( xarr )
    yave = np.mean( yarr )
    zave = np.mean( zarr )
    return ([xave, yave, zave])


# Get the x,y,z coordinates
def get_coor ( corfile ):
    xlist, ylist, zlist = [], [], []
    for line in corfile:
        items = string.split( line )
        xlist.append( float( items[4] ) )
        ylist.append( float( items[5] ) )
        zlist.append( float( items[6] ) )
    xarr = np.array( xlist )
    yarr = np.array( ylist )
    zarr = np.array( zlist )
    return ( np.array([xarr, yarr, zarr]) )


# Get the list of atom names
def get_names ( corfile, colnum ):
    namelist = []
    for line in corfile:
        items = string.split( line )
        namelist.append( str( items[ colnum ] ) )
    name_arr = np.array( namelist )
    return( name_arr )


# Generate a random coordinate centered around (0,0,0) with the given dimension
def gen_rancor ( dim ):
    xran = random.random()* dim - 0.5*dim
    yran = random.random()* dim - 0.5*dim
    zran = random.random()* dim - 0.5*dim
    return ( np.array([xran, yran, zran]) )


# Generate a rotation matrix with random rotation angle along x, y, z axes
def rot(): 
    a = random.random()*180
    b = random.random()*180
    c = random.random()*180
    rot_mat1 = np.matrix( [[1,0,0],[0,np.cos(a),np.sin(a)],[0,-np.sin(a),np.cos(a)]] )
    rot_mat2 = np.matrix( [[np.cos(b),0,-np.sin(b)],[0,1,0],[np.sin(b),0,np.cos(b)]] )
    rot_mat3 = np.matrix( [[np.cos(c),np.sin(c),0],[-np.sin(c),np.cos(c),0],[0,0,1]] )
    return( [rot_mat1, rot_mat2, rot_mat3] )


# Generate BGLC coordinates around the protein
def gen_newsys ( cor1, cor2 ):   # cor1: protein; cor2: BGLC unit
    
    # for rotation: calc relative coordinates of each BGLC atom wrt the BGLC center
    cen_bglc1 = get_center( cor2 )
    xyz_bglc1 = get_coor( cor2 )
    
    ox_bglc1 = xyz_bglc1[0] - cen_bglc1[0]
    ox_bglc1_mat = np.matrix( ox_bglc1 )
    oy_bglc1 = xyz_bglc1[1] - cen_bglc1[1]
    oy_bglc1_mat = np.matrix( oy_bglc1 )
    oz_bglc1 = xyz_bglc1[2] - cen_bglc1[2]
    oz_bglc1_mat = np.matrix( oz_bglc1 )
    o_xyz = np.vstack( (ox_bglc1_mat, oy_bglc1_mat, oz_bglc1_mat) )
    
    # for translation: get prot info
    cen_prot = get_center( cor1 )
    
    # get the list of names
    resname = get_names( cor2, 2 )
    atomname = get_names( cor2, 3 )
    segname = get_names( cor2, 7 )
    
    # start to generate BGLC coordinates
    xyz_prot = get_coor( cor1 )
    xyz_prot = np.array( zip( *xyz_prot ) ) # format: array([[x1,y1,z1],...])
    coor_accum = xyz_prot.copy()
    i = 0
    while i < bglc_num:
        print i
        # rotate bglc randomly
        rot_mat = rot()
        o_xyz_n1 = rot_mat[0]*o_xyz
        o_xyz_n2 = rot_mat[1]*o_xyz_n1
        o_xyz_n3 = rot_mat[2]*o_xyz_n2
        
        # translate bglc randomly
        cen_bglci = gen_rancor( boxsize ) + cen_prot
        cen_bglci_mat = np.matrix([[cen_bglci[0]],[cen_bglci[1]],[cen_bglci[2]]])  # 3X1 matrix
        o_xyz_n4 = o_xyz_n3 + cen_bglci_mat
        o_xyz_n5 = zip( *np.array( o_xyz_n4 ) )
        o_xyz_n5 = np.array( o_xyz_n5 )

        # check overlap using nearest neighbor algorithm (KDTree) with critera: mindist > 4 angstroms
        tree = cKDTree( coor_accum )
        distlist = []
        for atominfo in enumerate( o_xyz_n5 ):
            atom_xyz = atominfo[1]
            dist, index = tree.query( atom_xyz )
            distlist.append( dist )
        mindist = min( distlist )
        if mindist >= 4:
            coor_accum = np.vstack( (coor_accum, o_xyz_n5) )
            
            # write out coordinate
            for term in enumerate( o_xyz_n5 ):
                ind, xyz = term
                atomnum = i*atomnum_unit + ind + 1
                resid = i+1
                resn = resname[ind]
                atomn = atomname[ind]
                segn = segname[ind]
                x = float( xyz[0] )
                y = float( xyz[1] )
                z = float( xyz[2] )
                newline = '%5i'% int( atomnum ) + '%5i'% int(resid) + '%5s'% resn + ' ' + '%-3s'% atomn + '%11.5f'% x + '%10.5f'% y + '%10.5f' % z + '%4s' % segn + '%4i'% int(resid) + '      0.00000' + '\n'
                bglc_sys.write( str( newline ) )
            i = i + 1
        else:
            pass

    bglc_sys.close()


gen_newsys( prot, bglc_unit )

 
