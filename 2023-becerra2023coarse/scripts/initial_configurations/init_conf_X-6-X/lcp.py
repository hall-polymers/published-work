#####################################################
# Filename: kg_lc.py                                #
# Author: Diego Becerra, 2022                       #
#                                                   #
# Characteristics: Initial configuration of         #
# polymer chains with liquid crystals attached      #
#                                                   #
# Updated to April, 2022                            #
#####################################################
#!/usr/bin/env python

"""
# General functions
"""
import sys
import math
import numpy as np
from scipy.stats import uniform
import random
from scipy import stats
from scipy import constants
import seaborn as sns
from pyquaternion import Quaternion

class Status:
    def __init__(self, type, x, y, z, vx, vy, vz, angmomx, angmomy, angmomz,
                 quatw, quatx, quaty, quatz, shapex, shapey, shapez, molecule,
                 q, ellipsoidflag, density):
        self.type = type
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.angmomx = angmomx
        self.angmomy = angmomy
        self.angmomz = angmomz
        self.quatw = quatw
        self.quatx = quatx
        self.quaty = quaty
        self.quatz = quatz
        self.shapex = shapex
        self.shapey = shapey
        self.shapez = shapez
        self.molecule = molecule
        self.q = q
        self.ellipsoidflag = ellipsoidflag
        self.density = density

class Ellipsoid:
    def __init__(self, id, shapex, shapey, shapez, quatw, quatx, quaty, quatz):
        self.id = id
        self.shapex = shapex
        self.shapey = shapey
        self.shapez = shapez
        self.quatw = quatw
        self.quatx = quatx
        self.quaty = quaty
        self.quatz = quatz

class Bond:
    def __init__(self, type, a, b):
        self.type = type
        self.a = a
        self.b = b

class Angle:
    def __init__(self, type, a, b, c):
        self.type = type
        self.a = a
        self.b = b
        self.c = c

class Dihedral:
    def __init__(self, type, a, b, c, d):
        self.type = type
        self.a = a
        self.b = b
        self.c = c
        self.d = d

class LmpConf:
    def __init__(self):
        self.status = []
        self.ellipsoids = []
        self.velocities = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []


    def initialconfiguration(self):
        # Input parameters
        # nch = int(intput('Enter the number of polymer chains'))
        # nbbpc = int(intput('Enter the number of backbone beads per chain'))
        nch = 4                   # number of chains
        nbbpc = 100                # n of backbone beads per chain (even and >2)
        spacer = 6                  # spacer composed of #spacer_unit beads
        nlcpc = int(nbbpc / 2)      # number of liquid crystals per chain
        shape_bb = [1.0, 1.0, 1.0]
        shape_sc = [1.0, 1.0, 1.0]
        shape_lc = [1.0, 1.0, 3.0]
        nsbpc = spacer * nlcpc       # number of spacer beads per chain
        nmpc = nlcpc                # number of monomers per chain
        monomer_unit = spacer + 3   # monomers composed of #monomer_unit beads
        nbpc=int(nbbpc + nsbpc + nlcpc)    # number of beads per chain
        n_beads = int(nbpc * nch)
        n_backbone_beads = int(nbbpc * nch)
        n_spacers = int(nsbpc * nch)
        n_ellipsoids = int((nbbpc + nsbpc + nlcpc) * nch)
        n_bonds = (nbpc - 1) * nch
        n_angles = (
                   ((nbbpc - 2) + (spacer - 1) * nlcpc)
                   + ((nbbpc - 1) + (spacer - 1) * nlcpc)
                   + (nlcpc + 1)
                   + (nlcpc)
                   + (nlcpc)
                   + (nlcpc)
                   + (nlcpc)
                   ) * nch
        # n_angles = ((nbbpc - 2) + (spacer - 1) * nlcpc) * nch

        atomtypes = 3
        bondtypes = 2
        angletypes = 7

        mass_bead = 1.0
        volume_bead = float(4 / 3 * math.pi * (0.5 ** 3))
        mass_spacer = 1.0
        volume_spacer = float(4 / 3 * math.pi * (0.5 ** 3))
        mass_ellipsoid = 1.5
        volume_ellipsoid = float(4 / 3 * math.pi * (0.5 * 1.5 * 0.5))
        density_bead = float(mass_bead / volume_bead)
        density_spacer = float(mass_spacer / volume_spacer)
        density_ellipsoid = float(mass_ellipsoid / volume_ellipsoid)

        lb = 1.00                        # bond distance between beads
        lc_oblate = 2.00                        # bond distance lc-fixer
        lc_prolate = 1.0                        # bond distance lc-fixer
        theta1 = 18 * math.pi / 180    # angle between beads
        Nav = 6.022140857e+23              # Avogadro constant [mol^-1]
        rho = 0.00009                        # LCE density
        Mwc = (n_backbone_beads * mass_bead
               + n_spacers * mass_spacer
               + n_ellipsoids * mass_ellipsoid)
        N = int((nbbpc + nsbpc + nlcpc) * nch)
        TMw = Mwc * nch
        Tm = TMw / Nav
        #V=Tm/rho*1e+21
        V = float(N / rho)
        axis = V ** (float(1) / float(3))

        def rotation_matrix_from_vectors(vec1, vec2):
            """ Find the rotation matrix that aligns vec1 to vec2
            :param vec1: A 3d "source" vector
            :param vec2: A 3d "destination" vector
            :return mat: A transform matrix (3x3) which when applied
            to vec1, aligns it with vec2.
            """

            a, b = ((vec1 / np.linalg.norm(vec1)).reshape(3),
                    (vec2 / np.linalg.norm(vec2)).reshape(3))
            v = np.cross(a, b)
            if any(v):    # if not all zeros then
                c = np.dot(a, b)
                s = np.linalg.norm(v)
                kmat = np.array([[0, -v[2], v[1]], [v[2], 0,
                                  -v[0]], [-v[1], v[0], 0]])
                return np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
            else:
                return np.eye(3)


        orientation_lc = (int(input('Select: oblate[1] or prolate[2]?')))
        mol = 0
        for molec in range(0, nch):
            cont = 1
            mol += 1
            # First bead (with random position in the space inside the box)
            fbx = np.random.uniform(axis * 0.2, axis * 0.8)
            fby = np.random.uniform(axis * 0.2, axis * 0.8)
            fbz = np.random.uniform(axis * 0.2, axis * 0.8)
            first_bead = (fbx, fby, fbz)
            # Random vector in the space (from the origin)
            rv1x = np.random.uniform(-axis, axis)
            rv1y = np.random.uniform(-axis, axis)
            rv1z = np.random.uniform(-axis, axis)
            rv1 = (rv1x, rv1y, rv1z)
            # Normalized vector between rv1 and the first bead
            dx = np.subtract(rv1x, fbx)
            dy = np.subtract(rv1y, fby)
            dz = np.subtract(rv1z, fbz)
            norm = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
            rotated_vector1 = (dx / norm, dy / norm, dz / norm)
            # Second bead
            sbx = fbx + rotated_vector1[0] * lb
            sby = fby + rotated_vector1[1] * lb
            sbz = fbz + rotated_vector1[2] * lb
            second_bead = (sbx, sby, sbz)
            # obtaining the backbone quaternion
            uibbx = ((sbx - fbx))
            uibby = ((sby - fby))
            uibbz = ((sbz - fbz))
            uibb_norm = math.sqrt(uibbx ** 2
                                  + uibby ** 2
                                  + uibbz ** 2
                                  )
            uibb = (
                uibbx / uibb_norm,
                uibby / uibb_norm,
                uibbz / uibb_norm,
                )
            u0 = [0, 0, 1]
            u_bb = [uibb[0], uibb[1], uibb[2]]
            R = rotation_matrix_from_vectors(u0, u_bb)
            my_quaternionsbb = Quaternion(matrix=R)
            qbbw = my_quaternionsbb[0]
            qbbx = my_quaternionsbb[1]
            qbby = my_quaternionsbb[2]
            qbbz = my_quaternionsbb[3]
            self.status.append(Status(int(1),
                                      float(fbx), float(fby), float(fbz),
                                      float(0.0), float(0.0), float(0.0),
                                      float(0.0), float(0.0), float(0.0),
                                      float(qbbw),
                                      float(qbbx), float(qbby), float(qbbz),
                                      float(shape_bb[0]),
                                      float(shape_bb[1]),
                                      float(shape_bb[2]),
                                      int(mol), float(0.0), int(1),
                                      float(density_bead),
                                      )
                                      )
            while len(self.status) < nbpc * mol:
                #(step 1)
                rv2x = np.random.uniform(-axis, axis)
                rv2y = np.random.uniform(-axis, axis)
                rv2z = np.random.uniform(-axis, axis)
                rv2 = (rv2x, rv2y, rv2z)
                norm = math.sqrt(rv2x ** 2 + rv2y ** 2 + rv2z ** 2)
                rv2_norm = (rv2x / norm, rv2y / norm, rv2z / norm)
                # Cross Product (step 2)
                crossx = (
                    rotated_vector1[1] * rv2_norm[2]
                    - rotated_vector1[2] * rv2_norm[1]
                    )
                crossy = -(
                    rotated_vector1[0] * rv2_norm[2]
                    - rotated_vector1[2] * rv2_norm[0]
                    )
                crossz = (
                    rotated_vector1[0] * rv2_norm[1]
                    - rotated_vector1[1] * rv2_norm[0]
                    )
                crossvector = (crossx, crossy, crossz)
                # Rotated vector (step 3)
                rotated_vector2 = (
                    Quaternion(axis=crossvector, angle=theta1).
                    rotate(rotated_vector1)
                    )
                # third bead (step 4)
                tbx = sbx + rotated_vector2[0] * lb
                tby = sby + rotated_vector2[1] * lb
                tbz = sbz + rotated_vector2[2] * lb
                third_bead = (tbx, tby, tbz)
                # plane fb-sb-tb
                sfx = np.subtract(sbx, fbx)
                sfy = np.subtract(sby, fby)
                sfz = np.subtract(sbz, fbz)
                tsx = np.subtract(tbx, sbx)
                tsy = np.subtract(tby, sby)
                tsz = np.subtract(tbz, sbz)
                # spacer direction
                spacerx = sfy * tsz - sfz * tsy
                spacery = -(sfx * tsz - sfz * tsx)
                spacerz = sfx * tsy - sfy * tsx
                spacer_norm = math.sqrt(spacerx ** 2
                                        + spacery ** 2
                                        + spacerz ** 2)
                spacervector = (
                    spacerx / spacer_norm,
                    spacery / spacer_norm,
                    spacerz / spacer_norm,
                    )
                # spacer beads, syndiotactic 1
                spacer1_syn1x = sbx + spacervector[0] * lb
                spacer1_syn1y = sby + spacervector[1] * lb
                spacer1_syn1z = sbz + spacervector[2] * lb
                spacer2_syn1x = spacer1_syn1x + spacervector[0] * lb
                spacer2_syn1y = spacer1_syn1y + spacervector[1] * lb
                spacer2_syn1z = spacer1_syn1z + spacervector[2] * lb
                spacer3_syn1x = spacer2_syn1x + spacervector[0] * lb
                spacer3_syn1y = spacer2_syn1y + spacervector[1] * lb
                spacer3_syn1z = spacer2_syn1z + spacervector[2] * lb
                spacer4_syn1x = spacer3_syn1x + spacervector[0] * lb
                spacer4_syn1y = spacer3_syn1y + spacervector[1] * lb
                spacer4_syn1z = spacer3_syn1z + spacervector[2] * lb
                spacer5_syn1x = spacer4_syn1x + spacervector[0] * lb
                spacer5_syn1y = spacer4_syn1y + spacervector[1] * lb
                spacer5_syn1z = spacer4_syn1z + spacervector[2] * lb
                spacerl_syn1x = spacer5_syn1x + spacervector[0] * lb
                spacerl_syn1y = spacer5_syn1y + spacervector[1] * lb
                spacerl_syn1z = spacer5_syn1z + spacervector[2] * lb
                # spacer beads, syndiotactic 2
                spacer1_syn2x = sbx - spacervector[0] * lb
                spacer1_syn2y = sby - spacervector[1] * lb
                spacer1_syn2z = sbz - spacervector[2] * lb
                spacer2_syn2x = spacer1_syn2x - spacervector[0] * lb
                spacer2_syn2y = spacer1_syn2y - spacervector[1] * lb
                spacer2_syn2z = spacer1_syn2z - spacervector[2] * lb
                spacer3_syn2x = spacer2_syn2x - spacervector[0] * lb
                spacer3_syn2y = spacer2_syn2y - spacervector[1] * lb
                spacer3_syn2z = spacer2_syn2z - spacervector[2] * lb
                spacer4_syn2x = spacer3_syn2x - spacervector[0] * lb
                spacer4_syn2y = spacer3_syn2y - spacervector[1] * lb
                spacer4_syn2z = spacer3_syn2z - spacervector[2] * lb
                spacer5_syn2x = spacer4_syn2x - spacervector[0] * lb
                spacer5_syn2y = spacer4_syn2y - spacervector[1] * lb
                spacer5_syn2z = spacer4_syn2z - spacervector[2] * lb
                spacerl_syn2x = spacer5_syn2x - spacervector[0] * lb
                spacerl_syn2y = spacer5_syn2y - spacervector[1] * lb
                spacerl_syn2z = spacer5_syn2z - spacervector[2] * lb
                # Prolate normalized vector
                lc_oblx = spacervector[0]
                lc_obly = spacervector[1]
                lc_oblz = spacervector[2]
                oblate = (
                    spacervector[0],
                    spacervector[1],
                    spacervector[2],
                    )
                # Oblate normalized vector
                lc_prox = ((tbx - fbx))
                lc_proy = ((tby - fby))
                lc_proz = ((tbz - fbz))
                lc_pro_norm = math.sqrt(lc_prox ** 2
                                        + lc_proy ** 2
                                        + lc_proz ** 2
                                        )
                lc_pro = (
                    lc_prox / lc_pro_norm,
                    lc_proy / lc_pro_norm,
                    lc_proz / lc_pro_norm,
                    )
                prolate = (
                    lc_pro[0],
                    lc_pro[1],
                    lc_pro[2],
                    )
                if orientation_lc == 1:
                    lc_norm = oblate
                    lc = lc_oblate
                    lc_syn1x = spacerl_syn1x + spacervector[0] * 2 * lc
                    lc_syn1y = spacerl_syn1y + spacervector[1] * 2 * lc
                    lc_syn1z = spacerl_syn1z + spacervector[2] * 2 * lc
                    lc_syn1 = (lc_syn1x, lc_syn1y, lc_syn1z)
                    lc_syn2x = spacerl_syn2x - spacervector[0] * 2 * lc
                    lc_syn2y = spacerl_syn2y - spacervector[1] * 2 * lc
                    lc_syn2z = spacerl_syn2z - spacervector[2] * 2 * lc
                    lc_syn2 = (lc_syn2x, lc_syn2y, lc_syn2z)
                else:
                    lc_norm = prolate
                    lc = lc_prolate
                    lc_syn1x = spacerl_syn1x + spacervector[0] * lc
                    lc_syn1y = spacerl_syn1y + spacervector[1] * lc
                    lc_syn1z = spacerl_syn1z + spacervector[2] * lc
                    lc_syn1 = (lc_syn1x, lc_syn1y, lc_syn1z)
                    lc_syn2x = spacerl_syn2x - spacervector[0] * lc
                    lc_syn2y = spacerl_syn2y - spacervector[1] * lc
                    lc_syn2z = spacerl_syn2z - spacervector[2] * lc
                    lc_syn2 = (lc_syn2x, lc_syn2y, lc_syn2z)
                # Obtaining the backbone quaternion from (fi, vi, ui)
                uibbx = ((tbx - sbx))
                uibby = ((tby - sby))
                uibbz = ((tbz - sbz))
                uibb_norm = math.sqrt(uibbx ** 2
                                      + uibby ** 2
                                      + uibbz ** 2
                                      )
                uibb = (
                    uibbx / uibb_norm,
                    uibby / uibb_norm,
                    uibbz / uibb_norm,
                    )
                u0_bb = [0, 0, 1]
                u_bb = [uibb[0], uibb[1], uibb[2]]
                R1 = rotation_matrix_from_vectors(u0_bb, u_bb)
                f0_bb = [1, 0, 0]
                f_bb = Quaternion(matrix=R1).rotate(f0_bb)
                R2 = rotation_matrix_from_vectors(f_bb, spacervector)
                R2R1 = np.dot(R2, R1)
                my_quaternionsbb = Quaternion(matrix=R2R1)
                qbbw = my_quaternionsbb[0]
                qbbx = my_quaternionsbb[1]
                qbby = my_quaternionsbb[2]
                qbbz = my_quaternionsbb[3]
                # Obtaining the side chain quaternion from (fi, vi, ui)
                # syndiotactic 1
                u0_sc = [0, 0, 1]
                u_sc = [spacervector[0], spacervector[1], spacervector[2]]
                R1 = rotation_matrix_from_vectors(u0_sc, u_sc)
                f0_sc = [1, 0, 0]
                f_sc = Quaternion(matrix=R1).rotate(f0_sc)
                R2 = rotation_matrix_from_vectors(f_sc, lc_pro)
                R2R1 = np.dot(R2, R1)
                my_quaternionssc = Quaternion(matrix=R2R1)
                qscw_s1 = my_quaternionssc[0]
                qscx_s1 = my_quaternionssc[1]
                qscy_s1 = my_quaternionssc[2]
                qscz_s1 = my_quaternionssc[3]
                # syndiotactic 2
                u0_sc = [0, 0, 1]
                u_sc = [- spacervector[0], - spacervector[1], -spacervector[2]]
                R1 = rotation_matrix_from_vectors(u0_sc, u_sc)
                f0_sc = [1, 0, 0]
                f_sc = Quaternion(matrix=R1).rotate(f0_sc)
                neg_lc_pro = (
                    - lc_prox / lc_pro_norm,
                    - lc_proy / lc_pro_norm,
                    - lc_proz / lc_pro_norm,
                    )
                R2 = rotation_matrix_from_vectors(f_sc, neg_lc_pro)
                R2R1 = np.dot(R2, R1)
                my_quaternionssc = Quaternion(matrix=R2R1)
                qscw_s2 = my_quaternionssc[0]
                qscx_s2 = my_quaternionssc[1]
                qscy_s2 = my_quaternionssc[2]
                qscz_s2 = my_quaternionssc[3]
                # Obtaining the liquid crystal quaternion from (fi, vi, ui)
                # syndiotactic 1
                u0_lc = [0, 0, 1]
                u_lc = [lc_norm[0], lc_norm[1], lc_norm[2]]
                R1 = rotation_matrix_from_vectors(u0_lc, u_lc)
                if orientation_lc == 1:
                    v0_lc = [0, 1, 0]
                    v_lc = Quaternion(matrix=R1).rotate(v0_lc)
                    R2 = rotation_matrix_from_vectors(v_lc, lc_pro)
                else:
                    f0_lc = [1, 0, 0]
                    f_lc = Quaternion(matrix=R1).rotate(f0_lc)
                    R2 = rotation_matrix_from_vectors(f_lc, spacervector)
                R2R1 = np.dot(R2, R1)
                my_quaternionslc = Quaternion(matrix=R2R1)
                qlcw_s1 = my_quaternionslc[0]
                qlcx_s1 = my_quaternionslc[1]
                qlcy_s1 = my_quaternionslc[2]
                qlcz_s1 = my_quaternionslc[3]
                # syndiotactic 2
                u0_lc = [0, 0, 1]
                u_lc = [- lc_norm[0], - lc_norm[1], - lc_norm[2]]
                R1 = rotation_matrix_from_vectors(u0_lc, u_lc)
                if orientation_lc == 1:
                    v0_lc = [0, 1, 0]
                    v_lc = Quaternion(matrix=R1).rotate(v0_lc)
                    neg_lc_pro = (
                        - lc_prox / lc_pro_norm,
                        - lc_proy / lc_pro_norm,
                        - lc_proz / lc_pro_norm,
                        )
                    R2 = rotation_matrix_from_vectors(v_lc, neg_lc_pro)
                else:
                    f0_lc = [1, 0, 0]
                    f_lc = Quaternion(matrix=R1).rotate(f0_lc)
                    neg_spacervector = (
                        - spacerx / spacer_norm,
                        - spacery / spacer_norm,
                        - spacerz / spacer_norm,
                        )
                    R2 = rotation_matrix_from_vectors(f_lc, neg_spacervector)
                R2R1 = np.dot(R2, R1)
                my_quaternionslc = Quaternion(matrix=R2R1)
                qlcw_s2 = my_quaternionslc[0]
                qlcx_s2 = my_quaternionslc[1]
                qlcy_s2 = my_quaternionslc[2]
                qlcz_s2 = my_quaternionslc[3]

                # from scipy.spatial.transform import Rotation as Ro
                # r = Ro.from_quat([my_quaternionslc[1], my_quaternionslc[2], my_quaternionslc[3], my_quaternionslc[0]]) #obtaining rot matrix from quaternion
                # a = Quaternion(matrix=r.as_matrix()).rotate(u0_lc)
                #
                # print(u_lc,Quaternion(matrix=R1).rotate(u0_lc))
                # print('')
                # print(a)

                #Second way to obtain quaternions
                # R = my_quaternionslc.rotation_matrix    #obtaining rotation matrix from quaternions
                # ua = my_quaternionslc.axis #obtaining rotation axis from quaternions
                # thetaa = my_quaternionslc.angle #obtaining rotation angle from quaternions
                # my_quaternion2 = Quaternion(axis=ua, angle=thetaa) #obtiannig quaternion from axis and angle
                # my_vector2=Quaternion(axis=ua,angle=thetaa).rotate(f_lc) #obtaing rotated vector from reference vector with axis and angle
                # from scipy.spatial.transform import Rotation as Ro
                # r = Ro.from_quat([my_quaternionslc[1], my_quaternionslc[2], my_quaternionslc[3], my_quaternionslc[0]]) #obtaining rot matrix from quaternion
                # print(R,r.as_matrix())

                if 1 == 1:#0.1<tbx<(axis-0.1) and 0.1<tby<(axis-0.1) and 0.1<tbz<(axis-0.1) and 0.1<lcx_syn1<(axis-0.1) and 0.1<lcy_syn1<(axis-0.1) and 0.1<lcz_syn1<(axis-0.1):
                    cont+=1
                    print(cont)
                    rotated_vector1 = rotated_vector2
                    fbx = sbx
                    fby = sby
                    fbz = sbz
                    sbx = tbx
                    sby = tby
                    sbz = tbz
                    if (cont % 2) != 0:
                        self.status.append(Status(int(1),
                                           float(fbx), float(fby), float(fbz),
                                           float(0.0), float(0.0), float(0.0),
                                           float(0.0), float(0.0), float(0.0),
                                           float(qbbw),
                                           float(qbbx),
                                           float(qbby),
                                           float(qbbz),
                                           float(shape_bb[0]),
                                           float(shape_bb[1]),
                                           float(shape_bb[2]),
                                           int(mol), float(0.0), int(1),
                                           float(density_bead))
                                           )
                    if (cont % 2) == 0 and (cont % 4) != 0:
                        self.status.append(Status(int(1),
                                    float(fbx), float(fby), float(fbz),
                                    float(0.0), float(0.0), float(0.0),
                                    float(0.0), float(0.0), float(0.0),
                                    float(qbbw),
                                    float(qbbx), float(qbby), float(qbbz),
                                    float(shape_bb[0]),
                                    float(shape_bb[1]),
                                    float(shape_bb[2]),
                                    int(mol), float(0.0), int(1),
                                    float(density_bead)
                                    )
                                    )
                        self.status.append(Status(int(2),
                                    float(spacer1_syn1x), float(spacer1_syn1y),
                                    float(spacer1_syn1z),
                                    float(0.0), float(0.0), float(0.0),
                                    float(0.0), float(0.0), float(0.0),
                                    float(qscw_s1),
                                    float(qscx_s1),
                                    float(qscy_s1),
                                    float(qscz_s1),
                                    float(shape_sc[0]),
                                    float(shape_sc[1]),
                                    float(shape_sc[2]),
                                    int(mol), float(0.0), int(1),
                                    float(density_spacer),
                                    )
                                    )
                        self.status.append(Status(int(2),
                                    float(spacer2_syn1x), float(spacer2_syn1y),
                                    float(spacer2_syn1z),
                                    float(0.0), float(0.0), float(0.0),
                                    float(0.0), float(0.0), float(0.0),
                                    float(qscw_s1),
                                    float(qscx_s1),
                                    float(qscy_s1),
                                    float(qscz_s1),
                                    float(shape_sc[0]),
                                    float(shape_sc[1]),
                                    float(shape_sc[2]),
                                    int(mol), float(0.0), int(1),
                                    float(density_spacer),
                                    )
                                    )
                        self.status.append(Status(int(2),
                                    float(spacer3_syn1x), float(spacer3_syn1y),
                                    float(spacer3_syn1z),
                                    float(0.0), float(0.0), float(0.0),
                                    float(0.0), float(0.0), float(0.0),
                                    float(qscw_s1),
                                    float(qscx_s1),
                                    float(qscy_s1),
                                    float(qscz_s1),
                                    float(shape_sc[0]),
                                    float(shape_sc[1]),
                                    float(shape_sc[2]),
                                    int(mol), float(0.0), int(1),
                                    float(density_spacer),
                                    )
                                    )
                        self.status.append(Status(int(2),
                                    float(spacer4_syn1x), float(spacer4_syn1y),
                                    float(spacer4_syn1z),
                                    float(0.0), float(0.0), float(0.0),
                                    float(0.0), float(0.0), float(0.0),
                                    float(qscw_s1),
                                    float(qscx_s1),
                                    float(qscy_s1),
                                    float(qscz_s1),
                                    float(shape_sc[0]),
                                    float(shape_sc[1]),
                                    float(shape_sc[2]),
                                    int(mol), float(0.0), int(1),
                                    float(density_spacer),
                                    )
                                    )
                        self.status.append(Status(int(2),
                                    float(spacer5_syn1x), float(spacer5_syn1y),
                                    float(spacer5_syn1z),
                                    float(0.0), float(0.0), float(0.0),
                                    float(0.0), float(0.0), float(0.0),
                                    float(qscw_s1),
                                    float(qscx_s1),
                                    float(qscy_s1),
                                    float(qscz_s1),
                                    float(shape_sc[0]),
                                    float(shape_sc[1]),
                                    float(shape_sc[2]),
                                    int(mol), float(0.0), int(1),
                                    float(density_spacer),
                                    )
                                    )
                        self.status.append(Status(int(2),
                                    float(spacerl_syn1x), float(spacerl_syn1y),
                                    float(spacerl_syn1z),
                                    float(0.0), float(0.0), float(0.0),
                                    float(0.0), float(0.0), float(0.0),
                                    float(qscw_s1),
                                    float(qscx_s1),
                                    float(qscy_s1),
                                    float(qscz_s1),
                                    float(shape_sc[0]),
                                    float(shape_sc[1]),
                                    float(shape_sc[2]),
                                    int(mol), float(0.0), int(1),
                                    float(density_spacer),
                                    )
                                    )
                        self.status.append(Status(int(3),
                                    float(lc_syn1x), float(lc_syn1y),
                                    float(lc_syn1z),
                                    float(0.0), float(0.0), float(0.0),
                                    float(0.0), float(0.0), float(0.0),
                                    float(qlcw_s1),
                                    float(qlcx_s1),
                                    float(qlcy_s1),
                                    float(qlcz_s1),
                                    float(shape_lc[0]),
                                    float(shape_lc[1]),
                                    float(shape_lc[2]),
                                    int(mol), float(0.0), int(1),
                                    float(density_ellipsoid),
                                    )
                                    )
                    if (cont % 2) == 0 and (cont % 4) == 0:
                        self.status.append(Status(int(1),
                                    float(fbx), float(fby), float(fbz),
                                    float(0.0), float(0.0), float(0.0),
                                    float(0.0), float(0.0), float(0.0),
                                    float(qbbw),
                                    float(qbbx), float(qbby), float(qbbz),
                                    float(shape_bb[0]),
                                    float(shape_bb[1]),
                                    float(shape_bb[2]),
                                    int(mol), float(0.0), int(1),
                                    float(density_bead)
                                    )
                                    )
                        self.status.append(Status(int(2),
                                    float(spacer1_syn2x), float(spacer1_syn2y),
                                    float(spacer1_syn2z),
                                    float(0.0), float(0.0), float(0.0),
                                    float(0.0), float(0.0), float(0.0),
                                    float(qscw_s2),
                                    float(qscx_s2),
                                    float(qscy_s2),
                                    float(qscz_s2),
                                    float(shape_sc[0]),
                                    float(shape_sc[1]),
                                    float(shape_sc[2]),
                                    int(mol), float(0.0), int(1),
                                    float(density_spacer)
                                    )
                                    )
                        self.status.append(Status(int(2),
                                    float(spacer2_syn2x), float(spacer2_syn2y),
                                    float(spacer2_syn2z),
                                    float(0.0), float(0.0), float(0.0),
                                    float(0.0), float(0.0), float(0.0),
                                    float(qscw_s2),
                                    float(qscx_s2),
                                    float(qscy_s2),
                                    float(qscz_s2),
                                    float(shape_sc[0]),
                                    float(shape_sc[1]),
                                    float(shape_sc[2]),
                                    int(mol), float(0.0), int(1),
                                    float(density_spacer)
                                    )
                                    )
                        self.status.append(Status(int(2),
                                    float(spacer3_syn2x), float(spacer3_syn2y),
                                    float(spacer3_syn2z),
                                    float(0.0), float(0.0), float(0.0),
                                    float(0.0), float(0.0), float(0.0),
                                    float(qscw_s2),
                                    float(qscx_s2),
                                    float(qscy_s2),
                                    float(qscz_s2),
                                    float(shape_sc[0]),
                                    float(shape_sc[1]),
                                    float(shape_sc[2]),
                                    int(mol), float(0.0), int(1),
                                    float(density_spacer)
                                    )
                                    )
                        self.status.append(Status(int(2),
                                    float(spacer4_syn2x), float(spacer4_syn2y),
                                    float(spacer4_syn2z),
                                    float(0.0), float(0.0), float(0.0),
                                    float(0.0), float(0.0), float(0.0),
                                    float(qscw_s2),
                                    float(qscx_s2),
                                    float(qscy_s2),
                                    float(qscz_s2),
                                    float(shape_sc[0]),
                                    float(shape_sc[1]),
                                    float(shape_sc[2]),
                                    int(mol), float(0.0), int(1),
                                    float(density_spacer)
                                    )
                                    )
                        self.status.append(Status(int(2),
                                    float(spacer5_syn2x), float(spacer5_syn2y),
                                    float(spacer5_syn2z),
                                    float(0.0), float(0.0), float(0.0),
                                    float(0.0), float(0.0), float(0.0),
                                    float(qscw_s2),
                                    float(qscx_s2),
                                    float(qscy_s2),
                                    float(qscz_s2),
                                    float(shape_sc[0]),
                                    float(shape_sc[1]),
                                    float(shape_sc[2]),
                                    int(mol), float(0.0), int(1),
                                    float(density_spacer)
                                    )
                                    )
                        self.status.append(Status(int(2),
                                    float(spacerl_syn2x), float(spacerl_syn2y),
                                    float(spacerl_syn2z),
                                    float(0.0), float(0.0), float(0.0),
                                    float(0.0), float(0.0), float(0.0),
                                    float(qscw_s2),
                                    float(qscx_s2),
                                    float(qscy_s2),
                                    float(qscz_s2),
                                    float(shape_sc[0]),
                                    float(shape_sc[1]),
                                    float(shape_sc[2]),
                                    int(mol), float(0.0), int(1),
                                    float(density_spacer)
                                    )
                                    )
                        self.status.append(Status(int(3),
                                    float(lc_syn2x), float(lc_syn2y),
                                    float(lc_syn2z),
                                    float(0.0), float(0.0), float(0.0),
                                    float(0.0), float(0.0), float(0.0),
                                    float(qlcw_s2),
                                    float(qlcx_s2),
                                    float(qlcy_s2),
                                    float(qlcz_s2),
                                    float(shape_lc[0]),
                                    float(shape_lc[1]),
                                    float(shape_lc[2]),
                                    int(mol), float(0.0), int(1),
                                    float(density_ellipsoid)
                                    )
                                    )

        file = open('traj.lammpstrj','w')
        file.write("ITEM: TIMESTEP\n")
        file.write("%s\n" % int(0))
        file.write("ITEM: NUMBER OF ATOMS\n")
        file.write("%s\n" % (len(self.status)))
        file.write("ITEM: BOX BOUNDS pp pp pp\n")
        file.write("%s\t%s\n" % (float(0), axis))
        file.write("%s\t%s\n" % (float(0), axis))
        file.write("%s\t%s\n" % (float(0), axis))
        file.write("ITEM: ATOMS id type x y z vx vy vz angmomx angmomy angmomz  c_orient[1] c_orient[2] c_orient[3] c_orient[4] c_shape[1] c_shape[2] c_shape[3]\n")
        for j in range(len(self.status)):
            s = self.status[j]
            file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
            (j+1, s.type,
            s.x, s.y, s.z,
            s.vx, s.vy, s.vz,
            s.angmomx, s.angmomy, s.angmomz,
            s.quatx, s.quaty, s.quatz, s.quatw,
            s.shapex, s.shapey, s.shapez,)
            )
        file.close()

        file = open('in.lammps', 'w')
        file.write("LAMMPS configuration file modified using lmp_io\n\n")
        file.write("\t%ld atoms\n" % n_beads)
        file.write("\t%ld ellipsoids\n" % n_ellipsoids)
        file.write("\t%ld bonds\n" % n_bonds)
        file.write("\t%ld angles\n\n" % n_angles)
        file.write("\t%ld atom types\n" % atomtypes)
        file.write("\t%ld bond types\n" % bondtypes)
        file.write("\t%ld angle types\n\n" % angletypes)
        file.write("\t%lf\t%lf xlo xhi\n" % (0, axis))
        file.write("\t%lf\t%lf ylo yhi\n" % (0, axis))
        file.write("\t%lf\t%lf zlo zhi\n\n" % (0, axis))
        file.write("Masses\n\n")
        file.write("\t%ld\t%lf\n" % (1, mass_bead))
        file.write("\t%ld\t%lf\n" % (2, mass_spacer))
        file.write("\t%ld\t%lf\n" % (3, mass_ellipsoid))
        file.write("\n")
        file.write("Atoms\n\n")
        for k in range(len(self.status)):
            s = self.status[k]
            file.write("\t%ld\t%ld\t%lf\t%lf\t%lf\t%ld\t%lf\t%ld\t%lf\t%d\t%d\t%d\n" %
            (k+1, s.type,
            s.x, s.y, s.z,
            s.molecule, s.q, s.ellipsoidflag, s.density, 0, 0, 0)
            )
            if s.ellipsoidflag == 1:
                self.ellipsoids.append(
                    Ellipsoid(k+1, s.shapex, s.shapey, s.shapez,
                    s.quatw, s.quatx, s.quaty, s.quatz)
                    )
        file.write("\n")
        file.write("Ellipsoids\n\n")
        for k in range(len(self.ellipsoids)):
            e = self.ellipsoids[k]
            file.write("\t%ld\t%lf\t%lf\t%f\t%lf\t%lf\t%f\t%lf\n" %
            (e.id, e.shapex, e.shapey, e.shapez,
            e.quatw, e.quatx, e.quaty, e.quatz)
            )
        file.write("\n")
        # file.write("Velocities\n\n")
        # for k in range(len(self.status)):
        #     v = self.status[k]
        #     file.write("\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n" %
        #     (k+1, v.vx, v.vy, v.vz,
        #     v.angmomx, v.angmomy, v.angmomz)
        #     )
        # file.write("\n")
        file.write("Bonds\n\n")
        for i in range(0, nch):
            for j in range(0, nmpc):
                self.bonds.append(Bond(
                    int(1),
                    int((j * monomer_unit + 1) + (i * nbpc)),
                    int((j * monomer_unit + 2) + (i * nbpc)),
                    )
                    )
                for k in range(0, spacer):
                    self.bonds.append(Bond(
                        int(1),
                        int((k + 2) + (j * monomer_unit) + (i * nbpc)),
                        int((k + 3) + (j * monomer_unit) + (i * nbpc)),
                        )
                        )
                self.bonds.append(Bond(
                    int(2),
                    int((j * monomer_unit + 8) + (i * nbpc)),
                    int((j * monomer_unit + 9) + (i * nbpc)),
                    )
                    )
                if j == 0 or j % (nmpc-1) != 0:
                    self.bonds.append(Bond(
                        int(1),
                        int((j * monomer_unit + 2) + (i * nbpc)),
                        int((j * monomer_unit + 10) + (i * nbpc)),
                        )
                        )
        for k in range(len(self.bonds)):
            b = self.bonds[k]
            file.write("\t%ld\t%ld\t%ld\t%d\n" % (k+1, b.type, b.a, b.b))
        file.write("\n")
        file.write("Angles\n\n")
        for i in range(0, nch):
            for j in range(0, nmpc):
                if j == 0 or j % (nmpc-1) != 0:
                    if j == 0:
                        self.angles.append(Angle(
                            int(2),
                            int((j * monomer_unit + 1) + (i * nbpc)),
                            int((j * monomer_unit + 2) + (i * nbpc)),
                            int((j * monomer_unit + 2) + (i * nbpc)),
                            )
                            )
                    self.angles.append(Angle(
                        int(1),
                        int((j * monomer_unit + 1) + (i * nbpc)),
                        int((j * monomer_unit + 2) + (i * nbpc)),
                        int((j * monomer_unit + 10) + (i * nbpc)),
                        )
                        )
                    self.angles.append(Angle(
                        int(2),
                        int((j * monomer_unit + 2) + (i * nbpc)),
                        int((j * monomer_unit + 10) + (i * nbpc)),
                        int((j * monomer_unit + 10) + (i * nbpc)),
                        )
                        )
                    self.angles.append(Angle(
                        int(1),
                        int((j * monomer_unit + 2) + (i * nbpc)),
                        int((j * monomer_unit + 10) + (i * nbpc)),
                        int((j * monomer_unit + 11) + (i * nbpc)),
                        )
                        )
                    self.angles.append(Angle(
                        int(2),
                        int((j * monomer_unit + 10) + (i * nbpc)),
                        int((j * monomer_unit + 11) + (i * nbpc)),
                        int((j * monomer_unit + 11) + (i * nbpc)),
                        )
                        )
                if j == (nmpc-1):
                    self.angles.append(Angle(
                        int(3),
                        int((j * monomer_unit + 1) + (i * nbpc)),
                        int((j * monomer_unit + 2) + (i * nbpc)),
                        int((j * monomer_unit + 2) + (i * nbpc)),
                        )
                        )
                for k in range(0, spacer-1):
                    if k == 0:
                        self.angles.append(Angle(
                            int(4),
                            int((k + 2) + (j * monomer_unit) + (i * nbpc)),
                            int((k + 3) + (j * monomer_unit) + (i * nbpc)),
                            int((k + 3) + (j * monomer_unit) + (i * nbpc)),
                            )
                            )
                        self.angles.append(Angle(
                            int(5),
                            int((k + 2) + (j * monomer_unit) + (i * nbpc)),
                            int((k + 3) + (j * monomer_unit) + (i * nbpc)),
                            int((k + 3) + (j * monomer_unit) + (i * nbpc)),
                            )
                            )
                    self.angles.append(Angle(
                        int(1),
                        int((k + 2) + (j * monomer_unit) + (i * nbpc)),
                        int((k + 3) + (j * monomer_unit) + (i * nbpc)),
                        int((k + 4) + (j * monomer_unit) + (i * nbpc)),
                        )
                        )
                    self.angles.append(Angle(
                        int(2),
                        int((k + 3) + (j * monomer_unit) + (i * nbpc)),
                        int((k + 4) + (j * monomer_unit) + (i * nbpc)),
                        int((k + 4) + (j * monomer_unit) + (i * nbpc)),
                        )
                        )
                    if k == (spacer-2):
                        self.angles.append(Angle(
                            int(3),
                            int((k + 3) + (j * monomer_unit) + (i * nbpc)),
                            int((k + 4) + (j * monomer_unit) + (i * nbpc)),
                            int((k + 4) + (j * monomer_unit) + (i * nbpc)),
                            )
                            )
                        self.angles.append(Angle(
                            int(6),
                            int((k + 4) + (j * monomer_unit) + (i * nbpc)),
                            int((k + 5) + (j * monomer_unit) + (i * nbpc)),
                            int((k + 5) + (j * monomer_unit) + (i * nbpc)),
                            )
                            )
                        self.angles.append(Angle(
                            int(7),
                            int((k + 4) + (j * monomer_unit) + (i * nbpc)),
                            int((k + 5) + (j * monomer_unit) + (i * nbpc)),
                            int((k + 5) + (j * monomer_unit) + (i * nbpc)),
                            )
                            )
        for k in range(len(self.angles)):
            d = self.angles[k]
            file.write("\t%ld\t%ld\t%ld\t%d\t%d\n" % (k+1, d.type,
                                                      d.a, d.b, d.c))
        file.write("\n")
        file.close()

def get_command_line_args(args):
    if len(args) != 1:
        print("Usage: %s <input files> <output files>"  % args[0])
        sys.exit(1)
    return args

def main():
    args = get_command_line_args(sys.argv)
    lmpconf = LmpConf()
    lmpconf.initialconfiguration()

if __name__ == "__main__":
    main()
