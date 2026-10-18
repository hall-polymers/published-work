import numpy as np

N_polymer = 500 # number of the polymer chain    

mass=[0, 1.0, 0.15, 1.65, 1.45] #bead mass from type 1 to type n
LA = 100 # chain length of A block   
LB = 0 # chain length of B block
L_polymer = LA + LB

n_anion=0
n_cation=n_anion

Branchinterval_homo = 0 

bcpbackboneAty = 1 # bead type of A block   
bcpbackboneBty = 1 # ... of B block

cationty = 2
anionAty = 3
anionBty = 4

bcpbackboneAbty = 1 #backbone  type of A block   
bcpbackboneBbty = 1
anionbty = 2

bcpbackboneBaty = 1 #angle  type of A block
bcpbackboneAaty = 1

bondb = 1.0 # bond length of block A
bondc = 1.0 
density = 0.88 # if the bond length is 1.0, the default value should be 0.88; if it is 2.0, density should be 0.12

N_total =   + N_polymer * L_polymer + n_anion*3

volume = N_total / density
overall = N_total+1

L = volume**(1.0/3.0)

resid = 0
tmpcount=int(L_polymer)+1
cbcp = np.zeros((N_polymer+1,tmpcount, 3))

alength = np.zeros(N_polymer)
index = 0

def randomwalk(bond):
    theta = np.random.rand() * 2 * np.pi
    dz = 2 * np.random.rand() - 1
    dx = np.sqrt(1 - dz**2) * np.cos(theta)
    dy = np.sqrt(1 - dz**2) * np.sin(theta)
    r = np.sqrt(dx**2 + dy**2 + dz**2)
    scale = bond / r
    dx *= scale
    dy *= scale
    dz *= scale
    return dx, dy, dz



for i in range(1,N_polymer+1):             #constract backbone of BCP
    ind_branch = 0
    for j in range(1,LB+1):                #construct B block 
        dx, dy, dz = randomwalk(bondb)
        tmpind = LA + j
        if j == 1:
            cbcp[i][tmpind, :] = [0, 0, 0]
        if j > 1:
            cbcp[i][tmpind, :] = cbcp=[i][tmpind - 1, :] + [dx, dy, dz]

    for j in range(1,LA+1):                #construct A block
        dx, dy, dz = randomwalk(bondb)
        tmpind = LA - (j - 1)
        if j == 1:
            cbcp[i][tmpind, :] = [dx, dy, dz]
        if j > 1:
            cbcp[i][tmpind, :] = cbcp[i][tmpind + 1, :] + [dx, dy, dz]
            if LB == 0:
                tmpindex = 0
            else:
                tmpindex = 1
               


index = 0
# asign the position of polymer
for j in range(1,N_polymer+1):
    index += 1
    zplane = np.random.rand() * L
    xplane = np.random.rand() * L
    yplane = np.random.rand() * L
    tmp=cbcp[j]     
    cbcp[j][:,0]+=[i+xplane for i in np.zeros(len(tmp[:,1]))]
    cbcp[j][:,1]+=[i+yplane for i in np.zeros(len(tmp[:,1]))]
    cbcp[j][:,2]+=[i+zplane for i in np.zeros(len(tmp[:,1]))]
    #print(cbcp[j])




index = 0
overall=int(overall)
atomid = np.zeros(overall)
chainid = np.zeros(overall)
charge = np.zeros(overall)
typeid = np.zeros(overall)
position = np.zeros((overall, 3))
imagenumber = np.zeros((overall, 3),dtype = int)


for i in range(1,N_polymer+1):
    cbcptmp = cbcp[i]
    for j in range(1,L_polymer+1):
        index = index + 1
        atomid[index] = index
        chainid[index] = i
        charge[index] = 0
        if j <= LA:
            typeid[index] = bcpbackboneAty
        else:
            typeid[index] = bcpbackboneBty
        position[index, :] = cbcptmp[j, :]
        imagenumber[index, :] = [0, 0, 0]

indexa = index

#anionposition = anionposition[anionposition[:, 0] != 0]
for i in range(1,n_cation+1):
    index = index + 1
    atomid[index] = index
    chainid[index] = N_polymer  + i
    charge[index] = 1
    typeid[index] = cationty
    position[index, :] = [np.random.rand() * L,np.random.rand() * L,np.random.rand() * L]
    imagenumber[index, :] = [0, 0, 0]

for i in range(1,n_anion+1):
    index = index + 1
    atomid[index] = index
    chainid[index] = N_polymer + n_cation + i
    charge[index] = -0.87
    typeid[index] = anionAty
    position[index, :] = [np.random.rand() * L,np.random.rand() * L,np.random.rand() * L]
    imagenumber[index, :] = [0, 0, 0]    

for i in range(1,n_anion+1):
    [dx,dy,dz]=randomwalk(bondc)
    index = index + 1
    atomid[index] = index
    chainid[index] = N_polymer + n_cation + i
    charge[index] = -0.13
    typeid[index] = anionBty
    position[index, :] = position[index-n_anion,:]+ [dx, dy, dz] 
    imagenumber[index, :] = [0, 0, 0]


bi = 0
est = (L_polymer - 1) * N_polymer  +  n_anion + 1

A1 = np.zeros(est)
A2 = np.zeros(est)
B1 = np.zeros(est)
Btype = np.zeros(est)
check = np.zeros(est)
check2 = np.zeros(est)

for i in range(1,N_polymer+1):
    for j in range(1,1 + L_polymer - 1):
        bi = bi + 1
        if j == 1:
            A1[bi] = 1 + (i - 1) * L_polymer
            A2[bi] = A1[bi] + 1
        if j > 1:
            A1[bi] = A2[bi-1]
            A2[bi] = A1[bi] + 1
        check[bi] =     (position[int(A1[bi]), 0] - position[int(A2[bi]), 0])**2 + \
                        (position[int(A1[bi]), 1] - position[int(A2[bi]), 1])**2 + \
                        (position[int(A1[bi]), 2] - position[int(A2[bi]), 2])**2

        B1[bi] = bi
        if j > LA:
            Btype[bi] = bcpbackboneBbty
        else:
            Btype[bi] = bcpbackboneAbty


for i in range(1,n_anion):
    bi=bi+1
    A1[bi]=N_polymer*L_polymer+i+n_cation
    A2[bi]=N_polymer*L_polymer+i+n_cation+n_anion
    check[bi] =     (position[int(A1[bi]), 0] - position[int(A2[bi]), 0])**2 + \
                    (position[int(A1[bi]), 1] - position[int(A2[bi]), 1])**2 + \
                    (position[int(A1[bi]), 2] - position[int(A2[bi]), 2])**2
    B1[bi]=bi
    Btype[bi]=anionbty

for i in range(1,bi+1):
    #(check[i])
    if  abs(check[i] - bondb**2) > 0.1 and abs(check[i] - bondc**2) > 0.1:
        
        print('Warning: bond between %i and %i is wrong, length is: %.3f '% (int(A1[i]),int(A2[i]),check[i]))




ai = 0
estangle = overall
An1 = np.zeros(estangle)
An2 = np.zeros(estangle)
An3 = np.zeros(estangle)
ai1 = np.zeros(estangle)
at = np.zeros(estangle)

for i in range(1,N_polymer+1):
    for j in range(1,1 + L_polymer - 2):
        ai = ai + 1
        if j == 1:
            An1[ai] = 1 + (i - 1) * L_polymer
            An2[ai] = An1[ai] + 1
            An3[ai] = An2[ai] + 1
        if j > 1:
            An1[ai] = An2[ai-1]
            An2[ai] = An1[ai] + 1
            An3[ai] = An2[ai] + 1
        if j > LA:
            at[ai] = bcpbackboneBaty
        else:
            at[ai] = bcpbackboneAaty
        ai1[ai] = ai


cp_bcp_branch = ai
checkangle1=np.zeros(overall)
checkangle2=np.zeros(overall)
for i in range(1,1+ai):
    checkangle1[i] = (position[int(An1[i]), 0] - position[int(An2[i]), 0])**2 + \
                     (position[int(An1[i]), 1] - position[int(An2[i]), 1])**2 + \
                     (position[int(An1[i]), 2] - position[int(An2[i]), 2])**2
    checkangle2[i] = (position[int(An2[i]), 0] - position[int(An3[i]), 0])**2 + \
                     (position[int(An2[i]), 1] - position[int(An3[i]), 1])**2 + \
                     (position[int(An2[i]), 2] - position[int(An3[i]), 2])**2
    if  abs(checkangle1[i] - bondb**2) > 0.1 and abs(checkangle1[i] - bondc**2) > 0.1:
        print('Warning: angle between %i and %i is wrong, length is: %.3f '% (int(An1[i]),int(An2[i]),checkangle1[i]))



with open('bulk.lammps', 'w') as file:
    file.write('\n')
    file.write('{:>8}     atoms\n'.format(N_total))
    file.write('{:>8}     bonds\n'.format(bi))
    file.write('{:>8}     angles\n'.format(ai))
    file.write('{:>8}     dihedrals\n'.format(0))
    file.write('\n')

    file.write('{:>8}     atom types\n'.format(int(max(typeid))))
    file.write('{:>8}     bond types\n'.format(int(max(Btype))))
    file.write('{:>8}     angle types\n'.format(int(max(at))))
    file.write('{:>8}     dihedral types\n'.format(0))
    file.write('\n')

    file.write('{:>16.3f} {:>8.3f} xlo xhi\n'.format(0, float(L)))    
    file.write('{:>16.3f} {:>8.3f} ylo yhi\n'.format(0, float(L)))
    file.write('{:>16.3f} {:>8.3f} zlo zhi\n'.format(0, float(L)))
    file.write('\n')
    file.write('Atoms\n')
    file.write('\n')

    for i in range(1,overall):
        file.write('{:>8} {:>8} {:>8} {:>8.2f} {:>8.3f} {:>8.3f} {:>8.3f} {:>8d} {:>8d} {:>8d}\n'.format(int(atomid[i]),
                   int(chainid[i]), int(typeid[i]), charge[i], position[i, 0], position[i, 1], position[i, 2], imagenumber[i, 0],
                   imagenumber[i, 1], imagenumber[i, 2]))

    file.write('\n')
    file.write('Bonds\n')
    file.write('\n')

    for i in range(1,bi+1):
        file.write('{:>8} {:>8} {:>8} {:>8}\n'.format(int(B1[i]), int(Btype[i]), int(A1[i]), int(A2[i])))

    file.write('\n')
    file.write('Angles\n')
    file.write('\n')

    for i in range(1,ai+1):
        file.write('{:>8} {:>8} {:>8} {:>8} {:>8}\n'.format(int(ai1[i]), int(at[i]), int(An1[i]), int(An2[i]), int(An3[i])))
    
        
    file.write('\n')
    file.write('Masses\n')
    file.write('\n')
    for i in range(1,int(max(typeid))+1):
        file.write('     %i  %.2f\n'% (i,float(mass[i])))
