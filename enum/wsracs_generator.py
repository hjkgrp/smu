import sys, time, os,shutil,datetime, math
import numpy as np
import matplotlib.pylab as plt
from molSimplify.Classes import globalvars
from molSimplify.Classes import mol3D
from molSimplify.Informatics.autocorrelation import*
from molSimplify.Informatics.misc_descriptors import*
from molSimplify.Informatics.RACassemble import*
from molSimplify.Informatics.graph_analyze import*
from molSimplifyAD.ga_tools import* 
import pickle
from scipy import sparse
import time


# from IPython.core.display import display, HTML
# display(HTML("<style>.container { width:80% !important; }</style>"))


def delete_row_csr(mat, i):
    if not isinstance(mat, sparse.csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")
    n = mat.indptr[i+1] - mat.indptr[i]
    if n > 0:
        mat.data[mat.indptr[i]:-n] = mat.data[mat.indptr[i+1]:]
        mat.data = mat.data[:-n]
        mat.indices[mat.indptr[i]:-n] = mat.indices[mat.indptr[i+1]:]
        mat.indices = mat.indices[:-n]
    mat.indptr[i:-1] = mat.indptr[i+1:]
    mat.indptr[i:] -= n
    mat.indptr = mat.indptr[:-1]
    mat._shape = (mat._shape[0]-1, mat._shape[1])
    
    return mat

class quickClass:
    def __init__(self,name,descriptor_names,descriptors):
        self.name = name
        self.descriptor_names = descriptor_names
        self.descriptors = descriptors

def write_descriptor_csv_local(list_of_runs,name):
    with open(name +'.csv','w') as f:
        f.write('runs,')
        n_cols = len(list_of_runs[0].descriptor_names)
        
        for i, names in enumerate(list_of_runs[0].descriptor_names):
            if i < (n_cols - 1):
                f.write(names + ',')
            else:
                f.write(names + '\n')
        for runs in list_of_runs:
            f.write(runs.name)
            
            for properties in runs.descriptors:
                    f.write(','+str(properties))
            f.write('\n')
            
def SMILEs_to_liglist(smilesstr,conatoms):
    this_mol = mol3D()
    this_mol.getOBMol(smilesstr,'smistring')
    this_mol.convert2mol3D()
    this_lig  = ligand(mol3D(), [],len(conatoms))
    this_lig.mol = this_mol
    return(this_lig)

liglistM = list()
connectingAtomM = list()
with open("../enum/finalSmiMonodentate.txt",'r') as f:
    for line in f:
        liglistM.append(line.split()[0].replace('#4','#'))
        connectingAtomM.append(0)

liglistB = list()
connectingAtomB = list()
with open("../enum/finalSmiBidentate.txt",'r') as f:
    for line in f:
        liglistB.append(line.split()[0].replace('#4','#'))
        connectingAtomB.append(int(line.split()[2])-1)

connectingAtom = connectingAtomM + connectingAtomB
liglist = liglistM + liglistB

SMILEs_strings_M = liglistM
SMILEs_strings_B = liglistB
SMILEs_strings = liglist

ligsM = [SMILEs_to_liglist(s,[0]) for s in SMILEs_strings_M] 
ligsB = [SMILEs_to_liglist(s,[0]) for s in SMILEs_strings_B] 
ligs = [SMILEs_to_liglist(s,[0]) for s in SMILEs_strings] 


####################################

# # generate racs for wsMono
print('Monodentates:')
start_time = time.time()
list_of_runs = list()
block = int(sys.argv[1])
if block == -1:
    pass
else:
    print('Block Number: ' + str(block))
    wsProj = pickle.load( open( "../enum/weaksymMonodentatesProj/weaksymMonodentatesProj_" + str('{:03}'.format(block)) + ".p2", "r" ) )
    wsProj = delete_row_csr(wsProj, 0)
    # print(wsProj.getrow(1))
    print('Shape: ' + str(wsProj.get_shape()))

    for metal in ["Cr", "Mn", "Fe", "Co"]: 
        print('Starting with ' + metal)
        for i in range(0, wsProj.shape[0]): # wsProj.shape[0]
        	
            metal_mol = mol3D()
            metal_mol.addAtom(atom3D(metal)) 
            
            if int(wsProj[i].data[0]) == 4:
                eq = wsProj[i].indices[0]
                top = wsProj[i].indices[1]
                bot = wsProj[i].indices[2]
                  
            elif int(wsProj[i].data[1]) == 4:
                eq = wsProj[i].indices[1]
                top = wsProj[i].indices[0]
                bot = wsProj[i].indices[2]
            
            elif int(wsProj[i].data[2]) == 4:
                eq = wsProj[i].indices[2]
                top = wsProj[i].indices[1]
                bot = wsProj[i].indices[0]

            name = "_".join([metal, 'top', str(top),'bottom', str(bot), 'eq', str(eq)])
            eq_ligs = 4*[ligsM[eq]]
            ax_ligs = [ligsM[top], ligsM[bot]]
                
            eq_cons = 4*[[0]]
            ax_cons = 2*[[0]]
            
            custom_ligand_dict = {"eq_ligand_list":eq_ligs,
                                  "ax_ligand_list":ax_ligs,
                                  "eq_con_int_list":eq_cons,
                                  "ax_con_int_list":ax_cons}

            this_complex = assemble_connectivity_from_parts(metal_mol,custom_ligand_dict)
            con_mat  = this_complex.graph  
            descriptor_names, descriptors = get_descriptor_vector(this_complex,custom_ligand_dict)
            list_of_runs.append(quickClass(name,descriptor_names,descriptors))
            new_time = (time.time() - start_time)
            
    #         plt.figure(figsize=(5,5))
    #         plt.spy(con_mat,precision=0.01, markersize=8) 
            
    #         plt.xticks(np.arange(0,20))
    #         plt.yticks(np.arange(0,20))
    #         plt.grid()
    #         plt.title('smiles complex. eq:' + liglistM[eq] + ' Top:' + liglistM[top] + 'Bot: ' + liglistM[bot])
    #         plt.show()

    print("WS monodentates took: " + str(new_time) + ' s.')
    write_descriptor_csv_local(list_of_runs,"RAC_ws")

# ----

# Bidentates:
print('Bidentates:')
block = int(sys.argv[2])
if block == -1:
    pass
else:
    print('Block Number: ' + str(block))
    start_time = time.time()
    ws_b_Proj = pickle.load( open( "../enum/weaksymBidentatesProj/weaksymBidentatesProj_" + str('{:03}'.format(block)) + ".p2", "r" ) )
    ws_b_Proj = ws_b_Proj.tocsr()
    ws_b_Proj = delete_row_csr(ws_b_Proj, 0)
    # print(ws_b_Proj.getrow(0))
    print('Shape' + str(ws_b_Proj.get_shape()))
    list_of_runs2 =list()


    for metal in ["Cr", "Mn", "Fe", "Co"]:  
        print('Starting with ' + metal)
        for i in range(0, ws_b_Proj.shape[0]):  #ws_b_Proj.shape[0]
            metal_mol = mol3D()
            metal_mol.addAtom(atom3D(metal)) 
            
            if int(ws_b_Proj[i].indices[0]) > 404: # 
                eq = ws_b_Proj[i].indices[0]
                top = ws_b_Proj[i].indices[1]
                bot = ws_b_Proj[i].indices[2]

            elif int(ws_b_Proj[i].indices[1]) > 404: # 
                eq = ws_b_Proj[i].indices[1]
                top = ws_b_Proj[i].indices[0]
                bot = ws_b_Proj[i].indices[2]

            elif int(ws_b_Proj[i].indices[2]) > 404: # 
                eq = ws_b_Proj[i].indices[2]
                top = ws_b_Proj[i].indices[1]
                bot = ws_b_Proj[i].indices[0]
            
            name = "_".join([metal, 'top', str(top),'bottom', str(bot), 'eq', str(eq)])
            eq_ligs = 2*[ligs[eq]]
            ax_ligs = [ligs[top],ligs[bot]]
            
            eq_cons = 2*[[0, connectingAtom[eq]]]
            ax_cons = 2*[[0]]

            custom_ligand_dict = {"eq_ligand_list":eq_ligs,
                                  "ax_ligand_list":ax_ligs,
                                  "eq_con_int_list":eq_cons,
                                  "ax_con_int_list":ax_cons}
            
            this_complex = assemble_connectivity_from_parts(metal_mol, custom_ligand_dict)
            con_mat  = this_complex.graph  
            descriptor_names, descriptors = get_descriptor_vector(this_complex, custom_ligand_dict)
            list_of_runs2.append(quickClass(name, descriptor_names, descriptors))
            new_time = (time.time() - start_time)
            
    # #         plt.figure(figsize=(5,5))
    # #         plt.spy(con_mat,precision=0.01, markersize=8) 
            
    # #         plt.xticks(np.arange(0,20))
    # #         plt.yticks(np.arange(0,20))
    # #         plt.grid()
    # #         plt.title('smiles complex. eq:' + liglistB[eq] + ' Ax:' + liglistB[ax])
    # #         plt.show()


    print("WS Bidentates took: " + str(new_time) + ' s.')
    write_descriptor_csv_local(list_of_runs2,"RAC_ws_b")