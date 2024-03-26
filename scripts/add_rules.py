#!/data/apps/anaconda/2020.11/bin/python

import sys
sys.path.append("/data/home/scv8373/script/mimic_suite")
sys.path.append("/data/home/scv8373/script/ase_based_constraint_opt_suite")

from ase.io import write
from ase.io import read
from ase import neighborlist,Atoms
from ase import geometry
import numpy as np
from mimic_functions import *
from project_to_grid_functions import *
import os,time,re,random

from adddel_conditions_2_adddel_files import *


def gen_add_info(atoms,region,ele,natm_to_add=1):
  '''
  在区域内随机加natm_to_add个原子  
  '''
  type_map = get_type_map() 

  ele_lpj = list(type_map).index(ele) + 1  

  elelpj_cors = [] 

  for i in range(natm_to_add):

    x = np.random.uniform(region[0,0],region[0,1])
    y = np.random.uniform(region[1,0],region[1,1])
    z = np.random.uniform(region[2,0],region[2,1])
   
    new_cor = [x,y,z]
 
    #优化一下new_cor，让他不要离原有的原子太近
    new_cor = optimize_new_cor(atoms,new_cor)

    elelpj_cors = np.append(elelpj_cors,[ele_lpj,new_cor[0],new_cor[1],new_cor[2]])
 
  elelpj_cors = elelpj_cors.reshape(-1,4) 
 
  return elelpj_cors


def optimize_new_cor(atoms,new_cor):
   
  ds = 0.2
  min_dis_tresh = 0.8
 
  cors=atoms.get_positions()
  cell = atoms.get_cell()
  
  D, Dlen = geometry.get_distances(cors,new_cor,cell=cell,pbc=True)
  D = D[:,0,:]  

  #print("shape of D = ",np.shape(D))
  #print("shape of Dlen = ", np.shape(Dlen))
  
  vec_list=D/Dlen  
 
  #print("shape of vec_list = ", np.shape(vec_list)) 
  #print("vec_list = ",vec_list)

  char_vec = np.sum(vec_list,axis=0)
  char_vec = char_vec / np.linalg.norm(char_vec)
  #print("char_vec = ",char_vec)

  min_dis = np.min(Dlen)
  #print("min_dis =", min_dis)
  #print("new_cor = ",new_cor)

  count = 1 

  while min_dis < min_dis_tresh and count < 30:

    #print("###move iteration time,", count,"#########")
    count = count + 1 
    
    new_cor = new_cor + ds * char_vec
    #print("new_cor = ",new_cor)  
  
    D, Dlen = geometry.get_distances(cors,new_cor,cell=cell,pbc=True)
    D = D[:,0,:]

    vec_list=D/Dlen

    #print("shape of vec_list = ", np.shape(vec_list))
    #print("vec_list = ",vec_list)

    char_vec = np.sum(vec_list,axis=0)
    char_vec = char_vec / np.linalg.norm(char_vec)
    #print("char_vec = ",char_vec)

    min_dis = np.min(Dlen)

    #print("min_dis =", min_dis)

  return new_cor



