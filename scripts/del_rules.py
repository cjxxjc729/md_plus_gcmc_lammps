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

dir_4_ftmp="operation_center.tmp_files"


#编辑rules来增删原子

def gen_del_info(atoms, region, ele, n_del_tresh=1e9):

   chem_pot = -7

   f_ae_lpj = 'ae.ann.lammpstrj'
   f_ae_lpj_last = make_f_lpjlast(f_ae_lpj) 

   type_map = get_type_map()

   #ele_lpj = list(type_map).index(ele) + 1

   atm_ids_of_ele = [atom.index for atom in atoms if atom.symbol == ele]
   atm_ids_of_ele = np.array(atm_ids_of_ele)
   print("index of ", ele , "is ",atm_ids_of_ele)
   print("n_chose maximum =", n_del_tresh)

   edit_the_fae(f_ae_lpj_last)
   ae=np.loadtxt(f_ae_lpj_last)

   line_ids = np.where(ae[:,2]>chem_pot)[0]
   atm_ids_beyond_u = ae[line_ids,0]-1
   #print(atoms[atm_ids_beyond_u].get_chemical_symbols())

   del_atm_ids = np.intersect1d(atm_ids_of_ele, atm_ids_beyond_u)  #做交集，既要是ele，又要是ae大于u的

   print("index to del is ", del_atm_ids)

   if len(del_atm_ids) > n_del_tresh: #最大的只能取到n_del_tresh

     print("del_atm_ids is ",del_atm_ids,"larger then n_del_tresh ",n_del_tresh, "randomly pick some")
     del_atm_ids=np.random.choice(del_atm_ids,size=n_del_tresh,replace=False)     


   return del_atm_ids



def edit_the_fae(fae):
 
  lines_after = []  

  with open(fae, 'r') as file:
    # 读取文件内容到列表
    lines = file.readlines()
    # 初始化一个标记，指示是否找到特定行
    found = False
    for line in lines:
        # 如果找到了特定行，设置标记为True
        if "ITEM: ATOMS id type c_ae0" in line:
            found = True
            continue  # 跳过当前行，直接处理下一行
        # 如果已经找到了特定行，则保存后续所有行
        if found:
            lines_after.append(line)

  # 如果需要，将截取的内容写入新文件
  with open(fae, 'w') as file:
      file.writelines(lines_after)

