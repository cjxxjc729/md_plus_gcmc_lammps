#!/data/run01/scv8373/tools/deepmd-kit-2.2.6/bin/python3.10


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
import os,time,re,random,shutil

import mendeleev



ref_element=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si']

dir_delinfo="operation_center.del_info"
dir_addinfo="operation_center.add_info"
dir_4_ftmp="operation_center.tmp_files"



#-------------------------general rules -----------------------------------------------------#


def initialize():

  #initilize natm_to_add and natm_to_del. set them to be zero.
  #即使没有要增或删除的原子，也得写个0，避免报错

  np.savetxt('./'+dir_addinfo+'/natm_to_add.txt',[0],fmt='%d')
  np.savetxt('./'+dir_delinfo+'/natm_to_del.txt',[0],fmt='%d')


def mkdir(path):
  if os.path.exists(path):
    shutil.rmtree(path)
  os.makedirs(path, exist_ok=True)


def make_f_lpjlast(f_trj):
  
  '''
  将ann.lammpstrj的最后一步考出， 变成last.lammpstrj  
  '''

  # 定义源文件和目标文件的路径

  source_file_path = f_trj
  prefix = f_trj.split('ann.lammpstrj')[0]
  target_file_path = './operation_center.tmp_files/'+prefix+'last.lammpstrj'

  # 读取源文件内容
  with open(source_file_path, 'r') as file:
    lines = file.readlines()

  # 找到最后一个 "ITEM: TIMESTEP" 出现的位置
  index = len(lines) - 1  # 从文件末尾开始搜索
  for i in range(len(lines) - 1, -1, -1):  # 逆序遍历文件
    if "ITEM: TIMESTEP" in lines[i]:
        index = i  # 找到最后一个匹配的行
        break

  # 将从最后一个 "ITEM: TIMESTEP" 到文件末尾的内容写入目标文件
  with open(target_file_path, 'w') as target_file:
    for line in lines[index:]:
        target_file.write(line)

  return target_file_path 

def atoms_lpj_2_atoms_xyz(atoms_lpj, type_map):

  #1 get periodic table
  elements = mendeleev.get_all_elements()
  elementlist_ref=[]
  for element in elements:
    elementlist_ref.append(element.symbol)

  #2

  atoms_xyz = atoms_lpj.copy()
  eles_lpj=list(atoms_lpj.get_chemical_symbols())

  #eles_lpj -- > ele_xyz 
  eles_xyz = [type_map[elementlist_ref.index(x)] for x in eles_lpj]

  atoms_xyz.set_chemical_symbols(eles_xyz)

  return atoms_xyz


def get_type_map():

  if os.path.exists('uniq_element_list'):
    print("find uniq_element_list")
    type_map = np.loadtxt('uniq_element_list',ndmin=1,dtype=str)
    print("the elelist is ", type_map)

  else:

    from deepmd.infer import DeepPotential

    for file in os.listdir('.'):
      # 检查文件名是否以.pb结尾
      if file.endswith('.pb'):
        f_dp = os.path.join('.', file)
        break  # 找到第一个匹配的文件后退出循环    

    dp=DeepPotential(f_dp)
    type_map=dp.get_type_map()
    print(' '.join(type_map))
    with open ('uniq_element_list','w') as f:
      f.write(' '.join(type_map))  
   
    type_map = np.loadtxt('uniq_element_list',ndmin=1, dtype=str)

  return type_map


def gen_del_atm_files(del_atm_ids):

  '''
  生成删除原子序列文件
  '''
  n_del = 0
  for atm_id in del_atm_ids:

     n_del = n_del + 1
     atm_id_4_reader=atm_id+1

     np.savetxt('./'+dir_delinfo+'/index_to_del_'+str(n_del)+'.txt',[atm_id_4_reader],fmt='%d')

  np.savetxt('./'+dir_delinfo+'/natm_to_del.txt',[n_del],fmt='%d')


def gen_add_atm_files(add_eleidx_and_cors):

  '''
  生成增加原子序列及坐标文件
  '''

  n_add = 0
  for i in range(len(add_eleidx_and_cors)):

    n_add = n_add + 1

    np.savetxt('./'+dir_addinfo+'/eleidx_to_add_'+str(n_add)+'.txt',[add_eleidx_and_cors[i,0]],fmt='%d')
    np.savetxt('./'+dir_addinfo+'/x_to_add_'+str(n_add)+'.txt',[add_eleidx_and_cors[i,1]],fmt='%.3f')
    np.savetxt('./'+dir_addinfo+'/y_to_add_'+str(n_add)+'.txt',[add_eleidx_and_cors[i,2]],fmt='%.3f')
    np.savetxt('./'+dir_addinfo+'/z_to_add_'+str(n_add)+'.txt',[add_eleidx_and_cors[i,3]],fmt='%.3f')

  np.savetxt('./'+dir_addinfo+'/natm_to_add.txt',[n_add],fmt='%d')



if __name__ == "__main__":

  #*********************************#
  #input
  region = [131,181,145,190,80,140] #规则起效的区域 [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
  region = np.array(region).reshape(-1,2)
  natm_to_add = 6
  ele_idx = 3
  f_lpj = 'ann.lammpstrj'
  f_ae_lpj = 'ae.ann.lammpstrj' 
  #*********************************#
  '''
  data initialize
  '''
  mkdir('operation_center.tmp_files')
  mkdir('operation_center.add_info')
  mkdir('operation_center.del_info')

  f_lpj_last = make_f_lpjlast(f_lpj)  #将ann.lammpstrj的最后一步考出
  prefix = f_lpj_last.split('.lammpstrj')[0]

  atoms_lpj=read(f_lpj_last,format="lammps-dump-text")

  type_map = get_type_map()  #get type_map  
    
  atoms_xyz = atoms_lpj_2_atoms_xyz(atoms_lpj, type_map)
 

  '''
  增删原子规则：程序的主要部分，通过函数del_role 和 add_role 来产生指定增删原子规则，并由此产生需删除的原子序号和需增加的原子序号。可以是空，即没有增/删原子产生
  '''
  del_atm_ids = []
  add_eleidx_and_cors = [] 
  
  from del_rules import * 
  from add_rules import * 

  del_atm_ids = gen_del_info(atoms_xyz,region,ele='Na')     #del_infor 是一个需要删除原子的id_list
  add_eleidx_and_cors = gen_add_info(atoms_xyz, region, ele='Na', natm_to_add=7)
  print("add_eleidx_and_cors =", add_eleidx_and_cors)

  #特别的, 可以用和del_atm_ids=[]表示不增不减

  initialize()
  gen_del_atm_files(del_atm_ids)   #把要删除的原子信息写到del_tmp/中相关的文件
  gen_add_atm_files(add_eleidx_and_cors)    #把要增加的原子信息写到add_tmp/中相关的文件


