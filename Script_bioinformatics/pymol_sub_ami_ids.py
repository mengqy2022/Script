
from pymol import *
#from pymol.cgo import *

residue_id = []
print("####################################################################################################################")
print("# Usage: sele ligand ---> in 'A(action)' ---> modify ---> around ---> select residues with in 4 (or any numbers) A #")
print("# or select 4A_residues, byres sele around 4, and Protein-Protein interaction residues repeatet twice command      #")
print("# Note: the script 'list residue_id' was clear for many select command, so ...                                     #")
print("#####################################################################################################################")

def get_residue_id(selection = "(sele)"):
    
    cmd.hide("spheres")
    
    cmd.show("spheres", selection)
    
    # cmd.iterate obtain all residues id, but the id is repeated many times
    cmd.iterate(selection, "residue_id.append(int(resi))")
    # cmd.iterate(selection, "residue_id.append(resn+resi)")
    
    # new_ids and sorted_ids were obtain unique residue id
    new_ids = []
    new_ids = list(set(residue_id))
    print(sorted(new_ids))
    sorted_ids = sorted(new_ids)\
    
    print("all residues: ", len(sorted_ids))
    
    # the new_ids was clear and list covert to str
    new_ids = []
    for i in sorted_ids:
        
        new_ids.append(str(i))
        
    #print(new_ids)
    ids = "; ".join(new_ids)
    print("Please check your operation was correct!!!")
    print("copy the following line into the mmpbsa.in")
    print("; ".join(new_ids))
    
    
#    residues_ids = list(set(residues_ids))
#    for i in residues_ids[:-1]:
#        print(i, end = "    ")
    
    return

cmd.extend("get_residue_id", get_residue_id)

# 原因: 在使用gmx_MMPBSA.py这个程序的时候，偶尔会出现一些问题。特别是使用charmm力场的时候。
# 由此就老实的使用ambertools的MMPBSA.py程序，但该程序计算熵也没有gmx_MMPBSA程序比较完善，同时不能指定配体氨基酸附近的氨基酸，特别是如果是蛋白-蛋白体系，所以借助pymol获取蛋白-蛋白体系界面的氨基酸编号，并进行MMPBSA计算。

# 备注：
# 1. 由于residue_id列表写在函数里面一直报错，所以写在外面，这样容易引起该列表没有清空。所以，选错指定范围的氨基酸残基，可以不重新运行。但如果选错配体，请重新运行脚本。如果是蛋白-蛋白体系，选择配体，再选择范围，再运行该程序，就获取受体的氨基酸编号。同理，选中受体，就是获得配体的氨基酸编号。
# 2. 切记ambertools程序本身就习惯将体系编号从1开始，所以选中之前可以使用pdb4amber重新编号。
# 3. 代码19行，本意是获取氨基酸+编号的，但是没有试运行，有兴趣再修改。
# 4. 使用方法可以参考图片。