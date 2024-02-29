#unsupported : ibrav=5,-5,9~
#pymatgen,numpy,pandas are required

#memo
## https://qiita.com/ojiya/items/45b142d34e6fccf2ba1c
## https://pymatgen.org/pymatgen.io.html#module-pymatgen.io.pwscf
#https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm224
#https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-wp-list?gnum=74
#https://en.wikipedia.org/wiki/List_of_space_groups
def cif2qe_in(import_cif_path,calc,input_file_name,output_folder_path,UPF_folder_path,change_site=[],site2dummy_label=[]):
    #import_cif_path : Absolute path of the original pre-converted cif file
    #calc : 'scf','nscf','bands','relax','md','vc-relax','vc-md'
    #input_file_name : 368_FeAs_NaMnF_22212 like, without .relax.in
    #output_folder_path : Absolute path input_file_name is not needed
    #UPF_folder_path : Absolute path of the UPF folder
    #change_site : Description as ["Na2->Mg2", "F1Z->Xx1"], virtual atoms are written as Xx1
    #site2dummy_label : Description as ["Xx1->O2F3"], Used for UPF names

    #import_cif_path : 元となる変換まえのcifファイルの絶対path
    #calc : 'scf','nscf','bands','relax','md','vc-relax','vc-md'
    #input_file_name : 368_FeAs_NaMnF_22212 みたいな、.relax.inなし
    #output_folder_path : 絶対path,input_file_nameはいらない
    #UPF_folder_path : UPFフォルダーの絶対path
    #change_site : ["Na2->Mg2","F1Z->Xx1"]のように記述、仮想原子はXx1のように記述
    #site2dummy_label : ["Xx1->O2F3"]のように記述、仮想原子のラベル,UPFの名前に使用
    
    from pymatgen.io.cif import CifParser
    from pymatgen.core.periodic_table import Element
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pathlib import Path
    import numpy as np
    import warnings
    import os 
    import re
    import pandas as pd
    def drop_num(string):
        return "".join([i for i in string if not i.isdigit()])
        
    warnings.resetwarnings()
    
    UPF_FOLDER = Path(UPF_folder_path)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        mat = CifParser(import_cif_path).parse_structures(symmetrized=False)[0]
    mat_sym = SpacegroupAnalyzer(mat).get_symmetrized_structure()
    n_sites=len(mat.as_dict()["sites"])
    
    base_element_labels_list=list(dict.fromkeys([mat.as_dict()["sites"][i]['label'] for i in range(n_sites)]))
    base_element_labels_dict={i:j for i,j in enumerate(base_element_labels_list)}
    
    site_list=[mat.as_dict()["sites"][i]['label'] for i in range(n_sites)]
    labels_index_dict={}
    for i in base_element_labels_dict.values():
        labels_index_dict[i]=[k for k,l in enumerate(site_list) if l==i]
        
    if not np.all([(j in [k.split("->")[1] for k in change_site]) for j in [i.split("->")[0] for i in site2dummy_label]]):
        raise Exception("Atomic replacement, not consistent with label")
    if not np.all([i in base_element_labels_dict.values() for i in [i.split("->")[0] for i in change_site]]):
        if np.all([i in base_element_labels_dict.values() for i in [drop_num(i.split("->")[0]) for i in change_site]]):
            change_site_dict={drop_num(i.split("->")[0]):i.split("->")[1] for i in change_site}
            site2dummy_label_dict={i.split("->")[0]:i.split("->")[1] for i in site2dummy_label}
        else:
            raise Exception("Atomic replacement, You are trying to replace what is not there")
    else:
        change_site_keys=[i.split("->")[0] for i in change_site]
        change_site_values=[i.split("->")[1] for i in change_site]
        site2dummy_label_keys=[i.split("->")[0] for i in site2dummy_label]
        site2dummy_label_values=[i.split("->")[1] for i in site2dummy_label]
        for i,j in enumerate(set(site2dummy_label_values)):
            if len(set(site2dummy_label_values))==1:
                i=""
            same_index=np.where(j==np.array(site2dummy_label_values))[0]
            before_change=[site2dummy_label_keys[same_index[k]]for k in same_index]
            txt=re.sub(r"[^a-zA-Z]","",site2dummy_label_keys[same_index[0]])
            for k in same_index:
                site2dummy_label_keys[k]=f"{txt}{i}"
            same_index=np.array([np.where(k==np.array(change_site_values))[0] for k in before_change]).reshape(-1)
            for k in same_index:
                change_site_values[k]=f"{txt}{i}"
        change_site_dict={change_site_keys[i]:change_site_values[i] for i in range(len(change_site_values))}
        site2dummy_label_dict={site2dummy_label_keys[i]:site2dummy_label_values[i] for i in range(len(site2dummy_label_values))}
        
    changed_elements_dict={}
    for i,j in enumerate(base_element_labels_dict.values()):
        if not j in change_site_dict.keys():
            changed_elements_dict[i]=mat.as_dict()["sites"][labels_index_dict[j][0]]['species'][0]['element']
        elif change_site_dict[j] in site2dummy_label_dict.keys():
            changed_elements_dict[i]=change_site_dict[j]
        else:
            try:
                Element(drop_num(change_site_dict[j]))
                changed_elements_dict[i]=drop_num(change_site_dict[j])
            except:
                raise Exception(f"\{change_site_dict[j]} cannot get mass in pymatgen and is not set as a virtual atom")
    changed_elements=[]
    for i in list(dict.fromkeys(changed_elements_dict.values())):
        if i in site2dummy_label_dict.keys():
            changed_elements.append(site2dummy_label_dict[i])
        else:
            changed_elements.append(i)
            
    #control
    control_txt = "!---- starting input file ----\n"\
    "&control\n"\
    " calculation = '"+str(calc)+"'\n"\
    " prefix='"+str(input_file_name)+"',\n"\
    " tstress = .true.\n"\
    " tprnfor = .true.\n"\
    " pseudo_dir = '"+str(UPF_folder_path)+"',\n"\
    " outdir='./work/'\n"\
    " etot_conv_thr = 1.d-5\n"\
    " forc_conv_thr = 1.d-4\n"\
    " disk_io='low'\n"\
    " wf_collect=.true.\n"\
    "/\n"
    
    # system
    ibrav_symbol,ibrav=find_ibrav(mat)
    if None==ibrav_symbol:
        return None
    celldm_str=f" celldm(1) = {round(1.88972616*mat.lattice.a, 5)}\n"
    if ibrav in [8,9,91,10,11,12,-12,13,-13,14]:#celldm_2
        celldm_str+=f" celldm(2) = {round(mat.lattice.b/mat.lattice.a, 5)}\n"
    if ibrav in [4,6,7,8,9,91,10,11,12,-12,13,-13,14]:#celldm_3
        celldm_str+=f" celldm(3) = {round(mat.lattice.c/mat.lattice.a, 5)}\n"
    #if ibrav in []:#celldm_4~~
    
    nat=int(len(mat.sites)/2) if ibrav_symbol in ["C","I"] else len(mat.sites)  #Assumed to be the same as the number of columns in ATOMIC_POSITIONS,
                                                                                #ibrav_symbol=F,R not implemented
    
    ntyp= len(changed_elements)#Assumed to be the same as the number of columns in ATOMIC_SPECIES
    UPF_list=[]
    for i in changed_elements:
        UPF=list(UPF_FOLDER.glob(f"{i}.*.UPF"))
        if len(UPF)==0:
            raise Exception(f"UPF not found \n{i} has no UPF")
        UPF_list.append(UPF_folder_path+"/"+str(UPF[0].name))
    cutoff_wavefunctions=[]
    cutoff_charge_density=[]
    for i,UPF_path in enumerate(UPF_list):
        with open(UPF_path, 'r') as f:
            data = f.readlines()
            for j in data:
                if "cutoff for wavefunctions:" in j:
                    for k in j.split():
                        try:
                          cutoff_wavefunctions.append(float(k))
                        except:
                            None
                elif "cutoff for charge density:" in j:
                    for k in j.split():
                        try:
                          cutoff_charge_density.append(float(k))
                        except:
                            None
        if len(cutoff_wavefunctions)!=i+1:
            raise Exception(f"cutoff_wavefunctions not found\nElement:{changed_elements[i]} , UPF_path:{UPF_path}")
        elif len(cutoff_charge_density)!=i+1:
            raise Exception(f"cutoff_charge_density not found\nElement:{changed_elements[i]} , UPF_path:{UPF_path}")
    ecutwfc=max(cutoff_wavefunctions)
    ecutrho=max(cutoff_charge_density)
    
    system_txt = "&system\n"\
    " ibrav = "+str(ibrav)+",\n"\
    ""+celldm_str+""\
    " nat = "+str(nat)+",\n"\
    " ntyp = "+str(ntyp)+",\n"\
    " ecutwfc = "+str(ecutwfc)+",\n"\
    " ecutrho = "+str(ecutrho)+",\n"\
    " occupations = 'smearing'\n"\
    " smearing = 'm-p'\n"\
    " degauss = 0.01\n"\
    "/\n"
    
    #electrons
    electrons_txt = "&electrons\n"\
    " electron_maxstep = 100\n"\
    " mixing_beta = 0.3\n"\
    " conv_thr = 1.0d-8\n"\
    "/\n"
    
    #ions
    #ATOMIC_SPECIES
    df_atomic_mass=pd.read_csv(os.path.dirname(os.path.abspath("__file__"))+"/dummy_species.csv", index_col=0)
    ATOMIC_SPECIES_list=[]
    for i in range(len(changed_elements)):
        atom=str(list(dict.fromkeys(changed_elements_dict.values()))[i])
        try:#Setting mass to -1 automatically sets the standard atomic weight
            mass_num=Element(changed_elements[i]).data['Atomic mass']
        except ValueError as Ve:
            if changed_elements[i] in df_atomic_mass.index.values:
                mass_num=df_atomic_mass.loc[changed_elements[i],"Atomic_mass"]
            else:
                raise Exception(f"Atomic volume cannot be determined\n{i} cannot get mass in pymatgen and is not set as a virtual atom")
        mass=str(round(mass_num,3))#.rjust(7)
        UPF=str(UPF_list[i]).split("/")[-1]
        ATOMIC_SPECIES_list.append([atom,mass,UPF])
    max_len_atom=max([len(i[0]) for i in ATOMIC_SPECIES_list])
    max_len_mass=max([len(str(i[1])) for i in ATOMIC_SPECIES_list])
    ATOMIC_SPECIES=""
    for i in range(len(changed_elements)):
        atom=ATOMIC_SPECIES_list[i][0].ljust(max_len_atom)
        mass=ATOMIC_SPECIES_list[i][1].rjust(max_len_mass)
        UPF=ATOMIC_SPECIES_list[i][2]
        ATOMIC_SPECIES+=f" {atom} {mass}    {UPF}\n"
        
    #ATOMIC_POSITIONS
    ATOMIC_POSITIONS_list=[]
    lattice_abc=np.array([mat.lattice.a,mat.lattice.b,mat.lattice.c])/mat.lattice.a
    def judge_isclose(x):
        return ~np.any(np.isclose(x,[[2.**i,1-2.**i,0.5-2.**i,0.5+2.**i,0] for i in np.arange(-5,0)]))
    judge_isclose=np.frompyfunc(judge_isclose, 1, 1)

    for index in range(len(base_element_labels_dict)):
        label=base_element_labels_dict[index]
        atom=changed_elements_dict[index]
        atom_index=labels_index_dict[label]
        abc_list=np.array([[mat_sym.sites[j].a,mat_sym.sites[j].b,mat_sym.sites[j].c] for j in atom_index])
        if ibrav_symbol in ["C","I"]:
            abc_list=abc_drop(abc_list)
        for j in range(abc_list.shape[0]):
            can_move=judge_isclose(abc_list[j]).astype("bool")
            abc_list[j][can_move]=abc_list[j][can_move]*lattice_abc[can_move]
            abc=[str(round(k,5)) for k in abc_list[j]]
            can_move_01="".join(["1 " if k else "0 " for k in can_move])
            ATOMIC_POSITIONS_list.append([atom,*abc,can_move_01])
    max_len_atom=max([len(i[0]) for i in ATOMIC_POSITIONS_list])
    max_len_a=max([len(i[1]) for i in ATOMIC_POSITIONS_list])
    max_len_b=max([len(i[2]) for i in ATOMIC_POSITIONS_list])
    max_len_c=max([len(i[3]) for i in ATOMIC_POSITIONS_list])
    ATOMIC_POSITIONS=""
    for i in range(len(ATOMIC_POSITIONS_list)):
        atom=ATOMIC_POSITIONS_list[i][0].ljust(max_len_atom)
        a=ATOMIC_POSITIONS_list[i][1].rjust(max_len_a)
        b=ATOMIC_POSITIONS_list[i][2].rjust(max_len_b)
        c=ATOMIC_POSITIONS_list[i][3].rjust(max_len_c)
        can_move_01=ATOMIC_POSITIONS_list[i][4]
        ATOMIC_POSITIONS+=f" {atom}  {a}  {b}  {c}  {can_move_01}\n"

    ions_txt = "&ions\n"\
    "/\n"\
    "ATOMIC_SPECIES\n"+ATOMIC_SPECIES+""\
    "ATOMIC_POSITIONS {alat}\n"+ATOMIC_POSITIONS+""\
    "K_POINTS {automatic}\n"\
    " 12 12 4 0 0 0\n"\
    
    txt = control_txt+system_txt+electrons_txt+ions_txt+"\n"
    with open(output_folder_path+"/"+input_file_name+"."+calc+".in",'wb') as output_file:
        output_file.write(txt.encode('utf-8'))
    print("cif2qe_in : Done")

def abc_drop(abc_list):
    import numpy as np
    index_sort=np.argsort(np.linalg.norm(abc_list,axis=1,ord=2))[::-1]
    abc_list=abc_list[index_sort]
    abc_list_plushalf=(abc_list+0.5)%1
    drop_list=[]
    for i in range(abc_list.shape[0]):
        if i in drop_list:
            continue
        same_index=np.where(np.all(np.isclose(abc_list,abc_list_plushalf[i]),axis=1))[0]
        if len(same_index)!=0:
            drop_list.append(same_index[0])
    return np.delete(abc_list, drop_list, 0)

def find_ibrav(mat):
    short_name=mat.get_space_group_info()[0]#P4/nmm
    number=mat.get_space_group_info()[1]#129
    
    if number<=2:
        crystal_system="Triclinic"
    elif number<=15:
        crystal_system="Monoclinic"
    elif number<=74:
        crystal_system="Orthorhombic"
    elif number<=142:
        crystal_system="Tetragonal"
    elif number<=167:
        crystal_system="Trigonal"
    elif number<=194:
        crystal_system="Hexagonal"
    elif number<=230:
        crystal_system="Cubic"
    else:
        raise Exception(f"crystal structure {mat.get_space_group_info()} is not yet implemented")
        return None,None
        
    if "P" in short_name:
        ibrav_symbol="P"   #primitiv
    elif "I" in short_name:
        ibrav_symbol="I"   #body centered
    elif "F" in short_name:
        ibrav_symbol="F"   #face centered
    elif "A" in short_name:
        ibrav_symbol="A"   #centered on A faces only
    elif "B" in short_name:
        ibrav_symbol="B"   #centered on B faces only
    elif "C" in short_name:
        ibrav_symbol="C"   #centered on C faces only
    elif "R" in short_name:
        ibrav_symbol="R"   #rhombohedral
    else:
        raise Exception(f"crystal structure {mat.get_space_group_info()} is not yet implemented")
        return None,None
        
    ibrav_dict={"Cubic":{"P":1,"F":2,"I":3},
                "Hexagonal":{"P":4},
                #"Trigonal":{"R":5},#unsupported
                "Tetragonal":{"P":6,"I":7,},
                "Orthorhombic":{"P":8},}
                #"Monoclinic":{"P":12},#unsupported
                #"Triclinic":{"P":14}}#unsupported
    #unsupported : 5,-5,9~
    
    try:
        ibrav=ibrav_dict[crystal_system][ibrav_symbol]
    except Exception as e:
        raise Exception(f"crystal structure {mat.get_space_group_info()} is not yet implemented")
        return None,None
    return ibrav_symbol,ibrav