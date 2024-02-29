



def qe_out2cif(qe_input_path,qe_output_path,base_cif_path,output_cif_path):
    #qe_input_path : First inputfile, e.g. relax->vc-relax, if calculated in that order, use relax.in
    #qe_output_path : Final outputfile
    #base_cif_path : Absolute path of the original pre-converted cif file
    #output_cif_path : path and name of the cif to be made from the final outputfile
    #qe_input_path : 一番最初のinputfile, 例えばrelax->vc-relaxの順に計算したなら、relax.inを使用
    #qe_output_path : 最終的なoutputfile
    #base_cif_path : 元となる変換まえのcifファイルの絶対path
    #output_cif_path : 最終的なoutputfileから作るcifのpathとname
    import numpy as np
    import datetime
    import os
    base_positions,length_abc=base_cif2base_positions(base_cif_path)
    qe_input_num,ibrav=qe_input2num(qe_input_path,length_abc)
    cif2qe_index,plus_05_list=make_cif2qe_index(base_positions,qe_input_num)
    index2element,qe_out_num,end_length_abc=qe_output2num(qe_output_path,length_abc,ibrav)
    cif_name=os.path.basename(os.path.splitext(output_cif_path)[0])
    
    with open(base_cif_path,'r') as f:
        data = f.readlines()
        data[0]="data_"+cif_name+"\n"
        for i,readline in enumerate(data):
            if "_audit_creation_date" in readline:
                data[i]=f"_audit_creation_date              {datetime.datetime.now().strftime('%Y-%m-%d')}\n"
            if "_cell_length_a" in readline:
                data[i]=f"_cell_length_a                    {round(end_length_abc[0],5)}\n"
            if "_cell_length_b" in readline:
                data[i]=f"_cell_length_b                    {round(end_length_abc[1],5)}\n"
            if "_cell_length_c" in readline:
                data[i]=f"_cell_length_c                    {round(end_length_abc[2],5)}\n"
            if len(readline.split())==8:
               del data[i:]
        element_list=[]
        num_text_list=[]
        for tf_index,i in enumerate(cif2qe_index):
            element_list.append(index2element[i])
            if plus_05_list[tf_index]:
                qe_out_num[i]=(qe_out_num[i]+0.5)%1
                qe_out_num[i][np.isclose(qe_out_num[i],1)]=0
            num_text_list.append([f"{round(j,5):.5f}" for j in qe_out_num[i]])
        a_max_len=max([len(i[0])for i in num_text_list])
        b_max_len=max([len(i[1])for i in num_text_list])
        c_max_len=max([len(i[2])for i in num_text_list])
        element_name_list=element_list.copy()
        for i in list(set(element_list)):
            index_list=[j for j, x in enumerate(element_list) if x == i]
            for j,index in enumerate(index_list):
                element_name_list[index]=f"{element_list[index]}{j+1}"
        for i in range(len(element_list)):
            text=f"{element_name_list[i].ljust(3)}    {element_list[i].ljust(2)}    "
            text+=f"{num_text_list[i][0].ljust(a_max_len)}   {num_text_list[i][1].ljust(b_max_len)}   {num_text_list[i][2].ljust(c_max_len)}   "
            text+="0.00000  Uiso   1.00\n"
            data.append(text)
    with open(output_cif_path, 'w+') as f:
        f.writelines(data)
    print("qe_out2cif : Done")




def base_cif2base_positions(base_cif_path):
    import numpy as np
    with open(base_cif_path,'r') as f:
        base_cif_data = f.readlines()
        base_positions=[]
        for readline in base_cif_data:
            if "cell_length_a" in readline:
                length_a=float(readline.split()[1])
            if "cell_length_b" in readline:
                length_b=float(readline.split()[1])
            if "cell_length_c" in readline:
                length_c=float(readline.split()[1])
            if len(readline.split())==8:
                base_positions.append(readline.split())
    base_positions_num=[]
    for i in range(len(base_positions)):
        num_data=np.array([float(base_positions[i][j]) for j in range(2,5)])#index=2
        base_positions_num.append(num_data)
    return np.array(base_positions_num),np.array([length_a,length_b,length_c])


def qe_input2num(qe_input_path,length_abc):
    import numpy as np
    with open(qe_input_path,'r') as f:
        qe_input_data = f.readlines()
        qe_input_num=[]
        for readline in qe_input_data:
            if "ibrav" in readline:
                ibrav=int(readline.split()[2].replace(",",""))
            if len(readline.split())==7:
                qe_input_num.append(np.array([float(j) for j in readline.split()[1:4]]))
    return np.array(qe_input_num)/(length_abc/length_abc[0]),ibrav

def make_cif2qe_index(base_positions,qe_input_num):
    import numpy as np
    output_len=len(qe_input_num)
    cif2qe_index=[]
    plus_05_list=[]
    for i in base_positions:
        plus_05=False
        serch_index=np.all(np.isclose(qe_input_num-i,0,atol=0.0005),axis=1)
        if np.sum(serch_index)>1:
            raise Exception("make_cif2qe_index\nserch_index too many")
        if np.sum(serch_index)==0:
            qe_input_num_plushalf=(qe_input_num+0.5)%1
            qe_input_num_plushalf[np.isclose(qe_input_num_plushalf,1,atol=0.0005)]=0
            serch_index=np.all(np.isclose(qe_input_num_plushalf-i,0,atol=0.0005),axis=1)
            plus_05=True
            if np.sum(serch_index)==0:
                raise Exception("make_cif2qe_index\nserch_index not found")
        cif2qe_index.append(np.where(serch_index)[0][0])
        plus_05_list.append(plus_05)
    return cif2qe_index,plus_05_list


def qe_output2num(qe_output_path,length_abc,ibrav):
    import numpy as np
    with open(qe_output_path,'r') as f:
        qe_output_data = f.readlines()
        a_posi_index=[]
        c_axes_index=[]
        for i,readline in enumerate(qe_output_data):
            if "ATOMIC_POSITIONS" in readline:
                a_posi_index.append(i)
            if "End final coordinates" in readline:
                End_final_coordinates_index=i
            if "crystal axes" in readline:
                c_axes_index.append(i)
        end_data=qe_output_data[max(a_posi_index)+1:End_final_coordinates_index]
        if ibrav==7:
            crystal_axes_data=qe_output_data[max(c_axes_index)+1:max(c_axes_index)+4]
            crystal_axes_data=np.array([float(crystal_axes_data[i].split()[3+i]) for i in range(3)])
            crystal_axes_data*=2
        elif ibrav==6:
            crystal_axes_data=qe_output_data[max(c_axes_index)+1:max(c_axes_index)+4]
            crystal_axes_data=np.array([float(crystal_axes_data[i].split()[3+i]) for i in range(3)])
        else:
            raise Exception(f"unknown ibrav\nibrav : {ibrav}")
        index2element=[]
        qe_out_num=[]
        for i,data in enumerate(end_data):
            one_base_data=data.split()
            index2element.append(one_base_data[0])
            qe_out_num.append(np.array([float(one_base_data[j]) for j in range(1,4)]))
        qe_out_num=np.array(qe_out_num)/crystal_axes_data
        end_length_abc=length_abc[0]*crystal_axes_data
    return index2element,qe_out_num,end_length_abc
