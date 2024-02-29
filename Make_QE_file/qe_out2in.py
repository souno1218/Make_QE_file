
def qe_out2in(base_inputfile_path,base_outputfile_path,next_calc):
    #base_inputfile_path : path of inputfile of one previous calculation
    #base_outputfile_path : path of outputfile of one previous calculation
    #next_calc : Next calculation,'scf','nscf','bands','relax','md','vc-relax','vc-md'
    #base_inputfile_path : 一つ前の計算のinputfileのpath
    #base_outputfile_path : 一つ前の計算のoutputfileのpath
    #next_calc : 次の計算,'scf','nscf','bands','relax','md','vc-relax','vc-md'
    with open(base_outputfile_path, 'r') as f:
        data = f.readlines()
    JOB_DONE=False
    for readline in data:
        if "JOB DONE." in readline:
            JOB_DONE=True
    if not JOB_DONE:
        raise Exception("Input file cannot be created because the previous calculation has not been completed.")
    
    for i,readline in enumerate(data):
        if "ATOMIC_POSITIONS" in readline:
            final_ATOMIC_POSITIONS_index_start=i+1#最後のだけ保存
    
    if final_ATOMIC_POSITIONS_index_start==None:
        raise Exception("not found 'ATOMIC_POSITIONS' in output file")
        
    for i,readline in enumerate(data[final_ATOMIC_POSITIONS_index_start:]):#final_ATOMIC_POSITIONS_index_start以降
        if "End final coordinates" in readline:
            final_ATOMIC_POSITIONS_index_end=i+final_ATOMIC_POSITIONS_index_start
    if final_ATOMIC_POSITIONS_index_end==None:
        raise Exception("not found 'End final coordinates' in output file")
        
    ATOMIC_POSITIONS_list=[]
    for i in range(final_ATOMIC_POSITIONS_index_start,final_ATOMIC_POSITIONS_index_end):
        one_ATOMIC_POSITIONS=data[i].split()
        atom=one_ATOMIC_POSITIONS[0]
        abc=[str(round(float(k),5)) for k in one_ATOMIC_POSITIONS[1:4]]
        can_move_01=" ".join(one_ATOMIC_POSITIONS[4:])
        ATOMIC_POSITIONS_list.append([atom,*abc,can_move_01])
    max_len_atom=max([len(i[0]) for i in ATOMIC_POSITIONS_list])
    max_len_a=max([len(i[1]) for i in ATOMIC_POSITIONS_list])
    max_len_b=max([len(i[2]) for i in ATOMIC_POSITIONS_list])
    max_len_c=max([len(i[3]) for i in ATOMIC_POSITIONS_list])
    final_ATOMIC_POSITIONS=[]
    for i in range(len(ATOMIC_POSITIONS_list)):
        atom=ATOMIC_POSITIONS_list[i][0].ljust(max_len_atom)
        a=ATOMIC_POSITIONS_list[i][1].rjust(max_len_a)
        b=ATOMIC_POSITIONS_list[i][2].rjust(max_len_b)
        c=ATOMIC_POSITIONS_list[i][3].rjust(max_len_c)
        can_move_01=ATOMIC_POSITIONS_list[i][4]
        final_ATOMIC_POSITIONS.append(f" {atom}  {a}  {b}  {c}  {can_move_01}\n")
    
    with open(base_inputfile_path, 'r') as f:
        data = f.readlines()
    found_calc,found_ATOMIC_POSITIONS,found_K_POINTS=False,False,False
    for i,readline in enumerate(data):
        if "calculation" in readline:
            data[i]=f" calculation = '{next_calc}'\n"
        if "ATOMIC_POSITIONS" in readline:
            base_ATOMIC_POSITIONS_index_start=i+1
    if not found_calc:
        raise Exception("not found 'calculation' in input file")
    if not found_ATOMIC_POSITIONS:
        raise Exception("not found 'ATOMIC_POSITIONS' in input file")
    for i,readline in enumerate(data[base_ATOMIC_POSITIONS_index_start:]):
        if "K_POINTS" in readline:
            base_ATOMIC_POSITIONS_index_end=i+base_ATOMIC_POSITIONS_index_start
    if not found_K_POINTS:
        raise Exception("not found 'K_POINTS' in input file")
    
    data[base_ATOMIC_POSITIONS_index_start:base_ATOMIC_POSITIONS_index_end]=final_ATOMIC_POSITIONS

    next_input_path=base_inputfile_path.split(".")
    next_input_path[-2]=next_calc
    next_input_path=".".join(next_input_path)
    
    with open(next_input_path, 'w+') as f:
        f.writelines(data)
    print("qe_out2in : Done")