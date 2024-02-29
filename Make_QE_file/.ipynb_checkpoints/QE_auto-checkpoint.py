



def log_for_text(log_path,text):
    import datetime
    with open(log_path, 'a') as f:
        print(f"{datetime.datetime.now()} : {text}", file=f)


if __name__ == '__main__':
    import pandas as pd
    import os
    import subprocess
    import Make_input_file
    from contextlib import redirect_stdout
    excel_path="/home/hpc/QE_dir/PD_file/test_fluorides_struct-compositioin-feature_sisso.xlsx"
    test_fluorides=pd.read_excel(excel_path,index_col=0)
    list_Changed_Structure=test_fluorides[test_fluorides["M"]=="Fe"].index.values
    #output_file_path="/home/hpc/QE_dir/PD_file/inputfile"
    UPF_folder_path="/home/hpc/QE_dir/PD_file/UPF"
    cif_folder_path="/home/hpc/QE_dir/PD_file/cif MI_oxyarsenides_result"
    #list_Changed_Structure=[i for i in list_Changed_Structure if not "223223" in i]##325 skip
    log_path=f"/home/hpc/QE_dir/PD_file/logs/log_{datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')}.txt"
    if not os.path.exists("/home/hpc/QE_dir/PD_file/logs"):
        os.mkdir("/home/hpc/QE_dir/PD_file/logs")
    with redirect_stdout(open(log_path, 'a')):###
        for i in list_Changed_Structure:
            #log_for_text(log_path,f"{i}")
            print(f"{datetime.datetime.now()} : {i}")
            output_file_path=f"/home/hpc/QE_dir/{i}"
            if not os.path.exists(output_file_path):
                os.mkdir(output_file_path)
            os.chdir(output_file_path)
            print(f"{datetime.datetime.now()} : {i} : Parser_make_input_file(relax) : start")
            relax_input=Make_input_file.Parser_make_input_file(output_file_path,UPF_folder_path,cif_folder_path,i)
            relax_output=relax_input.split(".")
            relax_output[-1]="out"
            relax_output=".".join(relax_output)
            relax_done=False
            if os.path.exists(output_file_path+"/"+relax_output):
                with open(output_file_path+"/"+relax_output, 'r') as f:
                    data = f.readlines()
                for readline in data:
                    if "JOB DONE." in readline:
                        print(f"{datetime.datetime.now()} : {i} : relax : skip")
                        #log_for_text(log_path,f"{i} : relax : skip")
                        relax_done=True
            if not relax_done:
                print(f"{datetime.datetime.now()} : {i} : relax : start")
                #log_for_text(log_path,f"{i} : relax : start")
                command=f"mpirun -np 24 pw.x <{relax_input} >{relax_output}"
                ret = subprocess.run(command, shell=True, capture_output=True, text=True)
                if ret.returncode==0:
                    print(f"{datetime.datetime.now()} : {i} : relax : end")
                    #log_for_text(log_path,f"{i} : relax : end")
                with open(output_file_path+"/"+relax_output, 'r') as f:
                    data = f.readlines()
                next_ok=False
                for readline in data:
                    if "JOB DONE." in readline:
                        print(f"{datetime.datetime.now()} : {i} : relax : JOB DONE.")
                        #log_for_text(log_path,f"{i} : relax : JOB DONE.")
                        next_ok=True
                        continue
                if not next_ok:
                    print(f"{datetime.datetime.now()} : {i} : relax : Not JOB DONE.")
                    #log_for_text(log_path,f"{i} : relax : Not JOB DONE.")
                    continue
            path_in=output_file_path+"/"+relax_input
            path_out=output_file_path+"/"+relax_output
            next_calc="vc-relax"
            print(f"{datetime.datetime.now()} : {i} : Make_next_input(vc-relax) : start")
            vc_relax_input=Make_input_file.Make_next_input(path_in,path_out,next_calc)
            vc_relax_output=vc_relax_input.split(".")
            vc_relax_output[-1]="out"
            vc_relax_output=".".join(vc_relax_output)
            vc_relax_done=False
            if os.path.exists(output_file_path+"/"+vc_relax_output):
                with open(output_file_path+"/"+vc_relax_output, 'r') as f:
                    data = f.readlines()
                for readline in data:
                    if "JOB DONE." in readline:
                        print(f"{datetime.datetime.now()} : {i} : vc-relax : skip")
                        #log_for_text(log_path,f"{i} : vc-relax : skip")
                        vc_relax_done=True
            if not vc_relax_done:
                print(f"{datetime.datetime.now()} : {i} : vc-relax : start")
                #log_for_text(log_path,f"{i} : vc-relax : start")
                command=f"mpirun -np 24 pw.x <{vc_relax_input} >{vc_relax_output}"
                ret = subprocess.run(command, shell=True, capture_output=True, text=True)
                if ret.returncode==0:
                    print(f"{datetime.datetime.now()} : {i} : vc-relax : end")
                    #log_for_text(log_path,f"{i} : vc-relax : end")
                with open(output_file_path+"/"+vc_relax_output, 'r') as f:
                    data = f.readlines()
                next_ok=False
                for readline in data:
                    if "JOB DONE." in readline:
                        print(f"{datetime.datetime.now()} : {i} : vc-relax : JOB DONE.")
                        #log_for_text(log_path,f"{i} : vc-relax : JOB DONE.")
                        next_ok=True
                        continue
                if not next_ok:
                    print(f"{datetime.datetime.now()} : {i} : vc-relax : Not JOB DONE.")
                    #log_for_text(log_path,f"{i} : vc-relax : Not JOB DONE.")
            print(f"{datetime.datetime.now()} : ===================== Next ======================")
            #log_for_text(log_path,"===================== Next ======================")