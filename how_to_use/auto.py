from Make_QE_file import make_input, cif_to_params, check_JOB_DONE
from Make_QE_file import output_to_params, write_new_cif, output_to_nbnd
from Make_QE_file import get_EFermi, get_highest_occupied, plot_band, plot_pdos
import subprocess, os, datetime
import pandas as pd

# pip uninstall Make_0E_file -y
# pip install git+https://github.com/souno1218/Make_0E_file.git

os.chdir("/home/hpc/QE_dir")

df = pd. read_excel("base_data.xlsx")

list_structure_name = [i[:-2] if (i.split("_")[-1] in ["y", "n"]) else i for i in list(df ["name"].values)]
dict_phase = {list_structure_name[i]:df["phase"].values[i] for i in range(len(list_structure_name))}

ls_file_name = os.listdir("old_cif_dir")

dict_brilloin_zone_path = {212 : ["gG"],
                           325 : ["gG"],
                           426 : ["gG"]}

dict_color_dict = {}

for i in range(len(list_structure_name)):
    dict_color_dict[list_structure_name[i]] = {df.loc[i, "O"] : "c"}

dict_ibrav = {22212 : 7, 22325 : 7, 22426 : 6}

k_fineness_Magnification = 20

pseudo_dir = "/home/hpc/QE_dir/UPF"
temp_scf = "/home/hpc/QE_dir/templates/template_smearing.in"
temp_projwfc = "/home/hpc/QE_dir/templates/template.projwfc.in"
temp_dos = "/home/hpc/QE_dir/templates/template.dos.in"
temp_band_x = "/home/hpc/QE_dir/templates/template.band_x.in"

date = str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")).center(
    max([len(i) for i in list_structure_name]), "="
)
print(f"=================={date}==================", flush=True)
for structure_name in list_structure_name:
    str_structure_name = structure_name.center(max([len(i) for i in list_structure_name]), "=")
    print(f"=================={str_structure_name}==================", flush=True)

    structure_dir = f"/home/hpc/QE_dir/{structure_name}"
    if not os.path.exists(structure_dir):
        os.makedirs(structure_dir)

    import_cif_path = f"{structure_dir}/{structure_name}.cif"
    params_structure = cif_to_params(import_cif_path)
    prefix = structure_name
    ################################################################################################################
    if not os.path.exists(f"{structure_dir}/relax"):
        os.makedirs(f"{structure_dir}/relax")
    os.chdir(f"{structure_dir}/relax")
    ################################################################################################################
    calc = "relax"
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)
    input_path = f"{structure_dir}/relax/{structure_name}.relax.in"
    output_path = f"{structure_dir}/relax/{structure_name}.relax.out"
    if not check_JOB_DONE(calc, output_path):
        make_input(
            calc,
            input_path,
            prefix,
            pseudo_dir=pseudo_dir,
            params_structure=params_structure,
            k_fineness_Magnification=k_fineness_Magnification,
            template_path=temp_scf,
        )
        command = f"mpirun -np 32 pw.x <{input_path} >{output_path}"
        ret = subprocess.run(command, shell=True, capture_output=True, text=True)
    if check_JOB_DONE(calc, output_path):
        params_structure = output_to_params(calc, output_path, params_structure)
    else:
        print("停止", flush=True)
        continue
    ################################################################################################################
    calc = "vc-relax"
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)
    input_path = f"{structure_dir}/relax/{structure_name}.vc-relax.in"
    output_path = f"{structure_dir}/relax/{structure_name}.vc-relax.out"
    if not check_JOB_DONE(calc, output_path):
        make_input(
            calc,
            input_path,
            prefix,
            pseudo_dir=pseudo_dir,
            params_structure=params_structure,
            k_fineness_Magnification=k_fineness_Magnification,
            template_path=temp_scf,
        )
        command = f"mpirun -np 32 pw.x <{input_path} >{output_path}"
        ret = subprocess.run(command, shell=True, capture_output=True, text=True)
    if check_JOB_DONE(calc, output_path):
        params_structure = output_to_params(calc, output_path, params_structure)
        cif_out_path = f"{structure_dir}/relax/{structure_name}.cif"
        write_new_cif(cif_out_path, structure_name, params_structure)
    else:
        print("停止", flush=True)
        continue
    ################################################################################################################
    if not os.path.exists(f"{structure_dir}/nscf"):
        os.makedirs(f"{structure_dir}/nscf")
    os.chdir(f"{structure_dir}/nscf")
    ################################################################################################################
    calc = "scf"
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)
    input_path = f"{structure_dir}/nscf/{structure_name}.scf.in"
    output_path = f"{structure_dir}/nscf/{structure_name}.scf.out"
    if not check_JOB_DONE(calc, output_path):
        make_input(
            calc,
            input_path,
            prefix,
            pseudo_dir=pseudo_dir,
            params_structure=params_structure,
            k_fineness_Magnification=k_fineness_Magnification,
            template_path=temp_scf,
        )
        command = f"mpirun -np 32 pw.x <{input_path} >{output_path}"
        ret = subprocess.run(command, shell=True, capture_output=True, text=True)
    if check_JOB_DONE(calc, output_path):
        nbnd = 2 * output_to_nbnd(output_path)
    else:
        print("停止", flush=True)
        continue
    ################################################################################################################
    calc = "nscf"
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)
    input_path = f"{structure_dir}/nscf/{structure_name}.nscf.in"
    output_path = f"{structure_dir}/nscf/{structure_name}.nscf.out"
    if not check_JOB_DONE(calc, output_path):
        make_input(
            calc,
            input_path,
            prefix,
            pseudo_dir=pseudo_dir,
            params_structure=params_structure,
            k_fineness_Magnification=k_fineness_Magnification,
            nbnd=nbnd,
            template_path=temp_scf,
        )
        command = f"mpirun -np 32 pw.x <{input_path} >{output_path}"
        ret = subprocess.run(command, shell=True, capture_output=True, text=True)
    if not check_JOB_DONE(calc, output_path):
        print("停止", flush=True)
        continue
    ################################################################################################################
    calc = "projwfc"
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)
    input_path = f"{structure_dir}/nscf/{structure_name}.projwfc.in"
    output_path = f"{structure_dir}/nscf/{structure_name}.projwfc.out"
    if not check_JOB_DONE(calc, output_path):
        make_input(calc, input_path, prefix, template_path=temp_projwfc)
        command = f"projwfc.x <{input_path} >{output_path}"
        ret = subprocess.run(command, shell=True, capture_output=True, text=True)
    if not check_JOB_DONE(calc, output_path):
        print("停止", flush=True)
        continue
    ################################################################################################################
    calc = "dos"
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)
    input_path = f"{structure_dir}/nscf/{structure_name}.dos.in"
    output_path = f"{structure_dir}/nscf/{structure_name}.dos.out"
    if not check_JOB_DONE(calc, output_path):
        make_input(calc, input_path, prefix, template_path=temp_dos)
        command = f"dos.x <{input_path} >{output_path}"
        ret = subprocess.run(command, shell=True, capture_output=True, text=True)
    if not check_JOB_DONE(calc, output_path):
        print("停止", flush=True)
        continue
    ################################################################################################################
    if not os.path.exists(f"{structure_dir}/bands"):
        os.makedirs(f"{structure_dir}/bands")
    os.chdir(f"{structure_dir}/bands")
    ################################################################################################################
    calc = "scf"
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)
    input_path = f"{structure_dir}/bands/{structure_name}.scf.in"
    output_path = f"{structure_dir}/bands/{structure_name}.scf.out"
    if not check_JOB_DONE(calc, output_path):
        make_input(
            calc,
            input_path,
            prefix,
            pseudo_dir=pseudo_dir,
            params_structure=params_structure,
            k_fineness_Magnification=k_fineness_Magnification,
            template_path=temp_scf,
        )
        command = f"mpirun -np 32 pw.x <{input_path} >{output_path}"
        ret = subprocess.run(command, shell=True, capture_output=True, text=True)
    if check_JOB_DONE(calc, output_path):
        nbnd = 2 * output_to_nbnd(output_path)
    else:
        print("停止", flush=True)
        continue
    ################################################################################################################
    calc = "bands"
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)
    input_path = f"{structure_dir}/bands/{structure_name}.bands.in"
    output_path = f"{structure_dir}/bands/{structure_name}.bands.out"
    if not check_JOB_DONE(calc, output_path):
        brilloin_zone_path = dict_brilloin_zone_path[dict_phase[structure_name]]
        k_point_divisions = [20 for _ in range(len(brilloin_zone_path))]
        make_input(
            calc,
            input_path,
            prefix,
            pseudo_dir=pseudo_dir,
            params_structure=params_structure,
            brilloin_zone_path=brilloin_zone_path,
            k_point_divisions=k_point_divisions,
            nbnd=nbnd,
            ibrav=dict_ibrav[dict_phase[structure_name]],
            template_path=temp_scf,
        )
        command = f"mpirun -np 32 pw.x <{input_path} >{output_path}"
        ret = subprocess.run(command, shell=True, capture_output=True, text=True)
    if not check_JOB_DONE(calc, output_path):
        print("停止", flush=True)
        continue
    ################################################################################################################
    calc = "band_x"
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)
    input_path = f"{structure_dir}/bands/{structure_name}.band_x.in"
    output_path = f"{structure_dir}/bands/{structure_name}.band_x.out"
    if not check_JOB_DONE(calc, output_path):
        make_input(calc, input_path, prefix, template_path=temp_band_x)
        command = f"bands.x <{input_path} >{output_path}"
        ret = subprocess.run(command, shell=True, capture_output=True, text=True)
    if not check_JOB_DONE(calc, output_path):
        print("停止", flush=True)
        continue
    ################################################################################################################
    calc = "plot"
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)

    pdos_dir_path = f"{structure_dir}/nscf"
    color_dict = dict_color_dict[structure_name]

    gnu_path = f"{structure_dir}/bands/{structure_name}.gnu"
    brilloin_zone_path = dict_brilloin_zone_path[dict_phase[structure_name]]
    k_point_divisions = [20 for _ in range(len(brilloin_zone_path))]

    output_path = f"{structure_dir}/nscf/{structure_name}.dos"
    EFermi = get_EFermi(output_path)
    print(f"{structure_name} : EFermi = {EFermi}")
    try:
        nscf_out_path = f"{structure_dir}/nscf/{structure_name}.nscf.out"
        highest_occupied = get_highest_occupied(nscf_out_path)
        print(f"{structure_name} : highest_occupied = {highest_occupied}")
        savefig_path = f"{structure_dir}/{structure_name}.pdos_highest_occupied.png"
        plot_pdos (pdos_dir_path, highest_occupied=highest_occupied, savefig_path=savefig_path, is_save=True, color_dict=color_dict)
        savefig_path = f"{structure_dir}/{structure_name}.bands_highest_occupied.png"
        plot_band(gnu_path, k_point_divisions, brilloin_zone_path, highest_occupied=highest_occupied, is_save=True, savefig_path=savefig_path, ylim= [-5, 5])
    except:
        pass
    savefig_path = f"{structure_dir}/{structure_name}.pdos_efermi.png"
    plot_pdos (pdos_dir_path, EFermi=EFermi, savefig_path=savefig_path, is_save=True, color_dict=color_dict)
    savefig_path = f"{structure_dir}/{structure_name}.bands_efermi.png"
    plot_band(gnu_path, k_point_divisions, brilloin_zone_path, EFermi=EFermi, is_save=True, savefig_path=savefig_path, ylim=[-5, 5])
print("ALL DONE", flush=True)
