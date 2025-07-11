from Make_QE_file import make_input, cif_to_params, check_JOB_DONE, plot_relax_out
from Make_QE_file import output_to_params, write_new_cif, output_to_nbnd, MonitorPlotRelaxOut
from Make_QE_file import get_EFermi, get_highest_occupied, plot_band, plot_pdos
import subprocess, os, datetime, random
from decimal import Decimal, ROUND_HALF_UP
import pandas as pd
import numpy as np

# pip uninstall Make_QE_file -y
# pip install git+https://github.com/souno1218/Make_QE_file.git
# pip install git+https://github.com/souno1218/Make_QE_file.git@feature

PGID = os.getpgid(0)


# check "need change"
# ================== ↓ need change ↓ ==================

df = pd.read_excel("/hpc/QE_dir/all_data.xlsx")
list_structure_name = list(df["name"].values)

dict_phase = {list_structure_name[i]: df["phase"].values[i] for i in range(len(list_structure_name))}

dict_brilloin_zone_path = {
    212: ["gG", "X", "gS", "gG", "Z"],  # suggest : Γ—X—P—N—Γ—M—S|S0—Γ|X—R|G—M
    325: ["gG", "X", "gS", "gG", "Z"],  # suggest : Γ—X—P—N—Γ—M—S|S0—Γ|X—R|G—M
    426: ["gG", "X", "M", "gG", "Z"],  # suggest : Γ—X—M—Γ—Z—R—A—Z|X—R|M—A
}

dict_color_dict = {}
for i in range(len(list_structure_name)):
    dict_color_dict[list_structure_name[i]] = {
        df.loc[i, "M"]: "r",
        df.loc[i, "X"]: "b",
        df.loc[i, "A"]: "y",
        df.loc[i, "B"]: "g",
        df.loc[i, "O"]: "c",
    }

random.seed(314)
random.shuffle(list_structure_name)

dict_ibrav = {212: 7, 325: 7, 426: 6}

dict_K_point_Density_Product = {212: [20, 20], 325: [10, 20], 426: [20, 20]}
min_kpoint = 1
fix = True
figsize = (12, 4)
interval = 60  # seconds
is_show = False
save_fig = True


pseudo_dir = "/hpc/QE_dir/UPF"
temp_scf = "/hpc/QE_dir/templates/template_smearing.in"
temp_projwfc = "/hpc/QE_dir/templates/template.projwfc.in"
temp_dos = "/hpc/QE_dir/templates/template.dos.in"
temp_band_x = "/hpc/QE_dir/templates/template.band_x.in"

# ================== ↑ need change ↑ ==================

date = str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")).center(
    max([len(i) for i in list_structure_name]), "="
)
print(f"=================={date}==================", flush=True)
print(f"PGID = {PGID}")
for structure_name in list_structure_name:
    str_structure_name = structure_name.center(max([len(i) for i in list_structure_name]), "=")
    print(f"=================={str_structure_name}==================", flush=True)

    structure_dir = f"/hpc/QE_dir/calc/{structure_name}"  # need change
    if not os.path.exists(structure_dir):
        os.makedirs(structure_dir)

    import_cif_path = f"/hpc/QE_dir/old_cif_dir/{structure_name}.cif"  # need change
    params_structure = cif_to_params(import_cif_path)
    base_params_structure = params_structure.copy()
    prefix = structure_name
    ################################################################################################################
    if not os.path.exists(f"{structure_dir}/relax"):
        os.makedirs(f"{structure_dir}/relax")
    os.chdir(f"{structure_dir}/relax")
    ################################################################################################################
    calc = "relax"
    K_point_Density_Product = dict_K_point_Density_Product[dict_phase[structure_name]][0]
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)
    input_path = f"{structure_dir}/relax/{structure_name}.{calc}.in"
    output_path = f"{structure_dir}/relax/{structure_name}.{calc}.out"
    if save_fig:
        savefig_path = f"{structure_dir}/relax/{structure_name}.reduce.{calc}.png"
    else:
        savefig_path = None
    if check_JOB_DONE(calc, output_path):
        plot_relax_out(
            output_path, title=f"{structure_name}_{calc}", savefig_path=savefig_path, figsize=figsize, is_show=is_show
        )
    else:
        make_input(
            calc,
            input_path,
            prefix,
            pseudo_dir=pseudo_dir,
            params_structure=params_structure,
            K_point_Density_Product=K_point_Density_Product,
            min_kpoint=min_kpoint,
            fix=fix,
            template_path=temp_scf,
        )
        with MonitorPlotRelaxOut(
            output_path,
            savefig_path=savefig_path,
            is_show=is_show,
            figsize=figsize,
            interval=interval,
            title=f"{structure_name}_{calc}",
        ):
            command = f"nice -n 5 mpirun -np 32 pw.x <{input_path} >{output_path}"
            ret = subprocess.run(command, shell=True, capture_output=True, text=True, start_new_session=True)

    if check_JOB_DONE(calc, output_path):
        params_structure = output_to_params(calc, output_path, params_structure)
    else:
        print("停止", flush=True)
        continue
    ################################################################################################################
    calc = "print"
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)
    print("params_structure")
    base_ATOMIC_POSITIONS = base_params_structure["df_ATOMIC_POSITIONS"][["str_x", "str_y", "str_z"]].values.astype(
        "float"
    )
    diff = (
        params_structure["df_ATOMIC_POSITIONS"][["str_x", "str_y", "str_z"]].astype("float").copy()
        - base_ATOMIC_POSITIONS
    )
    over_plus_one = (diff >= 1).values
    index = np.arange(over_plus_one.shape[0])[np.any(over_plus_one, axis=1)]
    for i in index:
        for j in np.arange(3)[over_plus_one[i]]:
            change_index = params_structure["df_ATOMIC_POSITIONS"].index[i]
            change_col = ["str_x", "str_y", "str_z"][j]
            num = float(params_structure["df_ATOMIC_POSITIONS"].loc[change_index, change_col]) - 1
            params_structure["df_ATOMIC_POSITIONS"].loc[change_index, change_col] = f"{num:.5f}"
    diff = (
        params_structure["df_ATOMIC_POSITIONS"][["str_x", "str_y", "str_z"]].astype("float").copy()
        - base_ATOMIC_POSITIONS
    )
    cp_params_structure = params_structure["df_ATOMIC_POSITIONS"][["str_x", "str_y", "str_z"]].copy()
    cp_params_structure[["diff_x", "diff_y", "diff_z"]] = diff
    print(cp_params_structure[["str_x", "diff_x", "str_y", "diff_y", "str_z", "diff_z"]])
    ################################################################################################################
    calc = "vc-relax"
    K_point_Density_Product = dict_K_point_Density_Product[dict_phase[structure_name]][0]
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)
    input_path = f"{structure_dir}/relax/{structure_name}.{calc}.in"
    output_path = f"{structure_dir}/relax/{structure_name}.{calc}.out"
    if save_fig:
        savefig_path = f"{structure_dir}/relax/{structure_name}.reduce.{calc}.png"
    else:
        savefig_path = None
    if check_JOB_DONE(calc, output_path):
        plot_relax_out(
            output_path, title=f"{structure_name}_{calc}", savefig_path=savefig_path, figsize=figsize, is_show=is_show
        )
    else:
        make_input(
            calc,
            input_path,
            prefix,
            pseudo_dir=pseudo_dir,
            params_structure=params_structure,
            K_point_Density_Product=K_point_Density_Product,
            min_kpoint=min_kpoint,
            fix=fix,
            template_path=temp_scf,
        )
        with MonitorPlotRelaxOut(
            output_path,
            savefig_path=savefig_path,
            is_show=is_show,
            figsize=figsize,
            interval=interval,
            title=f"{structure_name}_{calc}",
        ):
            command = f"nice -n 5 mpirun -np 32 pw.x <{input_path} >{output_path}"
            ret = subprocess.run(command, shell=True, capture_output=True, text=True, start_new_session=True)

    if check_JOB_DONE(calc, output_path):
        params_structure = output_to_params(calc, output_path, params_structure)
    else:
        print("停止", flush=True)
        continue
    ################################################################################################################
    calc = "make_cif"
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)
    cif_out_path = f"{structure_dir}/relax/{structure_name}.cif"
    write_new_cif(cif_out_path, structure_name, params_structure)
    ################################################################################################################
    calc = "print"
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)
    round_half = lambda x: float(Decimal(str(x)).quantize(Decimal("1e-5"), ROUND_HALF_UP))
    for i in ["a", "b", "c"]:
        diff = round_half(params_structure[i] - base_params_structure[i])
        print(f"{i} : {params_structure[i]}, diff_{i} : {diff}, ")
    print("params_structure")
    base_ATOMIC_POSITIONS = base_params_structure["df_ATOMIC_POSITIONS"][["str_x", "str_y", "str_z"]].values.astype(
        "float"
    )
    diff = (
        params_structure["df_ATOMIC_POSITIONS"][["str_x", "str_y", "str_z"]].astype("float").copy()
        - base_ATOMIC_POSITIONS
    )
    over_plus_one = (diff >= 1).values
    index = np.arange(over_plus_one.shape[0])[np.any(over_plus_one, axis=1)]
    for i in index:
        for j in np.arange(3)[over_plus_one[i]]:
            change_index = params_structure["df_ATOMIC_POSITIONS"].index[i]
            change_col = ["str_x", "str_y", "str_z"][j]
            num = float(params_structure["df_ATOMIC_POSITIONS"].loc[change_index, change_col]) - 1
            params_structure["df_ATOMIC_POSITIONS"].loc[change_index, change_col] = f"{num:.5f}"
    diff = (
        params_structure["df_ATOMIC_POSITIONS"][["str_x", "str_y", "str_z"]].astype("float").copy()
        - base_ATOMIC_POSITIONS
    )
    cp_params_structure = params_structure["df_ATOMIC_POSITIONS"][["str_x", "str_y", "str_z"]].copy()
    cp_params_structure[["diff_x", "diff_y", "diff_z"]] = diff
    print(cp_params_structure[["str_x", "diff_x", "str_y", "diff_y", "str_z", "diff_z"]])
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
            K_point_Density_Product=K_point_Density_Product,
            min_kpoint=min_kpoint,
            template_path=temp_scf,
        )
        command = f"nice -n 5 mpirun -np 32 pw.x <{input_path} >{output_path}"
        ret = subprocess.run(command, shell=True, capture_output=True, text=True, start_new_session=True)
    if check_JOB_DONE(calc, output_path):
        nbnd = 2 * output_to_nbnd(output_path)
    else:
        print("停止", flush=True)
        continue
    ################################################################################################################
    calc = "nscf"
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)
    input_path = f"{structure_dir}/nscf/{structure_name}.{calc}.in"
    output_path = f"{structure_dir}/nscf/{structure_name}.{calc}.out"
    if not check_JOB_DONE(calc, output_path):
        make_input(
            calc,
            input_path,
            prefix,
            pseudo_dir=pseudo_dir,
            params_structure=params_structure,
            K_point_Density_Product=K_point_Density_Product,
            min_kpoint=min_kpoint,
            nbnd=nbnd,
            template_path=temp_scf,
        )
        command = f"nice -n 5 mpirun -np 32 pw.x <{input_path} >{output_path}"
        ret = subprocess.run(command, shell=True, capture_output=True, text=True, start_new_session=True)
    if not check_JOB_DONE(calc, output_path):
        print("停止", flush=True)
        continue
    ################################################################################################################
    calc = "projwfc"
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)
    input_path = f"{structure_dir}/nscf/{structure_name}.{calc}.in"
    output_path = f"{structure_dir}/nscf/{structure_name}.{calc}.out"
    if not check_JOB_DONE(calc, output_path):
        make_input(calc, input_path, prefix, template_path=temp_projwfc)
        command = f"nice -n 5 projwfc.x <{input_path} >{output_path}"
        ret = subprocess.run(command, shell=True, capture_output=True, text=True, start_new_session=True)
    if not check_JOB_DONE(calc, output_path):
        print("停止", flush=True)
        continue
    ################################################################################################################
    calc = "dos"
    str_calc = calc.center(max([len(i) for i in list_structure_name]), "-")
    print(f"------------------{str_calc}------------------", flush=True)
    input_path = f"{structure_dir}/nscf/{structure_name}.{calc}.in"
    output_path = f"{structure_dir}/nscf/{structure_name}.{calc}.out"
    if not check_JOB_DONE(calc, output_path):
        make_input(calc, input_path, prefix, template_path=temp_dos)
        command = f"nice -n 5 dos.x <{input_path} >{output_path}"
        ret = subprocess.run(command, shell=True, capture_output=True, text=True, start_new_session=True)
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
            K_point_Density_Product=K_point_Density_Product,
            min_kpoint=min_kpoint,
            template_path=temp_scf,
        )
        command = f"nice -n 5 mpirun -np 32 pw.x <{input_path} >{output_path}"
        ret = subprocess.run(command, shell=True, capture_output=True, text=True, start_new_session=True)
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
            min_kpoint=min_kpoint,
            nbnd=nbnd,
            ibrav=dict_ibrav[dict_phase[structure_name]],
            template_path=temp_scf,
        )
        command = f"nice -n 5 mpirun -np 32 pw.x <{input_path} >{output_path}"
        ret = subprocess.run(command, shell=True, capture_output=True, text=True, start_new_session=True)
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
        command = f"nice -n 5 bands.x <{input_path} >{output_path}"
        ret = subprocess.run(command, shell=True, capture_output=True, text=True, start_new_session=True)
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
        if save_fig:
            savefig_path = f"{structure_dir}/{structure_name}.pdos_highest_occupied.png"
        else:
            savefig_path = None
        plot_pdos(
            pdos_dir_path,
            highest_occupied=highest_occupied,
            savefig_path=savefig_path,
            title=f"{structure_name}_pdos_highest_occupied",
            is_save=True,
            color_dict=color_dict,
            is_show=is_show,
        )
        if save_fig:
            savefig_path = f"{structure_dir}/{structure_name}.bands_highest_occupied.png"
        else:
            savefig_path = None
        plot_band(
            gnu_path,
            k_point_divisions,
            brilloin_zone_path,
            highest_occupied=highest_occupied,
            is_show=is_show,
            title=f"{structure_name}_pdos_highest_occupied",
            is_save=True,
            savefig_path=savefig_path,
            ylim=[-5, 5],
        )
    except:
        pass
    if save_fig:
        savefig_path = f"{structure_dir}/{structure_name}.pdos_efermi.png"
    else:
        savefig_path = None
    plot_pdos(
        pdos_dir_path,
        EFermi=EFermi,
        savefig_path=savefig_path,
        is_show=is_show,
        title=f"{structure_name}_pdos_efermi",
        is_save=True,
        color_dict=color_dict,
    )
    if save_fig:
        savefig_path = f"{structure_dir}/{structure_name}.bands_efermi.png"
    else:
        savefig_path = None
    plot_band(
        gnu_path,
        k_point_divisions,
        brilloin_zone_path,
        EFermi=EFermi,
        is_show=is_show,
        title=f"{structure_name}_bands_efermi",
        is_save=True,
        savefig_path=savefig_path,
        ylim=[-5, 5],
    )

print("ALL DONE", flush=True)
