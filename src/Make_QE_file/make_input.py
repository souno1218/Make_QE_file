import numpy as np
import pandas as pd
import re, warnings, os
from pathlib import Path
from .utils import check_args, round_half
import importlib.resources as pkg_resources


def make_input(
    calc,
    save_path,
    prefix,
    template_path=None,
    pseudo_dir=None,  # only use in scf,nscf,relax,vc-relax,bands
    params_structure=None,  # only use in scf,nscf,relax,vc-relax,bands
    min_kpoint=None,  # only use in scf,nscf,relax,vc-relax
    K_point_Density_Product=None,  # only use in scf,nscf,relax,vc-relax
    fix=False,  # only use in scf,nscf,relax,vc-relax
    nbnd=None,  # only use in nscf, bands
    brilloin_zone_path=None,  # only use in bands
    k_point_divisions=None,  # only use in bands
    ibrav=None,  # only use in bands, Some versions　of QE may require ibrav
    emax=None,  # only use in projwfc and dos
    emin=None,  # only use in projwfc and dos
    deltae=None,  # only use in projwfc and dos
    # emax, emin, deltae
    # If provided, these values override; otherwise, values from the template file are used.
    # If not found in the template, default values (see below) are applied.
    # emax, emin, deltae = 50, -50, 0.01
):
    check_args(False, calc=calc, save_path=save_path, prefix=prefix)
    match calc:
        case t if t in ["scf", "nscf", "relax", "vc-relax", "bands"]:
            check_args(False, pseudo_dir=pseudo_dir, params_structure=params_structure)
            check_args(
                True,
                template_path=template_path,
                min_kpoint=min_kpoint,
                K_point_Density_Product=K_point_Density_Product,
                brilloin_zone_path=brilloin_zone_path,
                k_point_divisions=k_point_divisions,
                nbnd=nbnd,
                ibrav=ibrav,
            )
            sub_make_input(
                calc,
                save_path,
                prefix,
                pseudo_dir,
                params_structure,
                template_path=template_path,
                min_kpoint=min_kpoint,  # only use in scf,nscf,relax,vc-relax
                K_point_Density_Product=K_point_Density_Product,  # only use in scf,nscf,relax,vc-relax
                fix=fix,  # only use in scf,nscf,relax,vc-relax
                brilloin_zone_path=brilloin_zone_path,  # only use in bands
                k_point_divisions=k_point_divisions,  # only use in bands
                nbnd=nbnd,  # only use in nscf, bands
                ibrav=ibrav,  # only use in bands, Some versions　of QE may require ibrav
            )
        case "band_x":
            check_args(True, template_path=template_path)
            make_input_band_x(save_path, prefix, template_path=template_path)
        case "projwfc":
            check_args(True, emax=emax, emin=emin, deltae=deltae, template_path=template_path)
            make_input_projwfc(save_path, prefix, emax=emax, emin=emin, deltae=deltae, template_path=template_path)
        case "dos":
            check_args(True, emax=emax, emin=emin, deltae=deltae, template_path=template_path)
            make_input_dos(save_path, prefix, emax=emax, emin=emin, deltae=deltae, template_path=template_path)
        case _:
            raise Exception(f"not support this calc : {calc}")  #'md','vc-md' not support


def sub_make_input(
    calc,
    save_path,
    prefix,
    pseudo_dir,
    params_structure,
    template_path=None,
    min_kpoint=None,  # only use in scf,nscf,relax,vc-relax
    K_point_Density_Product=None,  # only use in scf,nscf,relax,vc-relax
    fix=False,  # only use in scf,nscf,relax,vc-relax
    brilloin_zone_path=None,  # only use in bands
    k_point_divisions=None,  # only use in bands
    nbnd=None,  # only use in nscf, bands
    ibrav=None,  # only use in bands, Some versions　of QE may require ibrav
):
    int_nan = -100
    df_ATOMIC_POSITIONS = params_structure["df_ATOMIC_POSITIONS"]
    df_pseudo_path_num_ecut = make_pseudo_path_num_ecut(pseudo_dir, df_ATOMIC_POSITIONS)
    if template_path is None:
        templates_path = str(pkg_resources.files("Make_QE_file").joinpath("templates"))
        template_path = f"{templates_path}/template.in"
    with open(template_path, "r") as f:
        template = f.readlines()
    base = template
    arr_n_flont_space = np.array([re.match(r" *", i).end() for i in base])
    n_flont_space = int(np.mean(arr_n_flont_space[arr_n_flont_space != 0]))
    list_title = [
        "&CONTROL",
        "&SYSTEM",
        "&ELECTRONS",
        "&IONS",
        "&CELL",
        "&FCP",
        "&RISM",
        "ATOMIC_SPECIES",
        "ATOMIC_POSITIONS",
        "K_POINTS",
        "ADDITIONAL_K_POINTS",
        "CELL_PARAMETERS",
        "CONSTRAINTS",
        "OCCUPATIONS",
        "ATOMIC_VELOCITIES",
        "ATOMIC_FORCES",
        "SOLVENTS",
        "HUBBARD",
    ]
    index_dict = {
        "&CONTROL": np.full((2), int_nan, dtype="int64"),
        "&SYSTEM": np.full((2), int_nan, dtype="int64"),
        "ATOMIC_SPECIES": np.full((2), int_nan, dtype="int64"),
        "ATOMIC_POSITIONS": np.full((2), int_nan, dtype="int64"),
        "K_POINTS": np.full((2), int_nan, dtype="int64"),
    }
    last_title = None
    for i in range(len(base)):
        line_title = None
        for j in list_title:
            if j in base[i].upper():
                if not "=" in base[i]:
                    line_title = j
        if line_title in index_dict.keys():
            index_dict[line_title][0] = i
            if last_title in index_dict.keys():
                index_dict[last_title][1] = i
            last_title = line_title
        elif (not line_title is None) and (last_title in index_dict.keys()):
            index_dict[last_title][1] = i
            last_title = line_title
    if last_title in index_dict.keys():
        index_dict[last_title][1] = len(base)
    # &CONTROL
    start_index = index_dict["&CONTROL"][0]
    last_index = index_dict["&CONTROL"][1]
    txt = " " * n_flont_space + f"calculation = {calc},\n"
    base.insert(start_index + 1, txt)
    txt = " " * n_flont_space + f"prefix = {prefix},\n"
    base.insert(start_index + 2, txt)
    txt = " " * n_flont_space + f"pseudo_dir = '{pseudo_dir}',\n"
    base.insert(start_index + 3, txt)
    n = 3
    i = start_index + 4
    drop_list = ["calculation", "prefix", "pseudo_dir"]
    while True:
        drop = False
        for j in drop_list:
            if j in base[i]:
                drop = True
                break
        if drop:
            base.pop(i)
            n -= 1
        else:
            i += 1
        if i == last_index + n:
            break
    for i in index_dict.keys():
        if last_index <= index_dict[i][0]:
            index_dict[i][:] += n

    # &SYSTEM
    start_index = index_dict["&SYSTEM"][0]
    last_index = index_dict["&SYSTEM"][1]
    txt = " " * n_flont_space + f"space_group = {params_structure['space_group_number']},\n"
    base.insert(start_index + 1, txt)
    A = "{:.05f}".format(round_half(params_structure["a"]))
    txt = " " * n_flont_space + f"A     = {A},\n"
    base.insert(start_index + 2, txt)
    B = "{:.05f}".format(round_half(params_structure["b"]))
    txt = " " * n_flont_space + f"B     = {B},\n"
    base.insert(start_index + 3, txt)
    C = "{:.05f}".format(round_half(params_structure["c"]))
    txt = " " * n_flont_space + f"C     = {C},\n"
    base.insert(start_index + 4, txt)
    cos_alpha = "{:.05f}".format(round_half(np.cos(np.pi * params_structure["alpha"] / 180)))
    txt = " " * n_flont_space + f"cosBC = {cos_alpha},\n"
    base.insert(start_index + 5, txt)
    cos_beta = "{:.05f}".format(round_half(np.cos(np.pi * params_structure["beta"] / 180)))
    txt = " " * n_flont_space + f"cosAC = {cos_beta},\n"
    base.insert(start_index + 6, txt)
    cos_gamma = "{:.05f}".format(round_half(np.cos(np.pi * params_structure["gamma"] / 180)))
    txt = " " * n_flont_space + f"cosAB = {cos_gamma},\n"
    base.insert(start_index + 7, txt)
    nat = df_ATOMIC_POSITIONS.shape[0]
    txt = " " * n_flont_space + f"nat = {nat},\n"
    base.insert(start_index + 8, txt)
    ntyp = np.unique(df_ATOMIC_POSITIONS["symbol"].values).shape[0]
    txt = " " * n_flont_space + f"ntyp = {ntyp},\n"
    base.insert(start_index + 9, txt)
    ecutwfc = np.max(df_pseudo_path_num_ecut["ecutwfc"])
    txt = " " * n_flont_space + f"ecutwfc = {ecutwfc},\n"
    base.insert(start_index + 10, txt)
    ecutrho = np.max(df_pseudo_path_num_ecut["ecutrho"])
    txt = " " * n_flont_space + f"ecutrho = {ecutrho},\n"
    base.insert(start_index + 11, txt)
    if calc in ["bands", "nscf"]:
        txt = " " * n_flont_space + f"nbnd = {nbnd},\n"
        base.insert(start_index + 12, txt)
        n = 12
        if calc == "bands":
            if not ibrav is None:  # only use in bands
                txt = " " * n_flont_space + f"ibrav = {ibrav},\n"
                base.insert(start_index + 13, txt)
                n = 13
            else:
                warnings.warn("Some versions of QE may require ibrav in Calculation : bands")
    else:
        n = 11
    i = start_index + n + 1
    drop_list = [
        "A ",
        "A=",
        "B ",
        "B=",
        "C ",
        "C=",
        "cosBC",
        "cosAC",
        "cosAB",
        "celldm(1)",
        "celldm(2)",
        "celldm(3)",
        "celldm(4)",
        "celldm(5)",
        "celldm(6)",
        "space_group",
        "ibrav",
        "nat",
        "ntyp",
        "ecutwfc",
        "ecutrho",
        "nbnd",
    ]
    while True:
        drop = False
        for j in drop_list:
            if j.upper() in base[i].upper():
                drop = True
                break
        if drop:
            base.pop(i)
            n -= 1
        else:
            i += 1
        if i == last_index + n:
            break
    for i in index_dict.keys():
        if last_index <= index_dict[i][0]:
            index_dict[i][:] += n

    # ATOMIC_SPECIES
    start_index = index_dict["ATOMIC_SPECIES"][0]
    last_index = index_dict["ATOMIC_SPECIES"][1]
    len_atom = max(len(i) for i in df_pseudo_path_num_ecut["atom"])
    n = 0
    for i in range(last_index - 1, start_index, -1):
        base.pop(i)
        n -= 1
    for i in range(df_pseudo_path_num_ecut.shape[0]):
        atom = df_pseudo_path_num_ecut["atom"].iloc[i]
        file_name = df_pseudo_path_num_ecut["file_name"].iloc[i]
        txt = " " * n_flont_space + f"{atom.rjust(len_atom)}   -1   {file_name}\n"
        base.insert(start_index + 1 + i, txt)
        n += 1
    for i in index_dict.keys():
        if last_index <= index_dict[i][0]:
            index_dict[i][:] += n
    # ATOMIC_POSITIONS
    start_index = index_dict["ATOMIC_POSITIONS"][0]
    last_index = index_dict["ATOMIC_POSITIONS"][1]
    base[start_index] = "ATOMIC_POSITIONS {crystal_sg}\n"
    n = 0
    for i in range(last_index - 1, start_index, -1):
        base.pop(i)
        n -= 1
    len_atom = max(len(i) for i in df_ATOMIC_POSITIONS["symbol"])
    len_x = max(len(i) for i in df_ATOMIC_POSITIONS["str_x"])
    len_y = max(len(i) for i in df_ATOMIC_POSITIONS["str_y"])
    len_z = max(len(i) for i in df_ATOMIC_POSITIONS["str_z"])
    for i in range(df_ATOMIC_POSITIONS.shape[0]):
        label = df_ATOMIC_POSITIONS["label"].iloc[i]
        symbol = df_ATOMIC_POSITIONS["symbol"][label]
        str_xyz = df_ATOMIC_POSITIONS[["str_x", "str_y", "str_z"]].loc[label].values
        xyz = str_xyz.astype("float")
        txt = " " * n_flont_space + f"{symbol.rjust(len_atom)}   "
        txt += f"{str_xyz[0].rjust(len_x)} {str_xyz[1].rjust(len_y)} {str_xyz[2].rjust(len_z)}"
        if fix:
            is_move = ((xyz * 2**5) % 1 != 0) * 1
            txt += f"   {is_move[0]} {is_move[1]} {is_move[2]}"
        txt += "\n"
        base.insert(start_index + 1 + i, txt)
        n += 1
    for i in index_dict.keys():
        if last_index <= index_dict[i][0]:
            index_dict[i][:] += n
    # K_POINTS
    start_index = index_dict["K_POINTS"][0]
    last_index = index_dict["K_POINTS"][1]
    n = 0
    for i in range(last_index - 1, start_index, -1):
        base.pop(i)
        n -= 1
    if calc != "bands":
        base[start_index] = "K_POINTS {automatic}\n"
    else:
        base[start_index] = "K_POINTS {tpiba_b}\n"
    lattice = np.array([params_structure["a"], params_structure["b"], params_structure["c"]])
    k_points = make_k_points(
        lattice,
        calc,
        min_kpoint=min_kpoint,
        K_point_Density_Product=K_point_Density_Product,
        brilloin_zone_path=brilloin_zone_path,
        k_point_divisions=k_point_divisions,
    )
    txt = " " * n_flont_space + k_points
    base.insert(start_index + 1 + i, txt)
    n += 1
    for i in index_dict.keys():
        if last_index <= index_dict[i][0]:
            index_dict[i][:] += n
    with open(save_path, "wb") as input:
        input.write("".join(base).encode("utf-8"))
    return base


def make_input_projwfc(save_path, prefix, emax=None, emin=None, deltae=None, template_path=None):
    # https://www.quantum-espresso.org/Doc/INPUT_DOS.html
    if template_path is None:
        templates_path = str(pkg_resources.files("Make_QE_file").joinpath("templates"))
        template_path = f"{templates_path}/template.projwfc.in"
    with open(template_path, "r") as f:
        template = f.readlines()
    base = template
    arr_n_flont_space = np.array([re.match(r" *", i).end() for i in base])
    n_flont_space = int(np.mean(arr_n_flont_space[arr_n_flont_space != 0]))
    start_index = 0
    for i in range(len(base)):
        if "&PROJWFC" in base[i].upper():
            if not "=" in base[i]:
                start_index = i
                break
    drop_list = ["PREFIX"]
    for i in range(len(base) - 1, start_index, -1):
        for j in drop_list:
            if j in base[i].upper():
                base.pop(i)

    check_exist_emax, check_exist_emin, check_exist_deltae = False, False, False
    for i in range(len(base) - 1, start_index, -1):
        if "EMAX" in base[i].upper():
            check_exist_emax = True
            if not emax is None:
                base.pop(i)
        if "EMIN" in base[i].upper():
            check_exist_emin = True
            if not emin is None:
                base.pop(i)
        if "DELTAE" in base[i].upper():
            check_exist_deltae = True
            if not deltae is None:
                base.pop(i)

    txt = " " * n_flont_space + f"prefix='{prefix}',\n"
    base.insert(start_index + 1, txt)
    n = 2
    if not emax is None:
        txt = " " * n_flont_space + f"emax = {emax},\n"
        base.insert(start_index + n, txt)
        n += 1
    elif not check_exist_emax:
        emax = 50.0
        txt = " " * n_flont_space + f"emax = {emax},\n"
        base.insert(start_index + n, txt)
        n += 1

    if not emin is None:
        txt = " " * n_flont_space + f"emin = {emin},\n"
        base.insert(start_index + n, txt)
        n += 1
    elif not check_exist_emin:
        emin = -50.0
        txt = " " * n_flont_space + f"emin = {emin},\n"
        base.insert(start_index + n, txt)
        n += 1

    if not deltae is None:
        txt = " " * n_flont_space + f"deltae = {deltae},\n"
        base.insert(start_index + n, txt)
        n += 1
    elif not check_exist_deltae:
        deltae = 0.01
        txt = " " * n_flont_space + f"deltae = {deltae},\n"
        base.insert(start_index + n, txt)
        n += 1

    with open(save_path, "wb") as input:
        input.write("".join(base).encode("utf-8"))
    # print("make_input_projwfc Done", flush=True)
    # print(save_path, flush=True)


def make_input_dos(save_path, prefix, emax=None, emin=None, deltae=None, template_path=None):
    # https://www.quantum-espresso.org/Doc/INPUT_DOS.html
    if template_path is None:
        templates_path = str(pkg_resources.files("Make_QE_file").joinpath("templates"))
        template_path = f"{templates_path}/template.dos.in"
    with open(template_path, "r") as f:
        template = f.readlines()
    base = template
    arr_n_flont_space = np.array([re.match(r" *", i).end() for i in base])
    n_flont_space = int(np.mean(arr_n_flont_space[arr_n_flont_space != 0]))
    start_index = 0
    for i in range(len(base)):
        if "&DOS" in base[i].upper():
            if not "=" in base[i]:
                start_index = i
                break
    drop_list = ["PREFIX", "FILDOS"]
    for i in range(len(base) - 1, start_index, -1):
        for j in drop_list:
            if j in base[i].upper():
                base.pop(i)

    check_exist_emax, check_exist_emin, check_exist_deltae = False, False, False
    for i in range(len(base) - 1, start_index, -1):
        if "EMAX" in base[i].upper():
            check_exist_emax = True
            if not emax is None:
                base.pop(i)
        if "EMIN" in base[i].upper():
            check_exist_emin = True
            if not emin is None:
                base.pop(i)
        if "DELTAE" in base[i].upper():
            check_exist_deltae = True
            if not deltae is None:
                base.pop(i)

    txt = " " * n_flont_space + f"prefix='{prefix}',\n"
    base.insert(start_index + 1, txt)
    txt = " " * n_flont_space + f"fildos='{prefix}.dos',\n"
    base.insert(start_index + 2, txt)
    n = 3
    if not emax is None:
        txt = " " * n_flont_space + f"emax = {emax},\n"
        base.insert(start_index + n, txt)
        n += 1
    elif not check_exist_emax:
        emax = 50.0
        txt = " " * n_flont_space + f"emax = {emax},\n"
        base.insert(start_index + n, txt)
        n += 1

    if not emin is None:
        txt = " " * n_flont_space + f"emin = {emin},\n"
        base.insert(start_index + n, txt)
        n += 1
    elif not check_exist_emin:
        emin = -50.0
        txt = " " * n_flont_space + f"emin = {emin},\n"
        base.insert(start_index + n, txt)
        n += 1

    if not deltae is None:
        txt = " " * n_flont_space + f"deltae = {deltae},\n"
        base.insert(start_index + n, txt)
        n += 1
    elif not check_exist_deltae:
        deltae = 0.01
        txt = " " * n_flont_space + f"deltae = {deltae},\n"
        base.insert(start_index + n, txt)
        n += 1

    with open(save_path, "wb") as input:
        input.write("".join(base).encode("utf-8"))
    # print("make_input_dos Done", flush=True)
    # print(save_path, flush=True)


def make_input_band_x(save_path, prefix, template_path=None):
    # https://www.quantum-espresso.org/Doc/INPUT_BANDS.html
    if template_path is None:
        templates_path = str(pkg_resources.files("Make_QE_file").joinpath("templates"))
        template_path = f"{templates_path}/template.band_x.in"
    with open(template_path, "r") as f:
        template = f.readlines()
    base = template
    arr_n_flont_space = np.array([re.match(r" *", i).end() for i in base])
    n_flont_space = int(np.mean(arr_n_flont_space[arr_n_flont_space != 0]))
    start_index = 0
    for i in range(len(base)):
        if "&BANDS" in base[i].upper():
            if not "=" in base[i]:
                start_index = i
                break
    drop_list = ["PREFIX", "FILBAND"]
    for i in range(len(base) - 1, start_index, -1):
        for j in drop_list:
            if j in base[i].upper():
                base.pop(i)
    txt = " " * n_flont_space + f"prefix='{prefix}',\n"
    base.insert(start_index + 1, txt)
    txt = " " * n_flont_space + f"filband = {prefix},\n"
    base.insert(start_index + 2, txt)
    with open(save_path, "wb") as input:
        input.write("".join(base).encode("utf-8"))
    # print("make_input_bands Done", flush=True)
    # print(save_path, flush=True)


def make_pseudo_path_num_ecut(pseudo_dir, df_ATOMIC_POSITIONS):
    df = pd.DataFrame(columns=["atom", "file_name", "ecutwfc", "ecutrho"])
    UPF_FOLDER = Path(pseudo_dir)
    for i in np.unique(df_ATOMIC_POSITIONS["symbol"].values):
        list_info = [i]
        list_UPF = list(UPF_FOLDER.glob(f"{i}.*.UPF"))
        if len(list_UPF) == 0:
            raise Exception(f"UPF not found \n{i} has no UPF in {pseudo_dir}")
        file_name = str(list_UPF[0].name)
        list_info.append(file_name)
        path = pseudo_dir + "/" + file_name
        with open(path, "r") as f:
            data = f.readlines()
            for j in data:
                if "cutoff for wavefunctions:" in j:
                    for k in j.split():
                        try:
                            list_info.append(float(k))
                            break
                        except:
                            None
                    else:
                        raise Exception(f"cutoff_wavefunctions not found\nElement:{i} , UPF_path:{file_name}")
                elif "cutoff for charge density:" in j:
                    for k in j.split():
                        try:
                            list_info.append(float(k))
                            break
                        except:
                            None
                    else:
                        raise Exception(f"cutoff_charge_density not found\nElement:{i} , UPF_path:{file_name}")
        df.loc[i] = list_info
    return df


def make_k_points(
    lattice, calc, min_kpoint=None, K_point_Density_Product=None, brilloin_zone_path=None, k_point_divisions=None
):
    # lattice = np.array([a, b, c])
    if not calc in ["relax", "vc-relax", "scf", "nscf", "bands"]:  #'md','vc-md' not support
        raise Exception(f"not support this calc : {calc}")
    match calc:
        case t if t in ["scf", "nscf", "relax", "vc-relax"]:
            if K_point_Density_Product is None:
                K_point_Density_Product = 20  # angs
            if min_kpoint is None:
                if calc == "nscf":
                    min_kpoint = 4
                else:
                    min_kpoint = 2
            k_points = np.ceil(K_point_Density_Product / lattice).astype("int")
            k_points[k_points < min_kpoint] = min_kpoint
            k_points_text = ""
            for i in k_points:
                k_points_text += f"{i} "
            k_points_text += "0 0 0\n"
            return k_points_text
        case "bands":
            is_None_brilloin_zone_path = brilloin_zone_path is None
            is_None_k_point_divisions = k_point_divisions is None
            if is_None_brilloin_zone_path or is_None_k_point_divisions:
                E_txt = ""
                if is_None_brilloin_zone_path:
                    E_txt += "brilloin_zone_path not found"
                if is_None_k_point_divisions:
                    E_txt += "k_point_divisions not found"
                raise Exception(E_txt)
            if len(brilloin_zone_path) != len(k_point_divisions):
                raise Exception(f"len(brilloin_zone_path) and len(k_point_divisions) are different")
            len_max_brilloin_zone_path = max([len(i) for i in brilloin_zone_path])
            k_points_text = ""
            k_points_text += f"{int(len(k_point_divisions))}\n"
            for i, j in zip(brilloin_zone_path, k_point_divisions):
                k_points_text += f"{str(i).ljust(len_max_brilloin_zone_path)}     {str(int(j))}\n"
            return k_points_text
