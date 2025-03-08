import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.symmetry.groups import SpaceGroup
import re, os, datetime, warnings
from pathlib import Path
from decimal import Decimal, ROUND_HALF_UP

# https://qiita.com/xa_member/items/a519c53b63df75a68894


def cif_to_params(import_cif_path):
    params_cif = dict()
    list_loop_index = []
    with open(import_cif_path, "r") as f:
        cif_data = f.readlines()
        for i in range(len(cif_data)):
            if "_cell_length_a" in cif_data[i]:
                for j in cif_data[i].split():
                    try:
                        params_cif["a"] = round_half(j)
                    except:
                        None
            if "_cell_length_b" in cif_data[i]:
                for j in cif_data[i].split():
                    try:
                        params_cif["b"] = round_half(j)
                    except:
                        None
            if "_cell_length_c" in cif_data[i]:
                for j in cif_data[i].split():
                    try:
                        params_cif["c"] = round_half(j)
                    except:
                        None
            if "_cell_angle_alpha" in cif_data[i]:
                for j in cif_data[i].split():
                    try:
                        params_cif["alpha"] = round_half(j)
                    except:
                        None
            if "_cell_angle_beta" in cif_data[i]:
                for j in cif_data[i].split():
                    try:
                        params_cif["beta"] = round_half(j)
                    except:
                        None
            if "_cell_angle_gamma" in cif_data[i]:
                for j in cif_data[i].split():
                    try:
                        params_cif["gamma"] = round_half(j)
                    except:
                        None
            if "_space_group_IT_number" in cif_data[i]:
                for j in cif_data[i].split():
                    try:
                        params_cif["space_group_number"] = int(j)
                    except:
                        None
            if "loop_" in cif_data[i]:
                list_loop_index.append(i)
        list_loop_index.append(len(cif_data))
        for i in range(len(list_loop_index) - 1):
            for j in range(list_loop_index[i], list_loop_index[i + 1]):
                if "_atom_site" in cif_data[j]:
                    loop_index = list_loop_index[i]
        label_ATOMIC_POSITIONS = []
        for i in range(loop_index + 1, len(cif_data)):
            if "_atom_site" in cif_data[i]:
                label_ATOMIC_POSITIONS.append(cif_data[i].split()[0])
        df_ATOMIC_POSITIONS = pd.DataFrame(columns=["label", "symbol", "str_x", "str_y", "str_z"])
        for i in range(loop_index + 1, len(cif_data)):
            if len(cif_data[i].split()) == len(label_ATOMIC_POSITIONS):
                one_data = cif_data[i].split()
                one_ATOMIC_POSITIONS = []
                index = label_ATOMIC_POSITIONS.index("_atom_site_label")
                one_ATOMIC_POSITIONS.append(one_data[index])
                index = label_ATOMIC_POSITIONS.index("_atom_site_type_symbol")
                one_ATOMIC_POSITIONS.append(one_data[index])
                index = label_ATOMIC_POSITIONS.index("_atom_site_fract_x")
                str_x = "{:.05f}".format(round_half(one_data[index]))
                one_ATOMIC_POSITIONS.append(str_x)
                index = label_ATOMIC_POSITIONS.index("_atom_site_fract_y")
                str_y = "{:.05f}".format(round_half(one_data[index]))
                one_ATOMIC_POSITIONS.append(str_y)
                index = label_ATOMIC_POSITIONS.index("_atom_site_fract_z")
                str_z = "{:.05f}".format(round_half(one_data[index]))
                one_ATOMIC_POSITIONS.append(str_z)
                df_ATOMIC_POSITIONS.loc[one_ATOMIC_POSITIONS[0]] = one_ATOMIC_POSITIONS
        params_cif["df_ATOMIC_POSITIONS"] = df_ATOMIC_POSITIONS
    # check
    key_check = ["a", "b", "c", "alpha", "beta", "gamma", "space_group_number", "df_ATOMIC_POSITIONS"]
    for i in key_check:
        if not i in params_cif.keys():
            raise Exception("not found param {i} in cif")
    check_params(params_cif)
    return params_cif


def check_params(params):
    key_check = ["a", "b", "c", "alpha", "beta", "gamma", "space_group_number", "df_ATOMIC_POSITIONS"]
    for i in key_check:
        if not i in params.keys():
            raise Exception("not found param {i} in this params")
    for i in ["a", "b", "c", "alpha", "beta", "gamma"]:
        if not isinstance(params[i], (int, float)):
            raise TypeError(f"{i} must be int or float")
    if not isinstance(params["space_group_number"], int):
        raise TypeError(f"space_group_number must be int")
    df_ATOMIC_POSITIONS = params["df_ATOMIC_POSITIONS"]
    if not isinstance(df_ATOMIC_POSITIONS, pd.core.frame.DataFrame):
        raise TypeError(f"df_ATOMIC_POSITIONS must be pd.core.frame.DataFrame")
    for i in ["label", "symbol", "str_x", "str_y", "str_z"]:
        if not i in df_ATOMIC_POSITIONS.columns:
            raise Exception("not found {i} in df_ATOMIC_POSITIONS.columns")
    try:
        df_ATOMIC_POSITIONS[["str_x", "str_y", "str_z"]].astype("float")
    except ValueError:
        raise ValueError("df_ATOMIC_POSITIONS[['str_x', 'str_y', 'str_z']] need str(float)")
    return True


def check_args(is_None_continue, **kwargs):
    for i in kwargs.keys():
        if is_None_continue:
            if kwargs[i] is None:
                continue
        match i:
            case "calc":
                if not isinstance(kwargs[i], str):
                    raise TypeError(f"The argument {i} must be str")
                support_calc = ["relax", "vc-relax", "scf", "nscf", "bands", "band_x", "projwfc", "dos"]
                if not kwargs[i] in support_calc:  #'md','vc-md' not support
                    raise TypeError(f"The argument {i} must be str")
            case "save_path":
                if not isinstance(kwargs[i], str):
                    raise TypeError(f"The argument {i} must be str")
            case "prefix":
                if not isinstance(kwargs[i], str):
                    raise TypeError(f"The argument {i} must be str")
            case "template_path":
                if not isinstance(kwargs[i], str):
                    raise TypeError(f"The argument {i} must be str")
            case "pseudo_dir":
                if not isinstance(kwargs[i], str):
                    raise TypeError(f"The argument {i} must be str")
            case "params_structure":
                if not isinstance(kwargs[i], pd.core.frame.DataFrame):
                    raise TypeError(f"The argument {i} must be pd.core.frame.DataFrame")
                check_params(kwargs[i])
            case "min_kpoint":
                if not isinstance(kwargs[i], int):
                    raise TypeError(f"The argument {i} must be str")
                if kwargs[i] < 0:
                    raise ValueError(f"The argument {i} must be str")
            case "k_fineness_Magnification":
                if not isinstance(kwargs[i], (int, float)):
                    raise TypeError(f"The argument {i} must be int or float")
            case "nbnd":
                if not isinstance(kwargs[i], int):
                    raise TypeError(f"The argument {i} must be str")
                if kwargs[i] < 1:
                    raise ValueError(f"The argument {i} must be str")
            case "brilloin_zone_path":
                try:
                    np.asarray(kwargs[i])
                except Exception:
                    raise TypeError(f"The argument {i} must be array_like(str)")
                for j in kwargs[i]:
                    if not isinstance(j, str):
                        raise TypeError(f"The argument {i} must be array_like(str)")
            case "k_point_divisions":
                try:
                    np.asarray(kwargs[i])
                except Exception:
                    raise TypeError(f"The argument {i} must be array_like(int)")
                for j in kwargs[i]:
                    if not isinstance(j, int):
                        raise TypeError(f"The argument {i} must be array_like(int)")
            case "ibrav":
                if not isinstance(kwargs[i], int):
                    raise TypeError(f"The argument {i} must be str")
                if kwargs[i] < 0:
                    raise ValueError(f"The argument {i} must be str")
            case "deltae":
                if not isinstance(kwargs[i], (int, float)):
                    raise TypeError(f"The argument {i} must be int or float")


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
    lattice, calc, min_kpoint=None, k_fineness_Magnification=None, brilloin_zone_path=None, k_point_divisions=None
):
    # lattice = np.array([a, b, c])
    if not calc in ["relax", "vc-relax", "scf", "nscf", "bands"]:  #'md','vc-md' not support
        raise Exception(f"not support this calc : {calc}")
    match calc:
        case t if t in ["scf", "nscf", "relax", "vc-relax"]:
            if k_fineness_Magnification is None:
                k_fineness_Magnification = 20  # angs
            if min_kpoint is None:
                if calc == "nscf":
                    min_kpoint = 4
                else:
                    min_kpoint = 2
            k_points = np.ceil(k_fineness_Magnification / lattice).astype("int")
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


def make_input(
    calc,
    save_path,
    prefix,
    template_path=None,
    pseudo_dir=None,  # only use in scf,nscf,relax,vc-relax,bands
    params_structure=None,  # only use in scf,nscf,relax,vc-relax,bands
    min_kpoint=None,  # only use in scf,nscf,relax,vc-relax
    k_fineness_Magnification=None,  # only use in scf,nscf,relax,vc-relax
    nbnd=None,  # only use in nscf, bands
    brilloin_zone_path=None,  # only use in bands
    k_point_divisions=None,  # only use in bands
    ibrav=None,  # only use in bands, Some versions　of QE may require ibrav
    deltae=None,  # only use in projwfc and dos, Some versions　of QE may require ibrav
):
    check_args(False, calc=calc, save_path=save_path, prefix=prefix)
    match calc:
        case t if t in ["scf", "nscf", "relax", "vc-relax", "bands"]:
            check_args(False, pseudo_dir=pseudo_dir, params_structure=params_structure)
            check_args(
                True,
                template_path=template_path,
                min_kpoint=min_kpoint,
                k_fineness_Magnification=k_fineness_Magnification,
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
                k_fineness_Magnification=k_fineness_Magnification,  # only use in scf,nscf,relax,vc-relax
                brilloin_zone_path=brilloin_zone_path,  # only use in bands
                k_point_divisions=k_point_divisions,  # only use in bands
                nbnd=nbnd,  # only use in nscf, bands
                ibrav=ibrav,  # only use in bands, Some versions　of QE may require ibrav
            )
        case "band_x":
            check_args(True, template_path=template_path)
            make_input_band_x(save_path, prefix, template_path=template_path)
        case "projwfc":
            check_args(True, deltae=deltae, template_path=template_path)
            make_input_projwfc(save_path, prefix, deltae=deltae, template_path=template_path)
        case "dos":
            check_args(True, deltae=deltae, template_path=template_path)
            make_input_dos(save_path, prefix, deltae=deltae, template_path=template_path)
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
    k_fineness_Magnification=None,  # only use in scf,nscf,relax,vc-relax
    brilloin_zone_path=None,  # only use in bands
    k_point_divisions=None,  # only use in bands
    nbnd=None,  # only use in nscf, bands
    ibrav=None,  # only use in bands, Some versions　of QE may require ibrav
):
    int_nan = -100
    df_ATOMIC_POSITIONS = params_structure["df_ATOMIC_POSITIONS"]
    df_pseudo_path_num_ecut = make_pseudo_path_num_ecut(pseudo_dir, df_ATOMIC_POSITIONS)
    if template_path is None:
        template_path = "templates/template.in"
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
        txt += f"{str_xyz[0].rjust(len_x)} {str_xyz[1].rjust(len_y)} {str_xyz[2].rjust(len_z)}   "
        is_move = ((xyz * 2**5) % 1 != 0) * 1
        txt += f"{is_move[0]} {is_move[1]} {is_move[2]}\n"
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
        k_fineness_Magnification=k_fineness_Magnification,
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


def make_input_projwfc(save_path, prefix, deltae=None, template_path=None):
    # https://www.quantum-espresso.org/Doc/INPUT_DOS.html
    if template_path is None:
        template_path = "templates/template.projwfc.in"
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
    drop_list = ["OUTDIR", "PREFIX", "DELTAE"]
    for i in range(len(base) - 1, start_index, -1):
        for j in drop_list:
            if j in base[i].upper():
                base.pop(i)
    txt = " " * n_flont_space + "outdir = './work/',\n"
    base.insert(start_index + 1, txt)
    txt = " " * n_flont_space + f"prefix='{prefix}',\n"
    base.insert(start_index + 2, txt)
    if deltae != None:
        deltae = 0.01
    txt = " " * n_flont_space + f"deltae = {deltae},\n"
    base.insert(start_index + 3, txt)
    with open(save_path, "wb") as input:
        input.write("".join(base).encode("utf-8"))
    # print("make_input_projwfc Done", flush=True)
    # print(save_path, flush=True)


def make_input_dos(save_path, prefix, deltae=None, template_path=None):
    # https://www.quantum-espresso.org/Doc/INPUT_DOS.html
    if template_path is None:
        template_path = "templates/template.projwfc.in"
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
    drop_list = ["OUTDIR", "PREFIX", "FILDOS"]
    for i in range(len(base) - 1, start_index, -1):
        for j in drop_list:
            if j in base[i].upper():
                base.pop(i)
    txt = " " * n_flont_space + "outdir = './work/',\n"
    base.insert(start_index + 1, txt)
    txt = " " * n_flont_space + f"prefix='{prefix}',\n"
    base.insert(start_index + 2, txt)
    txt = " " * n_flont_space + f"fildos='{prefix}.dos',\n"
    base.insert(start_index + 3, txt)
    if deltae != None:
        deltae = 0.01
    txt = " " * n_flont_space + f"deltae = {deltae},\n"
    base.insert(start_index + 4, txt)
    with open(save_path, "wb") as input:
        input.write("".join(base).encode("utf-8"))
    # print("make_input_dos Done", flush=True)
    # print(save_path, flush=True)


def make_input_band_x(save_path, prefix, template_path=None):
    # https://www.quantum-espresso.org/Doc/INPUT_BANDS.html
    if template_path is None:
        template_path = "templates/template.band_x.in"
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
    drop_list = ["OUTDIR", "PREFIX", "FILBAND"]
    for i in range(len(base) - 1, start_index, -1):
        for j in drop_list:
            if j in base[i].upper():
                base.pop(i)
    txt = " " * n_flont_space + "outdir = './work/',\n"
    base.insert(start_index + 1, txt)
    txt = " " * n_flont_space + f"prefix='{prefix}',\n"
    base.insert(start_index + 2, txt)
    txt = " " * n_flont_space + f"filband = {prefix},\n"
    base.insert(start_index + 3, txt)
    with open(save_path, "wb") as input:
        input.write("".join(base).encode("utf-8"))
    # print("make_input_bands Done", flush=True)
    # print(save_path, flush=True)


def output_to_params(calc, import_out_path, base_params):
    if calc == "relax":
        return relax_output_to_params(import_out_path, base_params)
    elif calc == "vc-relax":
        return vc_relax_output_to_params(import_out_path, base_params)
    else:
        raise Exception("not support this calc : {calc}")


def relax_output_to_params(import_out_path, base_params):
    int_nan = -100
    return_params = base_params.copy()
    with open(import_out_path, "r") as f:
        out_data = f.readlines()
        crystal_axes = None
        Cartesian_axes = None
        for i in range(len(out_data)):
            if "crystal axes: " in out_data[i]:
                crystal_axes = np.array(
                    [out_data[i + 1].split()[3:6], out_data[i + 2].split()[3:6], out_data[i + 3].split()[3:6]]
                ).astype("float64")
            if "Cartesian axes" in out_data[i]:
                if Cartesian_axes is None:
                    Cartesian_axes = pd.DataFrame(columns=["atom", "x", "y", "z"])
                    j = 3
                    while not "number of k points" in out_data[i + j + 1]:
                        one_line_split = out_data[i + j].split()
                        Cartesian_axes.loc[int(one_line_split[0])] = [
                            one_line_split[1],
                            one_line_split[6],
                            one_line_split[7],
                            one_line_split[8],
                        ]
                        j += 1
            if "End final coordinates" in out_data[i]:
                for j in range(i - 1, -1, -1):
                    if "ATOMIC_POSITIONS" in out_data[j]:
                        if "alat" in out_data[j]:
                            raise Exception("cannot convert alat data, sorry")
                        ATOMIC_POSITIONS = pd.DataFrame(columns=["atom", "x", "y", "z"])
                        for k in range(1, i - (j)):
                            one_line_split = out_data[j + k].split()
                            ATOMIC_POSITIONS.loc[k] = [
                                one_line_split[0],
                                one_line_split[1],
                                one_line_split[2],
                                one_line_split[3],
                            ]
                        break
    xyz_conv = np.dot(ATOMIC_POSITIONS[["x", "y", "z"]].values.astype("float"), crystal_axes)
    ATOMIC_POSITIONS[["x", "y", "z"]] = round_half(xyz_conv)
    Cartesian_axes["y"] = round_half(Cartesian_axes["y"].values.astype("float") / base_params["b"] * base_params["a"])
    Cartesian_axes["z"] = round_half(Cartesian_axes["z"].values.astype("float") / base_params["c"] * base_params["a"])
    ATOMIC_POSITIONS["x"] = round_half(ATOMIC_POSITIONS["x"].values.astype("float"))
    ATOMIC_POSITIONS["y"] = round_half(
        ATOMIC_POSITIONS["y"].values.astype("float") * (base_params["a"] / (base_params["b"]))
    )
    ATOMIC_POSITIONS["z"] = round_half(
        ATOMIC_POSITIONS["z"].values.astype("float") * (base_params["a"] / (base_params["c"]))
    )
    df_ATOMIC_POSITIONS = return_params["df_ATOMIC_POSITIONS"].copy()
    for i in range(df_ATOMIC_POSITIONS.shape[0]):
        for j in range(Cartesian_axes.shape[0]):
            target_label = df_ATOMIC_POSITIONS["label"].iloc[i]
            target_symbol = df_ATOMIC_POSITIONS["symbol"].iloc[i]
            target_xyz = df_ATOMIC_POSITIONS[["str_x", "str_y", "str_z"]].iloc[i].values.astype("float")
            if target_symbol == Cartesian_axes["atom"].iloc[j]:
                before_xyz = Cartesian_axes[["x", "y", "z"]].iloc[j].values.astype("float")
                after_xyz = ATOMIC_POSITIONS[["x", "y", "z"]].iloc[j].values.astype("float")
                if np.linalg.norm(target_xyz - before_xyz, ord=2) <= 1e-4:
                    df_ATOMIC_POSITIONS.loc[target_label, "str_x"] = "{:.05f}".format(round_half(after_xyz[0]))
                    df_ATOMIC_POSITIONS.loc[target_label, "str_y"] = "{:.05f}".format(round_half(after_xyz[1]))
                    df_ATOMIC_POSITIONS.loc[target_label, "str_z"] = "{:.05f}".format(round_half(after_xyz[2]))
                    break
                before_xyz = (before_xyz + 0.5) % 1
                after_xyz = (after_xyz + 0.5) % 1
                if np.linalg.norm(target_xyz - before_xyz, ord=2) <= 1e-4:
                    df_ATOMIC_POSITIONS.loc[target_label, "str_x"] = "{:.05f}".format(round_half(after_xyz[0]))
                    df_ATOMIC_POSITIONS.loc[target_label, "str_y"] = "{:.05f}".format(round_half(after_xyz[1]))
                    df_ATOMIC_POSITIONS.loc[target_label, "str_z"] = "{:.05f}".format(round_half(after_xyz[2]))
                    break
        else:
            print(target_label, target_symbol, target_xyz)
            raise
    return_params["df_ATOMIC_POSITIONS"] = df_ATOMIC_POSITIONS
    return return_params


def vc_relax_output_to_params(import_out_path, base_params):
    int_nan = -100
    return_params = base_params.copy()
    with open(import_out_path, "r") as f:
        out_data = f.readlines()
        crystal_axes = None
        reciprocal_axes = None
        Cartesian_axes = None
        CELL_PARAMETERS = None
        length_A = int_nan
        for i in range(len(out_data)):
            if "celldm(1)" in out_data[i]:
                line_celldm1 = out_data[i].split()
                for i in range(len(line_celldm1)):
                    if "celldm(1)" in line_celldm1[i]:
                        try:
                            if length_A == int_nan:
                                length_A = float(line_celldm1[i + 1])
                            elif not np.isclose(length_A, float(line_celldm1[i + 1])):
                                raise
                        except:
                            raise
            if "crystal axes: " in out_data[i]:
                if crystal_axes is None:
                    crystal_axes = np.array(
                        [out_data[i + 1].split()[3:6], out_data[i + 2].split()[3:6], out_data[i + 3].split()[3:6]]
                    ).astype("float64")
            if "reciprocal axes:" in out_data[i]:
                if reciprocal_axes is None:
                    reciprocal_axes = np.array(
                        [out_data[i + 1].split()[3:6], out_data[i + 2].split()[3:6], out_data[i + 3].split()[3:6]]
                    ).astype("float64")
            if "Cartesian axes" in out_data[i]:
                if Cartesian_axes is None:
                    Cartesian_axes = pd.DataFrame(columns=["atom", "x", "y", "z"])
                    j = 3
                    while not "number of k points" in out_data[i + j + 1]:
                        one_line_split = out_data[i + j].split()
                        Cartesian_axes.loc[int(one_line_split[0])] = [
                            one_line_split[1],
                            one_line_split[6],
                            one_line_split[7],
                            one_line_split[8],
                        ]
                        j += 1
            if "End final coordinates" in out_data[i]:
                for j in range(i - 1, -1, -1):
                    if "ATOMIC_POSITIONS" in out_data[j]:
                        last_ATOMIC_POSITIONS_index = j
                        if "alat" in out_data[j]:
                            raise Exception("cannot convert alat data, sorry")
                        ATOMIC_POSITIONS = pd.DataFrame(columns=["atom", "x", "y", "z"])
                        for k in range(1, i - j):
                            one_line_split = out_data[j + k].split()
                            ATOMIC_POSITIONS.loc[k] = [
                                one_line_split[0],
                                one_line_split[1],
                                one_line_split[2],
                                one_line_split[3],
                            ]
                        break
                for j in range(last_ATOMIC_POSITIONS_index - 1, -1, -1):
                    if "CELL_PARAMETERS" in out_data[j]:
                        split_ = out_data[j].split()
                        for k in range(len(split_)):
                            if "alat" in split_[k]:
                                new_A = float(re.search(r"-?\d+\.\d+", split_[k + 1]).group())
                                if not np.isclose(length_A, new_A):
                                    raise
                        CELL_PARAMETERS = np.array(
                            [out_data[j + 1].split(), out_data[j + 2].split(), out_data[j + 3].split()]
                        ).astype("float64")

    conventional_before2after = np.dot(CELL_PARAMETERS.T, reciprocal_axes)
    after_a = np.dot(conventional_before2after, [1, 0, 0])
    after_b = np.dot(conventional_before2after, [0, 1, 0])
    after_c = np.dot(conventional_before2after, [0, 0, 1])
    norm_after_a = np.linalg.norm(after_a, ord=2)
    norm_after_b = np.linalg.norm(after_b, ord=2)
    norm_after_c = np.linalg.norm(after_c, ord=2)
    return_params["a"] = round_half(base_params["a"] * norm_after_a)
    return_params["b"] = round_half(base_params["b"] * norm_after_b)
    return_params["c"] = round_half(base_params["c"] * norm_after_c)
    return_params["alpha"] = round_half(
        180 * np.arccos(np.inner(after_b, after_c) / (norm_after_b * norm_after_c)) / np.pi
    )
    return_params["beta"] = round_half(
        180 * np.arccos(np.inner(after_a, after_c) / (norm_after_a * norm_after_c)) / np.pi
    )
    return_params["gamma"] = round_half(
        180 * np.arccos(np.inner(after_a, after_b) / (norm_after_a * norm_after_b)) / np.pi
    )
    Cartesian_axes["y"] = round_half(Cartesian_axes["y"].values.astype("float") / base_params["b"] * base_params["a"])
    Cartesian_axes["z"] = round_half(Cartesian_axes["z"].values.astype("float") / base_params["c"] * base_params["a"])
    ATOMIC_POSITIONS["x"] = round_half(ATOMIC_POSITIONS["x"].values.astype("float") / norm_after_a)
    ATOMIC_POSITIONS["y"] = round_half(
        ATOMIC_POSITIONS["y"].values.astype("float") * (base_params["a"] / (base_params["b"] * norm_after_b))
    )
    ATOMIC_POSITIONS["z"] = round_half(
        ATOMIC_POSITIONS["z"].values.astype("float") * (base_params["a"] / (base_params["c"] * norm_after_c))
    )
    df_ATOMIC_POSITIONS = return_params["df_ATOMIC_POSITIONS"]
    for i in range(df_ATOMIC_POSITIONS.shape[0]):
        for j in range(Cartesian_axes.shape[0]):
            target_label = df_ATOMIC_POSITIONS["label"].iloc[i]
            target_symbol = df_ATOMIC_POSITIONS["symbol"].iloc[i]
            target_xyz = df_ATOMIC_POSITIONS[["str_x", "str_y", "str_z"]].iloc[i].values.astype("float")
            if target_symbol == Cartesian_axes["atom"].iloc[j]:
                before_xyz = Cartesian_axes[["x", "y", "z"]].iloc[j].values.astype("float")
                after_xyz = Cartesian_axes[["x", "y", "z"]].iloc[j].values.astype("float")
                if np.linalg.norm(target_xyz - before_xyz, ord=2) <= 1e-4:
                    df_ATOMIC_POSITIONS.loc[target_label, "str_x"] = "{:.05f}".format(round_half(after_xyz[0]))
                    df_ATOMIC_POSITIONS.loc[target_label, "str_y"] = "{:.05f}".format(round_half(after_xyz[1]))
                    df_ATOMIC_POSITIONS.loc[target_label, "str_z"] = "{:.05f}".format(round_half(after_xyz[2]))
                    break
                before_xyz = (before_xyz + 0.5) % 1
                after_xyz = (after_xyz + 0.5) % 1
                if np.linalg.norm(target_xyz - before_xyz, ord=2) <= 1e-4:
                    df_ATOMIC_POSITIONS.loc[target_label, "str_x"] = "{:.05f}".format(round_half(after_xyz[0]))
                    df_ATOMIC_POSITIONS.loc[target_label, "str_y"] = "{:.05f}".format(round_half(after_xyz[1]))
                    df_ATOMIC_POSITIONS.loc[target_label, "str_z"] = "{:.05f}".format(round_half(after_xyz[2]))
                    break
        else:
            raise
    return_params["df_ATOMIC_POSITIONS"] = df_ATOMIC_POSITIONS
    return return_params


def round_half(num_or_arr):
    if "ndarray" in str(type(num_or_arr)):
        return_arr = np.empty_like(num_or_arr)
        for i in range(len(num_or_arr.flat)):
            return_arr.flat[i] = float(Decimal(str(float(num_or_arr.flat[i]))).quantize(Decimal("1e-5"), ROUND_HALF_UP))
        return return_arr
    else:
        return float(Decimal(str(num_or_arr)).quantize(Decimal("1e-5"), ROUND_HALF_UP))


def output_to_nbnd(import_out_path):
    with open(import_out_path, "r") as f:
        out_data = f.readlines()
        for i in range(len(out_data)):
            if "number of Kohn-Sham states=" in out_data[i]:
                return int(out_data[i].split()[4])


def check_JOB_DONE(calc, QE_out_path):
    if os.path.isfile(QE_out_path):
        with open(QE_out_path, "r") as f:
            data = f.readlines()
        JOB_DONE = False
        if calc in ["vc-relax", "relax"]:
            End_final_coordinates = False
            for readline in data:
                if "End final coordinates" in readline:
                    End_final_coordinates = True
        for readline in data:
            if "JOB DONE." in readline:
                JOB_DONE = True
        if JOB_DONE:
            if calc in ["vc-relax", "relax"]:
                if End_final_coordinates:
                    print("JOB DONE., End final coordinates", flush=True)
                    return True
                else:
                    print("not End final coordinates", flush=True)
                    return False
            else:
                print("JOB DONE.", flush=True)
                return True
        else:  # not JOB_DONE:
            print("not JOB DONE.", flush=True)
            return False
    else:
        return False


def write_new_cif(out_path, material_name, params_structure):
    date = str(datetime.datetime.now().strftime("%Y-%m-%d"))
    time = str(datetime.datetime.now().strftime("%H:%M:%S"))
    space_group_number = params_structure["space_group_number"]
    hm_symbol = SpaceGroup.from_int_number(space_group_number).symbol  # Hermann-Mauguin 記号を取得
    with open(out_path, "w") as f:
        f.write("#======================================================================)\n")
        f.write("# CRYSTAL DATA\n")
        f.write("#----------------------------------------------------------------------\n")
        f.write("\n")
        f.write(f"data_{material_name}\n")
        f.write(f"_audit_creation_date    {date}\n")
        f.write(f"_audit_creation_time    {time}\n")
        # f.write(f"_audit_creation_method  \"{14:30:15}\"\n")
        f.write(f"_symmetry_space_group_number {space_group_number}\n")
        f.write(f"_space_group_name_H-M_alt '{hm_symbol}'\n")
        f.write(f"_cell_length_a     {params_structure['a']}\n")
        f.write(f"_cell_length_b     {params_structure['b']}\n")
        f.write(f"_cell_length_c     {params_structure['c']}\n")
        f.write(f"_cell_angle_alpha  {params_structure['alpha']}\n")
        f.write(f"_cell_angle_beta   {params_structure['beta']}\n")
        f.write(f"_cell_angle_gamma  {params_structure['gamma']}\n")
        f.write("\n")
        # In progress, unusable
        # f.write("loop_\n")
        # f.write("_symmetry_equiv_pos_as_xyz\n")
        # sg = SpaceGroup.from_int_number(params_structure["space_group_number"])
        # symmetry_xyz = [symmop_to_xyz(op) for op in sg.symmetry_ops]
        # for i in range(len(symmetry_xyz)):
        #     f.write(f"   '{symmetry_xyz[i]}'\n")
        f.write("\n")
        f.write("loop_\n")
        f.write("_atom_site_label\n")
        f.write("_atom_site_type_symbol\n")
        f.write("_atom_site_fract_x\n")
        f.write("_atom_site_fract_y\n")
        f.write("_atom_site_fract_z\n")
        df_ATOMIC_POSITIONS = params_structure["df_ATOMIC_POSITIONS"]
        for i in range(df_ATOMIC_POSITIONS.shape[0]):
            label = df_ATOMIC_POSITIONS["label"].iloc[i]
            symbol = df_ATOMIC_POSITIONS["symbol"].iloc[i]
            x = float(df_ATOMIC_POSITIONS["str_x"].iloc[i])
            y = float(df_ATOMIC_POSITIONS["str_y"].iloc[i])
            z = float(df_ATOMIC_POSITIONS["str_z"].iloc[i])
            f.write(f"{label} {symbol} {x:.5f} {y:.5f} {z:.5f}\n")


# In progress, unusable
def symmop_to_xyz(symmop):  # SymmOpを_space_group_symop_operation_xyzの形式に変換
    # https://qiita.com/ojiya/items/e98a9dd5cb6cd7ad38ca
    matrix = symmop.rotation_matrix
    translation = symmop.translation_vector
    xyz = []
    for i, axis in enumerate(["x", "y", "z"]):
        terms = []
        for j, symbol in enumerate(["x", "y", "z"]):
            coeff = matrix[i, j]
            if coeff == 1:
                terms.append(symbol)
            elif coeff == -1:
                terms.append(f"-{symbol}")
            elif coeff != 0:
                terms.append(f"{coeff:.1f}*{symbol}")

        if translation[i] != 0:
            match translation[i]:
                case 0.5:
                    terms.append(f"1/2")
                case 0.25:
                    terms.append(f"1/4")
                case 0.75:
                    terms.append(f"3/4")
                case _:
                    terms.append(f"{translation[i]:.3f}")
        xyz.append(" + ".join(terms) if terms else "0")
    return f"{xyz[0]}, {xyz[1]}, {xyz[2]}"


def check_output(import_out_path):
    with open(import_out_path, "r") as f:
        out_data = f.readlines()
        list_total_energy = []
        list_Total_force = []
        list_P = []
        for i in range(len(out_data)):
            line_split = out_data[i].split()
            if "!" in out_data[i] and "total energy" in out_data[i]:
                for j in range(len(line_split)):
                    try:
                        list_total_energy.append(float(line_split[j]))
                        break
                    except:
                        None
            if "Total force =" in out_data[i]:
                for j in range(len(line_split)):
                    try:
                        list_Total_force.append(float(line_split[j]))
                        break
                    except:
                        None
            if "(kbar)" in out_data[i] and "P" in out_data[i]:
                for j in range(len(line_split)):
                    if "P" in line_split[j]:
                        list_P.append(float(line_split[j + 1]))
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    plt.plot(range(len(list_P)), list_P)
    axes[0].set_title("total_energy")
    axes[0].plot(range(len(list_total_energy)), list_total_energy, label="total_energy", color="r")
    axes[1].set_title("Total_force")
    axes[1].plot(range(len(list_Total_force)), list_Total_force, label="Total_force", color="g")
    axes[2].set_title("P")
    axes[2].plot(range(len(list_P)), list_P, label="P", color="b")
    for ax in axes:
        ax.legend()
        ax.grid(True)  # グリッドを表示
    plt.tight_layout()
    plt.show()
