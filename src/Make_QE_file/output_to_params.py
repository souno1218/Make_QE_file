import pandas as pd
import numpy as np
import re
from .utils import check_args, round_half, make_angle

# https://qiita.com/xa_member/items/a519c53b63df75a68894


def output_to_params(calc, import_out_path, base_params):
    check_args(False, calc=calc, save_path=import_out_path, params_structure=base_params)
    if calc == "relax":
        return relax_output_to_params(import_out_path, base_params)
    elif calc == "vc-relax":
        return vc_relax_output_to_params(import_out_path, base_params)
    else:
        raise Exception("not support this calc : {calc}")

# relax
# Cartesian_axes     ： （拡張だけした）座標変換前の初期位置, alat units
# ATOMIC_POSITIONS   : 座標変換後の（複数ある最後が)最終位置, crystal
# crystal axes(0)    : , (v)alat units
# reciprocal axes(0) : , 2 pi/alat

# vc-relax
# Cartesian_axes(0)  : （拡張だけした）座標変換前の初期位置, （初期）alat units
# Cartesian_axes(1)  : （拡張だけした）座標変換前の最終位置, （初期）alat units
# ATOMIC_POSITIONS   : 座標変換後の（複数ある最後が)最終位置, crystal
# crystal axes(0)    : , （初期）alat units
# crystal axes(1)    : , （初期）alat units
# reciprocal axes(0) : , 2 pi/ （c）alat
# reciprocal axes(1) : , 2 pi/ （初期）alat
# CELL_PARAMETERS    : ちょっと正確なcrystal axes（1）, （初期）alat units

# alay_to_crystal = [[1, 0,   0  ], # (= a2c)
#                    [0, b/a, 0  ],
#                    [0, 0,   c/a]]

# -> 計算前, 座標変換前, 拡張後のcrystal
# Cartesian_axes(0) @ a2c
# -> 計算後, 座標変換前, 拡張後のcrystal
# ATOMIC_POSITIONS @ crystal axes(0) @ a2c
# -> 計算前後で各軸ののび
# diag(CELL_PARAMETERS.T @ reciprocal axes(0))

def relax_output_to_params(import_out_path, base_params):
    int_nan = -100
    return_params = base_params.copy()
    with open(import_out_path, "r") as f:
        out_data = f.readlines()
        ibrav = None
        crystal_axes = None
        Cartesian_axes = None
        for i in range(len(out_data)):
            if "bravais-lattice index" in out_data[i]:
                if ibrav is None:
                    ibrav = int(out_data[i].split()[3])
            if "crystal axes: " in out_data[i]:
                if crystal_axes is None:
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
                            float(one_line_split[6]),
                            float(one_line_split[7]),
                            float(one_line_split[8]),
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
                                float(one_line_split[1]),
                                float(one_line_split[2]),
                                float(one_line_split[3]),
                            ]
                        break
    alat_to_crystal = np.array([[1, 0, 0],
                                [0, base_params["a"]/ base_params["b"], 0],
                                [0, 0, base_params["a"]/ base_params["c"]]])

    # Cartesian_axes(0) @ a2c
    Cartesian_axes[["x", "y", "z"]] @= alat_to_crystal
    # ATOMIC_POSITIONS @ crystal axes(0) @ a2c
    ATOMIC_POSITIONS[["x", "y", "z"]] @= (crystal_axes @ alat_to_crystal)
    df_ATOMIC_POSITIONS = return_params["df_ATOMIC_POSITIONS"].copy()
    for i in range(df_ATOMIC_POSITIONS.shape[0]):
        for j in range(Cartesian_axes.shape[0]):
            target_label = df_ATOMIC_POSITIONS["label"].iloc[i]
            target_symbol = df_ATOMIC_POSITIONS["symbol"].iloc[i]
            target_xyz = df_ATOMIC_POSITIONS[["str_x", "str_y", "str_z"]].iloc[i].values.astype("float")
            if target_symbol == Cartesian_axes["atom"].iloc[j]:
                before_xyz = Cartesian_axes[["x", "y", "z"]].iloc[j].values
                after_xyz = ATOMIC_POSITIONS[["x", "y", "z"]].iloc[j].values
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
        ibrav = None
        # celldm = {1:[], 2:[], 3:[], 4:[], 5:[], 6:[]}
        for i in range(len(out_data)):
            #if "celldm(1)" in out_data[i]:
            #    celldm[1].append(float(out_data[i].split()[1]))
            #    celldm[2].append(float(out_data[i].split()[3]))
            #    celldm[3].append(float(out_data[i].split()[5]))
            #if "celldm(4)" in out_data[i]:
            #    celldm[4].append(float(out_data[i].split()[1]))
            #    celldm[5].append(float(out_data[i].split()[3]))
            #    celldm[6].append(float(out_data[i].split()[5]))
            if "bravais-lattice index" in out_data[i]:
                if ibrav is None:
                    ibrav = int(out_data[i].split()[3])
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
                            float(one_line_split[6]),
                            float(one_line_split[7]),
                            float(one_line_split[8]),
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
                                float(one_line_split[1]),
                                float(one_line_split[2]),
                                float(one_line_split[3]),
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
                        break
    alat_to_crystal = np.array([[1, 0, 0],
                                [0, base_params["a"]/ base_params["b"], 0],
                                [0, 0, base_params["a"]/ base_params["c"]]])
    # CELL_PARAMETERS.T @ reciprocal axes(0)
    ratio_abc = np.linalg.norm(CELL_PARAMETERS.T @ reciprocal_axes, axis=0, ord=2)
    return_params["a"] = round_half(base_params["a"] * ratio_abc[0])
    return_params["b"] = round_half(base_params["b"] * ratio_abc[1])
    return_params["c"] = round_half(base_params["c"] * ratio_abc[2])
    match ibrav:
        case 5 | -5 | 12 | -12 | 12 | 13 | -13 | 14:
            return_params["alpha"] = make_angle(CELL_PARAMETERS[1] @ alat_to_crystal, CELL_PARAMETERS[2] @ alat_to_crystal)
            return_params["beta"] = make_angle(CELL_PARAMETERS[2] @ alat_to_crystal, CELL_PARAMETERS[0] @ alat_to_crystal)
            return_params["gamma"] = make_angle(CELL_PARAMETERS[0] @ alat_to_crystal, CELL_PARAMETERS[1] @ alat_to_crystal)
        case _:
            return_params["alpha"] = base_params["alpha"]
            return_params["beta"] = base_params["beta"]
            return_params["gamma"] = base_params["gamma"]
    for k in range(2):
        match k:
            case 0:
                # Cartesian_axes(0) @ a2c
                before_arr = (Cartesian_axes[["x", "y", "z"]].values @ alat_to_crystal).copy()
                # ATOMIC_POSITIONS @ crystal axes(0) @ a2c
                after_arr = (ATOMIC_POSITIONS[["x", "y", "z"]].values @ crystal_axes @ alat_to_crystal).copy()
            case 1:
                before_arr = (Cartesian_axes[["x", "y", "z"]].values @ reciprocal_axes.T).copy()
                after_arr = ATOMIC_POSITIONS[["x", "y", "z"]].values.copy()
        found_all = True
        df_ATOMIC_POSITIONS = return_params["df_ATOMIC_POSITIONS"].copy()
        for i in range(df_ATOMIC_POSITIONS.shape[0]):
            for j in range(Cartesian_axes.shape[0]):
                target_label = df_ATOMIC_POSITIONS["label"].iloc[i]
                target_symbol = df_ATOMIC_POSITIONS["symbol"].iloc[i]
                target_xyz = df_ATOMIC_POSITIONS[["str_x", "str_y", "str_z"]].iloc[i].values.astype("float")
                if target_symbol == Cartesian_axes["atom"].iloc[j]:
                    before_xyz = before_arr[j]
                    after_xyz = after_arr[j]
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
                        found_all = False
                        break
            if not found_all:
                break
        if found_all:
            break
    if not found_all:
        raise
    return_params["df_ATOMIC_POSITIONS"] = df_ATOMIC_POSITIONS
    return return_params


def output_to_nbnd(import_out_path):
    check_args(False, save_path=import_out_path)
    with open(import_out_path, "r") as f:
        out_data = f.readlines()
        for i in range(len(out_data)):
            if "number of Kohn-Sham states=" in out_data[i]:
                return int(out_data[i].split()[4])


def get_EFermi(output_path):
    with open(output_path, "r") as output:
        data = output.readlines()
        for i in range(len(data)):
            #  E (eV)   dos(E)     Int dos(E) EFermi =   10.292 eV
            # the Fermi energy is    11.7660 ev
            if "the Fermi energy is" in data[i]:
                _split = data[i].replace("the Fermi energy is", "").split()
                for j in range(len(_split)):
                    try:
                        EFermi = float(_split[j])
                        break
                    except:
                        None
                return EFermi
            elif "EFermi =" in data[i]:
                _split = data[i].split()
                EFermi_index = -1
                for j in range(len(_split)):
                    if "EFermi" in _split[j]:
                        EFermi_index = j
                        break
                for j in range(EFermi_index + 1, len(_split)):
                    try:
                        EFermi = float(_split[j])
                        break
                    except:
                        None
                return EFermi
    raise ValueError(f"Error: Fermi energy not found in {output_path}.")


def get_highest_occupied(output_path):
    with open(output_path, "r") as output:
        data = output.readlines()
        for i in range(len(data)):
            # highest occupied, lowest unoccupied level (ev):    10.0904   10.6895
            # highest occupied level (ev):    10.0899
            if "highest occupied, lowest unoccupied level" in data[i]:
                _split = data[i].replace("highest occupied, lowest unoccupied level", "").split()
                for j in range(len(_split)):
                    try:
                        highest_occupied_level = float(_split[j])
                        break
                    except:
                        None
                return highest_occupied_level
            elif "highest occupied level" in data[i]:
                _split = data[i].replace("highest occupied level", "").split()
                for j in range(len(_split)):
                    try:
                        highest_occupied_level = float(_split[j])
                        break
                    except:
                        None
                return highest_occupied_level
    raise ValueError(f"Error: highest occupied level not found in {output_path}.")
