import pandas as pd
import numpy as np
import re
from .utils import check_args, round_half

# https://qiita.com/xa_member/items/a519c53b63df75a68894


def output_to_params(calc, import_out_path, base_params):
    check_args(False, calc=calc, save_path=import_out_path, params_structure=base_params)
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
