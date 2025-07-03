import pandas as pd
import numpy as np
from decimal import Decimal, ROUND_HALF_UP
import os, glob, re

# https://qiita.com/xa_member/items/a519c53b63df75a68894


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
                if not isinstance(kwargs[i], dict):
                    raise TypeError(f"The argument {i} must be dict")
                check_params(kwargs[i])
            case "min_kpoint":
                if not isinstance(kwargs[i], int):
                    raise TypeError(f"The argument {i} must be str")
                if kwargs[i] < 0:
                    raise ValueError(f"The argument {i} must be str")
            case "K_point_Density_Product":
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
            case "emax":
                if not isinstance(kwargs[i], (int, float)):
                    raise TypeError(f"The argument {i} must be int or float")
            case "emin":
                if not isinstance(kwargs[i], (int, float)):
                    raise TypeError(f"The argument {i} must be int or float")
            case "deltae":
                if not isinstance(kwargs[i], (int, float)):
                    raise TypeError(f"The argument {i} must be int or float")


def round_half(num_or_arr):
    if "ndarray" in str(type(num_or_arr)):
        return_arr = np.empty_like(num_or_arr)
        for i in range(len(num_or_arr.flat)):
            return_arr.flat[i] = float(Decimal(str(float(num_or_arr.flat[i]))).quantize(Decimal("1e-5"), ROUND_HALF_UP))
        return return_arr
    else:
        return float(Decimal(str(num_or_arr)).quantize(Decimal("1e-5"), ROUND_HALF_UP))


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


def check_scf_out(import_out_path):
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
                    except ValueError:  # floatへの変換エラーを具体的に捕獲
                        pass
            if "Total force =" in out_data[i]:
                for j in range(len(line_split)):
                    try:
                        list_Total_force.append(float(line_split[j]))
                        break
                    except ValueError:  # floatへの変換エラーを具体的に捕獲
                        pass
            if "(kbar)" in out_data[i] and "P" in out_data[i]:
                for j in range(len(line_split)):
                    if "P" in line_split[j]:
                        try:
                            list_P.append(float(line_split[j + 1]))
                            break  # 値が見つかったらループを抜ける
                        except (ValueError, IndexError):  # 複数エラーを捕獲
                            pass
    if len(list_total_energy) == 1:
        if len(list_Total_force) == 1:
            if len(list_P) == 1:
                return list_total_energy[0], list_Total_force[0], list_P[0]
    # print(len(list_total_energy), len(list_Total_force), len(list_P))


def make_angle(vec0, vec1):
    norm_vec0 = np.linalg.norm(vec0, ord=2)
    norm_vec1 = np.linalg.norm(vec1, ord=2)
    return np.degrees(np.arccos((vec0 @ vec1) / (norm_vec0 * norm_vec1)))


def bohr2angs(x):
    return x * 0.5291772083
