import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, glob, re
from decimal import Decimal, ROUND_HALF_UP

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


def read_EFermi(output_path):
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


def read_highest_occupied(output_path):
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


def band_plot(
    gnu_path,
    k_point_divisions,
    brilloin_zone_path,
    EFermi=None,
    highest_occupied=None,
    is_save=False,
    is_plot=False,
    savefig_path=None,
    ylim=[-5, 5],
):
    plt.clf()
    if is_save and (savefig_path is None):
        raise ValueError(f"Error: if save pic, need savefig_path")
    with open(gnu_path, "r") as bands_gnu:
        data = bands_gnu.readlines()
        separate_index = [-1]
        separate_index += [i for i, j in enumerate(data) if len(j.split()) == 0]
        x = np.array([float(data[i].split()[0]) for i in range(separate_index[0] + 1, separate_index[1])])
        y = np.empty((len(separate_index) - 1, x.shape[0]), dtype="float64")
        for t in range(len(separate_index) - 1):
            for i, j in enumerate(range(separate_index[t] + 1, separate_index[t + 1])):
                y[t, i] = float(data[j].split()[1])
        y_ = np.sort(y[y - EFermi > ylim[0]])
        highest_occupied = y_[np.argmax(np.diff(y_))]
        if not highest_occupied is None:
            plt.plot(x, y.T - highest_occupied)
            plt.ylabel("E - highest_occupied_level (eV)", fontsize="xx-large")
        elif not EFermi is None:
            plt.plot(x, y.T - EFermi)
            plt.ylabel("E - EFermi (eV)", fontsize="xx-large")
        else:
            raise ValueError(f"Error: need set highest_occupied or EFermi")
    index = 0
    for i in range(len(brilloin_zone_path)):
        plt.vlines(x[index], ylim[0], ylim[1], color="black")
        text = brilloin_zone_path[i]
        if text == "gG":
            text = r"$\Gamma$"
        elif text == "gS":
            text = r"$\Sigma$"
        plt.text(x[index], ylim[0] - 0.5, text, va="top", ha="center", fontsize="xx-large")
        index += k_point_divisions[i]
    plt.xlim(np.min(x), np.max(x))
    plt.ylim(ylim)
    plt.tick_params(labelbottom=False, bottom=False)
    if is_save:
        plt.savefig(savefig_path, dpi=200, bbox_inches="tight")
    if is_plot:
        plt.show()


def plot_pdos(
    pdos_dir_path,
    highest_occupied=None,
    EFermi=None,
    plot_list=["pdos"],
    xlim=[-10, 10],
    savefig_path=None,
    ylim=None,
    is_save=False,
    is_plot=False,
    color_dict=None,
):
    plt.clf()
    if is_save and (savefig_path is None):
        raise ValueError(f"Error: if save pic, need savefig_path")
    if not highest_occupied is None:
        border = highest_occupied
        plt.xlabel("E - highest_occupied_level (eV)")
    elif not EFermi is None:
        border = EFermi
        plt.xlabel("E - EFermi (eV)")
    else:
        raise ValueError(f"Error: need set highest_occupied or EFermi")
    y_max = -1
    # plot_list => dos, pdos, tot_pdos, tot_pdos, tot_dos
    if "dos" in plot_list:
        dos_path = f"{pdos_dir_path}/*.dos"
        dos_files = glob.glob(dos_path)
        if len(dos_files) == 0:
            raise
        elif len(dos_files) != 1:
            raise
        with open(dos_files[0], "r") as dos:
            data = dos.readlines()
            x = np.array([float(data[i].split()[0]) for i in range(1, len(data) - 1)])
            y = np.array([float(data[i].split()[1]) for i in range(1, len(data) - 1)])
            integral_y = np.array([float(data[i].split()[2]) for i in range(1, len(data) - 1)])
        TF = (xlim[0] < x - border) & (x - border < xlim[1])
        plt.plot(x - border, y.T, label="dos")
        y_max = max(y_max, np.max(y.T[TF]))
        plt.plot(x - border, integral_y.T, label="integral dos")
        y_max = max(y_max, np.max(integral_y.T[TF]))
    if "pdos" in plot_list:
        files = glob.glob(f"{pdos_dir_path}/*_wfc*")
        elements = [
            re.search(r"\((\D*)\)", file.split("/")[-1].split(".")[-1].split("#")[-2]).group(1) for file in files
        ]
        # site_elements=[int(re.search(r"(\d*)\(",file.split("/")[-1].split(".")[-1].split("#")[-2]).group(1)) for file in files]
        # orbit=[file.split("/")[-1].split(".")[-1].split("#")[-1] for file in files]
        with open(files[0]) as pdos:
            data = pdos.readlines()
            x = np.array([float(data[i].split()[0]) for i in range(1, len(data) - 1)])
        TF = (xlim[0] < x - border) & (x - border < xlim[1])
        dict_y = {}
        for i, j in enumerate(files):
            with open(j, "r") as pdos:
                data = pdos.readlines()
                if elements[i] in dict_y.keys():
                    dict_y[elements[i]] += np.array([float(data[t].split()[1]) for t in range(1, len(data) - 1)])
                else:
                    dict_y[elements[i]] = np.array([float(data[t].split()[1]) for t in range(1, len(data) - 1)])
        if isinstance(color_dict, dict):
            for i in dict_y.keys():
                plt.plot(x - border, dict_y[i].T, label=i, c=color_dict[i])
                y_max = max(y_max, np.max(dict_y[i].T[TF]))
        else:
            for i in dict_y.keys():
                plt.plot(x - border, dict_y[i].T, label=i)
                y_max = max(y_max, np.max(dict_y[i].T[TF]))
    if "tot_pdos" in plot_list:
        tot_pdos_path = f"{pdos_dir_path}/*.pdos_tot"
        tot_pdos_files = glob.glob(tot_pdos_path)
        if len(tot_pdos_files) == 0:
            raise
        elif len(tot_pdos_files) != 1:
            raise
        with open(tot_pdos_files[0], "r") as dos:
            data = dos.readlines()
            x = np.array([float(data[i].split()[0]) for i in range(1, len(data) - 1)])
            y = np.array([float(data[i].split()[2]) for i in range(1, len(data) - 1)])
        plt.plot(x - border, y.T, label="tot pdos")
        TF = (xlim[0] < x - border) & (x - border < xlim[1])
        y_max = max(y_max, np.max(y.T[TF]))
    if "tot_dos" in plot_list:
        tot_pdos_path = f"{pdos_dir_path}/*.pdos_tot"
        tot_pdos_files = glob.glob(tot_pdos_path)
        if len(tot_pdos_files) == 0:
            raise
        elif len(tot_pdos_files) != 1:
            raise
        with open(tot_pdos_files[0], "r") as dos:
            data = dos.readlines()
            x = np.array([float(data[i].split()[0]) for i in range(1, len(data) - 1)])
            y = np.array([float(data[i].split()[1]) for i in range(1, len(data) - 1)])
        plt.plot(x - border, y.T, label="tot dos")
        TF = (xlim[0] < x - border) & (x - border < xlim[1])
        y_max = max(y_max, np.max(y.T[TF]))
    plt.ylabel("PDoS")
    plt.xlim(xlim)
    if not ylim is None:
        plt.ylim(ylim)
    else:
        plt.ylim((-1, y_max * 1.1))
    plt.legend()
    if is_save:
        plt.savefig(savefig_path, dpi=200, bbox_inches="tight")
    if is_plot:
        plt.show()
