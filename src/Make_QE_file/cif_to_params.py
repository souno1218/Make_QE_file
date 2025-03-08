import pandas as pd
from .utils import round_half, check_params

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
