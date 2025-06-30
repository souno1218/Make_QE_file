import pandas as pd
from .utils import round_half, check_params

# https://qiita.com/xa_member/items/a519c53b63df75a68894


def cif_to_params(import_cif_path):
    params_cif = dict()
    list_loop_index = []
    with open(import_cif_path, "r") as f:
        cif_data = f.readlines()
        for i in range(len(cif_data)):
            one_line = cif_data[i].split("#")[0]
            if "_cell_length_a" in cif_data[i]:
                for j in one_line.split():
                    try:
                        params_cif["a"] = round_half(j)
                    except:
                        None
            if "_cell_length_b" in cif_data[i]:
                for j in one_line.split():
                    try:
                        params_cif["b"] = round_half(j)
                    except:
                        None
            if "_cell_length_c" in cif_data[i]:
                for j in one_line.split():
                    try:
                        params_cif["c"] = round_half(j)
                    except:
                        None
            if "_cell_angle_alpha" in cif_data[i]:
                for j in one_line.split():
                    try:
                        params_cif["alpha"] = round_half(j)
                    except:
                        None
            if "_cell_angle_beta" in cif_data[i]:
                for j in one_line.split():
                    try:
                        params_cif["beta"] = round_half(j)
                    except:
                        None
            if "_cell_angle_gamma" in cif_data[i]:
                for j in one_line.split():
                    try:
                        params_cif["gamma"] = round_half(j)
                    except:
                        None
            if "_space_group_IT_number" in cif_data[i]:
                for j in one_line.split():
                    try:
                        params_cif["space_group_number"] = int(j)
                    except:
                        None
            if "_symmetry_Int_Tables_number" in cif_data[i]:
                for j in one_line.split():
                    try:
                        params_cif["space_group_number"] = int(j)
                    except:
                        None
            if "loop_" in one_line:
                list_loop_index.append(i)
        list_loop_index.append(len(cif_data))
        for i in range(len(list_loop_index) - 1):
            for j in range(list_loop_index[i], list_loop_index[i + 1]):
                if "_atom_site" in cif_data[i].split("#")[0]:
                    loop_index = list_loop_index[i]
        loop_label_2_index = {}
        len_labels = 0
        for i in range(loop_index + 1, len(cif_data)):
            one_line = cif_data[i].split("#")[0]
            if "_atom_site_label" in one_line:
                loop_label_2_index["label"] = len_labels
            elif "_atom_site_type_symbol" in one_line:
                loop_label_2_index["symbol"] = len_labels
            elif "_atom_site_fract_x" in one_line:
                loop_label_2_index["str_x"] = len_labels
            elif "_atom_site_fract_y" in one_line:
                loop_label_2_index["str_y"] = len_labels
            elif "_atom_site_fract_z" in one_line:
                loop_label_2_index["str_z"] = len_labels
            if "_atom_site" in one_line:
                len_labels += 1
        labels = []
        for i in range(loop_index + 1, len(cif_data)):
            one_data = cif_data[i].split("#")[0].split()
            if len(one_data) == len_labels:
                labels.append(one_data[loop_label_2_index["label"]])
        n = 0
        for i in range(len(labels)):
            while labels[i] in [*labels[:i], *labels[i+1:]]:
                labels[i] = f'{labels[i]}_{n}'
                n += 1

        df_ATOMIC_POSITIONS = pd.DataFrame(columns=["label", "symbol", "str_x", "str_y", "str_z"])
        n = 0
        for i in range(loop_index + 1, len(cif_data)):
            one_data = cif_data[i].split("#")[0].split()
            if len(one_data) == len_labels:
                one_ATOMIC_POSITIONS = [labels[n]]
                one_ATOMIC_POSITIONS.append(one_data[loop_label_2_index["symbol"]])
                str_x = "{:.05f}".format(round_half(one_data[loop_label_2_index["str_x"]]))
                one_ATOMIC_POSITIONS.append(str_x)
                str_y = "{:.05f}".format(round_half(one_data[loop_label_2_index["str_y"]]))
                one_ATOMIC_POSITIONS.append(str_y)
                str_z = "{:.05f}".format(round_half(one_data[loop_label_2_index["str_z"]]))
                one_ATOMIC_POSITIONS.append(str_z)
                df_ATOMIC_POSITIONS.loc[labels[n]] = one_ATOMIC_POSITIONS
                n += 1
        params_cif["df_ATOMIC_POSITIONS"] = df_ATOMIC_POSITIONS

    # check
    key_check = ["a", "b", "c", "alpha", "beta", "gamma", "space_group_number", "df_ATOMIC_POSITIONS"]
    for i in key_check:
        if not i in params_cif.keys():
            raise Exception(f"not found param {i} in cif")
    check_params(params_cif)
    return params_cif
