from moyopy import SpaceGroupType
import datetime

# from pymatgen.symmetry.groups import SpaceGroup


def write_new_cif(out_path, material_name, params_structure):
    date = str(datetime.datetime.now().strftime("%Y-%m-%d"))
    time = str(datetime.datetime.now().strftime("%H:%M:%S"))
    space_group_number = params_structure["space_group_number"]
    hm_symbol = SpaceGroupType(space_group_number).hm_short  # Hermann-Mauguin 記号を取得
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
