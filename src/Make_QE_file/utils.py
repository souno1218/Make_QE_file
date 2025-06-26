import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
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


def check_output(import_out_path, prefix=None, figsize=(12, 4), save_path=None): # save_path引数を追加
    if not os.path.exists(import_out_path):
        print(f"Error: Input file not found: {import_out_path}")
        return None
    
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
                    except ValueError: # floatへの変換エラーを具体的に捕獲
                        pass
            if "Total force =" in out_data[i]:
                for j in range(len(line_split)):
                    try:
                        list_Total_force.append(float(line_split[j]))
                        break
                    except ValueError: # floatへの変換エラーを具体的に捕獲
                        pass
            if "(kbar)" in out_data[i] and "P" in out_data[i]:
                for j in range(len(line_split)):
                    if "P" in line_split[j]:
                        try:
                            list_P.append(float(line_split[j + 1]))
                            break # 値が見つかったらループを抜ける
                        except (ValueError, IndexError): # 複数エラーを捕獲
                            pass
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    if not prefix is None:
        fig.suptitle(prefix, fontsize=16)
    data_sets = [
        (list_total_energy, "Total Energy (Ry)", "r"),
        (list_Total_force, "Total Force", "g"),
        (list_P, "Pressure (kbar)", "b")
    ]


    # テキスト表示のオフセット量（データ座標単位）
    # Y軸方向のオフセットをグラフ高さの2%とする
    text_offset_y_factor = 0.02
    # X軸方向のオフセットをグラフ幅の1%とする
    text_offset_x_factor = 0.01

    for idx, (data_list, title, color) in enumerate(data_sets):
        ax = axes[idx]
        ax.set_title(title)

        # データがない場合はプロットせずにメッセージを表示
        if not data_list:
            ax.text(0.5, 0.5, "No data found", transform=ax.transAxes,
                    ha='center', va='center', fontsize=12, color='gray')
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        ax.plot(range(len(data_list)), data_list, color=color)

        # X軸の目盛りを整数にする設定
        ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

        # 最初の点と最後の点に数値をプロット
        if len(data_list) > 0:

            # Y軸の表示範囲を取得してYオフセット量を動的に計算
            y_min, y_max = ax.get_ylim()
            current_y_offset = (y_max - y_min) * text_offset_y_factor
            
            # X軸の表示範囲を取得してXオフセット量を動的に計算
            x_min, x_max = ax.get_xlim()
            current_x_offset = (x_max - x_min) * text_offset_x_factor

            # データリストの平均値を計算 (NaNや空リストを考慮)
            data_mean = np.mean(data_list) if data_list else 0
            
            # 最初の点
            first_idx, first_val = 0, data_list[0]
            
            # Y値が平均より上か下かでテキストの垂直位置を決定
            if first_val >= data_mean: # 平均以上なら下にオフセット
                text_first_y_offset = -current_y_offset
                first_va = 'top' # テキストの上端をデータポイントに合わせる
            else: # 平均より下なら上にオフセット
                text_first_y_offset = current_y_offset
                first_va = 'bottom' # テキストの下端をデータポイントに合わせる

            # X軸方向のオフセットは右に
            ax.text(first_idx + current_x_offset, first_val + text_first_y_offset,
                    f'{first_val:.6f}',
                    ha='left', va=first_va, fontsize=8, color='blue',
                    bbox=dict(boxstyle='round,pad=0.2', fc='yellow', ec='none', alpha=0.7))
            ax.scatter(first_idx, first_val, color=color, s=20, zorder=5)

            # 最後の点
            last_idx, last_val = len(data_list) - 1, data_list[-1]
            
            # Y値が平均より上か下かでテキストの垂直位置を決定
            if last_val >= data_mean: # 平均以上なら下にオフセット
                text_last_y_offset = -current_y_offset
                last_va = 'top'
            else: # 平均より下なら上にオフセット
                text_last_y_offset = current_y_offset
                last_va = 'bottom'

            # X軸方向のオフセットは左に
            ax.text(last_idx - current_x_offset, last_val + text_last_y_offset,
                    f'{last_val:.6f}',
                    ha='right', va=last_va, fontsize=8, color='blue',
                    bbox=dict(boxstyle='round,pad=0.2', fc='yellow', ec='none', alpha=0.7))
            ax.scatter(last_idx, last_val, color=color, s=20, zorder=5)

            # データの差をグラフ内に表示
            data_range = last_val - first_val

            # Y軸のデータ範囲に基づいてテキスト位置を決定
            y_min, y_max = ax.get_ylim()
            y_span = y_max - y_min

            # 最後のデータポイントがサブプロットのY軸のどこにあるか確認
            # 例えば、下1/3にいるならテキストを上に、上1/3にいるならテキストを下に配置
            relative_last_y = (last_val - y_min) / y_span

            if relative_last_y < 0.33: # グラフの下1/3に最後のデータがある場合
                # テキストを上に配置 (va='top')
                text_va = 'top'
                text_y_pos = 0.98 # Axes座標のY軸の上端近く
            elif relative_last_y > 0.66: # グラフの上1/3に最後のデータがある場合
                # テキストを下に配置 (va='bottom')
                text_va = 'bottom'
                text_y_pos = 0.02 # Axes座標のY軸の下端近く
            else: # 中間の場合 (デフォルトの右下)
                text_va = 'bottom'
                text_y_pos = 0.02

            ax.text(0.98, text_y_pos, # Xは右端に固定、Yは動的に決定
                    f'Diff: {data_range:.6f}',
                    transform=ax.transAxes,
                    ha='right', va=text_va, fontsize=9,
                    bbox=dict(boxstyle='round,pad=0.3', fc='wheat', ec='k', lw=0.5, alpha=0.7))

        ax.grid(True)
        ax.set_xlabel("Step") # 共通のX軸ラベル

    plt.tight_layout(rect=[0, 0, 1, 0.96]) # 全体タイトルと重ならないように調整

    # 画像の保存ロジック
    if save_path: # save_pathがNoneでない場合（つまり、保存したい場合）
        if not isinstance(save_path, str) or not save_path: # 文字列でない、または空文字列の場合
            print("Error: Invalid save_path provided. Please provide a valid string path.")
        else:
            try:
                # ディレクトリが存在しない場合は作成
                output_dir = os.path.dirname(save_path)
                if output_dir and not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                plt.savefig(save_path)
                print(f"Plot saved successfully to: {save_path}")
            except Exception as e:
                print(f"Error saving plot to {save_path}: {e}")

    plt.show() # 保存する場合もしない場合も、画面には表示する
