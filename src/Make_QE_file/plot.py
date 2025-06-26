import numpy as np
import matplotlib.pyplot as plt
import glob, re, os

def plot_band(
    gnu_path,
    k_point_divisions,
    brilloin_zone_path,
    EFermi=None,
    highest_occupied=None,
    is_save=False,
    is_plot=False,
    savefig_path=None,
    ylim=[-5, 5],
    figsize=(10, 7),
):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111) # Axesオブジェクトを明示的に取得

    # グラフ保存のチェックを最初に行う
    if is_save:
        if savefig_path is None or not isinstance(savefig_path, str) or not savefig_path.strip():
            raise ValueError("Error: If 'is_save' is True, 'savefig_path' must be a valid non-empty string.")

    # gnuファイルが見つかるかチェック
    if not os.path.exists(gnu_path):
        raise FileNotFoundError(f"Error: GNU file not found: {gnu_path}")

    # フェルミ準位または最高被占準位の基準を設定
    if not highest_occupied is None:
        energy_border = highest_occupied
        ax.set_ylabel("E - highest_occupied_level (eV)", fontsize="x-large")
    elif not EFermi is None:
        energy_border = EFermi
        ax.set_ylabel("E - EFermi (eV)", fontsize="x-large")
    else:
        raise ValueError("Error: Either 'highest_occupied' or 'EFermi' must be set.")

    try:
        with open(gnu_path, "r") as bands_gnu:
            data_lines = bands_gnu.readlines()
            
            # 空行を区切りとしてバンドデータを抽出
            # `strip()`を使って、空白のみの行も適切に区切りとして処理
            separate_indices = [-1] + [i for i, line in enumerate(data_lines) if not line.strip()]
            
            # データが全くない場合や、空行だけの場合を考慮
            if len(separate_indices) < 2:
                raise ValueError(f"Error: Invalid or insufficient band data in GNU file: {gnu_path}")

            x_coords = []
            y_bands = [] # 各バンドのY値を格納するリストのリスト

            # 最初のバンドのX軸データを抽出
            # separate_indices[0] + 1 は最初のバンドデータの開始行
            # separate_indices[1] は最初のバンドデータの終了行（次の区切り）
            if separate_indices[0] + 1 >= separate_indices[1]:
                raise ValueError(f"Error: First band data block is empty in GNU file: {gnu_path}")
            
            for line_idx in range(separate_indices[0] + 1, separate_indices[1]):
                parts = data_lines[line_idx].split()
                if len(parts) >= 2: # XとYのデータがあることを確認
                    try:
                        x_coords.append(float(parts[0]))
                    except ValueError:
                        # 数値でない場合はスキップ
                        print(f"Warning: Skipping malformed line in GNU file: '{data_lines[line_idx].strip()}'")
                        continue
                
            x_data_np = np.array(x_coords)
            
            # X軸データが空の場合はエラー
            if len(x_data_np) == 0:
                raise ValueError(f"Error: No X-axis data parsed from GNU file: {gnu_path}")

            # 各バンドのY軸データを抽出
            for t in range(len(separate_indices) - 1): # バンドの数だけループ
                band_y_vals = []
                start_line = separate_indices[t] + 1
                end_line = separate_indices[t+1]
                
                for line_idx in range(start_line, end_line):
                    parts = data_lines[line_idx].split()
                    if len(parts) >= 2: # XとYのデータがあることを確認
                        try:
                            band_y_vals.append(float(parts[1]))
                        except ValueError:
                            # 数値でない場合はスキップ。データの整合性が重要ならraiseでも良い
                            print(f"Warning: Skipping malformed Y-data line in GNU file: '{data_lines[line_idx].strip()}'")
                            continue
                
                # バンドデータがX軸データと長さが一致するかチェック
                if len(band_y_vals) == len(x_data_np):
                    y_bands.append(np.array(band_y_vals))
                else:
                    print(f"Warning: Band {t+1} data length ({len(band_y_vals)}) does not match X-axis data length ({len(x_data_np)}). Skipping this band.")
            
            if not y_bands:
                raise ValueError(f"Error: No valid Y-axis band data parsed from GNU file: {gnu_path}")
            
            y_data_np = np.array(y_bands) # YデータをNumpy配列に変換

    except (FileNotFoundError, ValueError) as e:
        print(f"Skipping band plot: {e}")
        # プロット処理を続行できないので、ここで関数を終了
        plt.close(fig) # エラーで終了する前にFigureを閉じる
        return

    # バンド構造をプロット
    ax.plot(x_data_np, y_data_np.T - energy_border, color='blue', linewidth=0.8)

    # ブリリアンゾーンパスの垂直線とラベル
    current_x_idx = 0
    for i in range(len(brilloin_zone_path)):
        # インデックスが範囲外にならないようにガード
        if current_x_idx >= len(x_data_np):
            print(f"Warning: Brillouin Zone Path length exceeds data points. Stopping line/label placement.")
            break

        # 垂直線
        ax.axvline(x_data_np[current_x_idx], color="black", linestyle='--', linewidth=0.8)

        # ラベルの整形
        label_text = brilloin_zone_path[i]
        if label_text == "gG":
            label_text = r"$\Gamma$"
        elif label_text == "gS":
            label_text = r"$\Sigma$"
        
        # Y軸の範囲を考慮してテキストの位置を調整 (下から5%の位置)
        text_y_position = ylim[0] + (ylim[1] - ylim[0]) * 0.05
        
        ax.text(x_data_np[current_x_idx], text_y_position, label_text,
                va="bottom", ha="center", fontsize="large",
                bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.7))

        current_x_idx += k_point_divisions[i]
        
    # 最後のk_point_divisionの後に最後の垂直線を引く（もし必要なら）
    # データポイントの最大値の位置に線が引かれているか確認し、引かれていなければ引く
    if len(x_data_np) > 0 and (current_x_idx - k_point_divisions[-1] < len(x_data_np) - 1 or len(brilloin_zone_path) == 0):
        # 最後の点に線が引かれていない、またはパスが指定されていないがデータがある場合
        if len(x_data_np) > 1 and x_data_np[-1] not in [ax.lines[j].get_xdata()[0] for j in range(len(ax.lines)) if ax.lines[j].get_xdata() is not None]:
             ax.axvline(x_data_np[-1], color="black", linestyle='--', linewidth=0.8)


    # X軸の範囲を設定
    ax.set_xlim(np.min(x_data_np), np.max(x_data_np))
    # Y軸の範囲を設定
    ax.set_ylim(ylim)
    
    # X軸の目盛りラベルを非表示に（ブリリアンゾーンパスのラベルを使うため）
    ax.tick_params(labelbottom=False, bottom=False)
    
    # グリッドを追加
    ax.grid(True, linestyle=':', alpha=0.7)

    # 画像の保存ロジック
    if is_save:
        output_dir = os.path.dirname(savefig_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        fig.savefig(savefig_path, dpi=300, bbox_inches="tight") # figオブジェクトから保存
        print(f"Plot saved to: {savefig_path}")
    if is_plot:
        plt.show()
    plt.close(fig)


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
    figsize=(8, 6),
):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111) # Axesオブジェクトを明示的に取得

    # グラフ保存のチェックを最初に行う
    if is_save:
        if savefig_path is None or not isinstance(savefig_path, str) or not savefig_path.strip():
            raise ValueError("Error: If 'is_save' is True, 'savefig_path' must be a valid non-empty string.")

    # フェルミ準位または最高被占準位の基準を設定
    if not highest_occupied is None:
        energy_border = highest_occupied
        ax.set_xlabel("E - highest_occupied_level (eV)")
    elif not EFermi is None:
        energy_border = EFermi
        ax.set_xlabel("E - EFermi (eV)")
    else:
        raise ValueError("Error: Either 'highest_occupied' or 'EFermi' must be set.")

    # y軸の最大値を初期化
    y_max_plot = -np.inf # プロットされるデータの最大値を追跡

    # ヘルパー関数: ファイルの読み込みとエラーチェックを共通化
    def _read_data_from_file(file_pattern, single_file=True):
        files = glob.glob(file_pattern)
        if len(files) == 0:
            raise FileNotFoundError(f"Error: No files found matching pattern '{file_pattern}'")
        if single_file and len(files) != 1:
            raise ValueError(f"Error: Expected exactly one file for pattern '{file_pattern}', but found {len(files)}.")
        return files[0] if single_file else files # 単一ファイルならそのパス、複数ならパスのリストを返す

    # dosのプロット
    if "dos" in plot_list:
        try:
            dos_file = _read_data_from_file(f"{pdos_dir_path}/*.dos")
            with open(dos_file, "r") as f:
                data_lines = f.readlines()
                # データ行が十分にあるかチェック
                if len(data_lines) < 2:
                    raise ValueError(f"Error: DOS file '{dos_file}' has insufficient data.")
                
                parsed_data = []
                for line_num in range(1, len(data_lines)): # 1行目から解析を開始（ヘッダーをスキップ）
                    parts = data_lines[line_num].split()
                    if len(parts) >= 3: # x, y, integral_yの3列以上を期待
                        try:
                            parsed_data.append((float(parts[0]), float(parts[1]), float(parts[2])))
                        except ValueError:
                            # 解析できない行はスキップ
                            print(f"Warning: Skipping malformed line in DOS file: '{data_lines[line_num].strip()}'")
                            continue
                
                if not parsed_data:
                    raise ValueError(f"Error: No valid data parsed from DOS file '{dos_file}'.")

                x_data = np.array([d[0] for d in parsed_data])
                y_dos = np.array([d[1] for d in parsed_data])
                y_integral_dos = np.array([d[2] for d in parsed_data])

            # xlim範囲でデータをフィルタリング
            x_filtered_indices = (xlim[0] < x_data - energy_border) & (x_data - energy_border < xlim[1])
            
            if np.any(x_filtered_indices): # フィルタリング結果が空でないことを確認
                ax.plot(x_data - energy_border, y_dos, label="DOS")
                y_max_plot = max(y_max_plot, np.max(y_dos[x_filtered_indices]))
                
                ax.plot(x_data - energy_border, y_integral_dos, label="Integral DOS")
                y_max_plot = max(y_max_plot, np.max(y_integral_dos[x_filtered_indices]))
            else:
                print(f"Warning: No DOS data found within xlim range {xlim}.")

        except (FileNotFoundError, ValueError) as e:
            print(f"Skipping DOS plot: {e}")

    # pdosのプロット
    if "pdos" in plot_list:
        try:
            pdos_files = _read_data_from_file(f"{pdos_dir_path}/*_wfc*", single_file=False)
            
            if not pdos_files: # _wfc_ファイルが見つからない場合を処理
                print(f"Skipping PDOS plot: No _wfc_ files found in {pdos_dir_path}.")
            else:
                # 最初のファイルからxデータを読み込む
                with open(pdos_files[0]) as f:
                    first_file_data = f.readlines()
                    if len(first_file_data) < 2:
                        raise ValueError(f"Error: PDOS file '{pdos_files[0]}' has insufficient data for x-axis.")
                    x_data_pdos = np.array([float(first_file_data[i].split()[0]) for i in range(1, len(first_file_data))])

                x_filtered_indices = (xlim[0] < x_data_pdos - energy_border) & (x_data_pdos - energy_border < xlim[1])
                
                if not np.any(x_filtered_indices):
                    print(f"Warning: No PDOS data found within xlim range {xlim}.")
                else:
                    element_y_data = {} # 元素ごとのYデータを集計
                    for file in pdos_files:
                        element_match = re.search(r"\((\D*)\)", os.path.basename(file))
                        if not element_match:
                            print(f"Warning: Could not extract element from filename: {file}. Skipping.")
                            continue
                        element = element_match.group(1)

                        with open(file, "r") as f:
                            data = f.readlines()
                            if len(data) < 2: # 十分な行があるかチェック
                                print(f"Warning: PDOS file '{file}' has insufficient data. Skipping.")
                                continue
                            
                            y_vals = np.array([float(data[t].split()[1]) for t in range(1, len(data))])
                            
                            if element in element_y_data:
                                element_y_data[element] += y_vals
                            else:
                                element_y_data[element] = y_vals
                    
                    if isinstance(color_dict, dict):
                        for el, y_data in element_y_data.items():
                            ax.plot(x_data_pdos - energy_border, y_data, label=el, c=color_dict.get(el, 'black'))
                            y_max_plot = max(y_max_plot, np.max(y_data[x_filtered_indices]))
                    else:
                        for el, y_data in element_y_data.items():
                            ax.plot(x_data_pdos - energy_border, y_data, label=el)
                            y_max_plot = max(y_max_plot, np.max(y_data[x_filtered_indices]))
        except (FileNotFoundError, ValueError) as e:
            print(f"Skipping PDOS plot: {e}")

    # tot_pdos および tot_dos のプロット（ファイル読み込みを共通化）
    if "tot_pdos" in plot_list or "tot_dos" in plot_list:
        try:
            tot_pdos_file = _read_data_from_file(f"{pdos_dir_path}/*.pdos_tot")
            with open(tot_pdos_file, "r") as f:
                data_lines = f.readlines()
                if len(data_lines) < 2:
                    raise ValueError(f"Error: Total PDOS/DOS file '{tot_pdos_file}' has insufficient data.")
                
                parsed_data = []
                for line_num in range(1, len(data_lines)):
                    parts = data_lines[line_num].split()
                    if len(parts) >= 3: # x, tot_dos, tot_pdosの3列以上を期待
                        try:
                            parsed_data.append((float(parts[0]), float(parts[1]), float(parts[2])))
                        except ValueError:
                            print(f"Warning: Skipping malformed line in Total PDOS/DOS file: '{data_lines[line_num].strip()}'")
                            continue
                
                if not parsed_data:
                    raise ValueError(f"Error: No valid data parsed from Total PDOS/DOS file '{tot_pdos_file}'.")

                x_data_tot = np.array([d[0] for d in parsed_data])
                y_total_dos = np.array([d[1] for d in parsed_data])
                y_total_pdos = np.array([d[2] for d in parsed_data])

                x_filtered_indices = (xlim[0] < x_data_tot - energy_border) & (x_data_tot - energy_border < xlim[1])
                
                if np.any(x_filtered_indices):
                    if "tot_pdos" in plot_list:
                        ax.plot(x_data_tot - energy_border, y_total_pdos, label="Total PDOS")
                        y_max_plot = max(y_max_plot, np.max(y_total_pdos[x_filtered_indices]))
                    if "tot_dos" in plot_list:
                        ax.plot(x_data_tot - energy_border, y_total_dos, label="Total DOS")
                        y_max_plot = max(y_max_plot, np.max(y_total_dos[x_filtered_indices]))
                else:
                    print(f"Warning: No Total PDOS/DOS data found within xlim range {xlim}.")

        except (FileNotFoundError, ValueError) as e:
            print(f"Skipping Total PDOS/DOS plot: {e}")

    ax.set_ylabel("Density of States") # Y軸ラベル
    ax.set_xlim(xlim)

    if not ylim is None:
        ax.set_ylim(ylim)
    elif y_max_plot > -np.inf: # 何らかのデータがプロットされた場合のみ自動Y軸範囲を設定
        # Y軸の下限を0より少し下にする（y_max_plotが正の場合）
        ax.set_ylim((-0.05 * y_max_plot if y_max_plot > 0 else -1, y_max_plot * 1.1))
    else: # データがプロットされず、ylimも指定されていない場合のデフォルト
        ax.set_ylim((-1, 10)) # 妥当なデフォルト範囲

    # フェルミ準位/最高被占準位を示す垂直線を追加
    ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.8, label='Fermi/Highest Occupied Level')
    ax.legend()
    ax.grid(True, linestyle=':', alpha=0.7) # グリッド線を追加

    # 画像の保存ロジック
    if is_save:
        output_dir = os.path.dirname(savefig_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        fig.savefig(savefig_path, dpi=300, bbox_inches="tight") # Figureオブジェクトから保存
        print(f"Plot saved to: {savefig_path}")
    if is_plot:
        plt.show()
    plt.close(fig)
