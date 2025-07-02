import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import glob, re, os


def plot_band(
    gnu_path,
    k_point_divisions,
    brilloin_zone_path,
    EFermi=None,
    highest_occupied=None,
    title=None,
    is_save=False,
    is_plot=False,
    savefig_path=None,
    ylim=[-5, 5],
    figsize=(10, 7),
):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)  # Axesオブジェクトを明示的に取得

    if not title is None:
        fig.suptitle(title, fontsize=16)

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
            y_bands = []  # 各バンドのY値を格納するリストのリスト

            # 最初のバンドのX軸データを抽出
            # separate_indices[0] + 1 は最初のバンドデータの開始行
            # separate_indices[1] は最初のバンドデータの終了行（次の区切り）
            if separate_indices[0] + 1 >= separate_indices[1]:
                raise ValueError(f"Error: First band data block is empty in GNU file: {gnu_path}")

            for line_idx in range(separate_indices[0] + 1, separate_indices[1]):
                parts = data_lines[line_idx].split()
                if len(parts) >= 2:  # XとYのデータがあることを確認
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
            for t in range(len(separate_indices) - 1):  # バンドの数だけループ
                band_y_vals = []
                start_line = separate_indices[t] + 1
                end_line = separate_indices[t + 1]

                for line_idx in range(start_line, end_line):
                    parts = data_lines[line_idx].split()
                    if len(parts) >= 2:  # XとYのデータがあることを確認
                        try:
                            band_y_vals.append(float(parts[1]))
                        except ValueError:
                            # 数値でない場合はスキップ。データの整合性が重要ならraiseでも良い
                            print(
                                f"Warning: Skipping malformed Y-data line in GNU file: '{data_lines[line_idx].strip()}'"
                            )
                            continue

                # バンドデータがX軸データと長さが一致するかチェック
                if len(band_y_vals) == len(x_data_np):
                    y_bands.append(np.array(band_y_vals))
                else:
                    print(
                        f"Warning: Band {t+1} data length ({len(band_y_vals)}) does not match X-axis data length ({len(x_data_np)}). Skipping this band."
                    )

            if not y_bands:
                raise ValueError(f"Error: No valid Y-axis band data parsed from GNU file: {gnu_path}")

            y_data_np = np.array(y_bands)  # YデータをNumpy配列に変換

    except (FileNotFoundError, ValueError) as e:
        print(f"Skipping band plot: {e}")
        # プロット処理を続行できないので、ここで関数を終了
        plt.close(fig)  # エラーで終了する前にFigureを閉じる
        return

    # バンド構造をプロット
    ax.plot(x_data_np, y_data_np.T - energy_border, color="blue", linewidth=0.8)

    # ブリリアンゾーンパスの垂直線とラベル
    current_x_idx = 0
    for i in range(len(brilloin_zone_path)):
        # インデックスが範囲外にならないようにガード
        if current_x_idx >= len(x_data_np):
            print(f"Warning: Brillouin Zone Path length exceeds data points. Stopping line/label placement.")
            break

        # 垂直線
        ax.axvline(x_data_np[current_x_idx], color="black", linestyle="--", linewidth=0.8)

        # ラベルの整形
        label_text = brilloin_zone_path[i]
        if label_text == "gG":
            label_text = r"$\Gamma$"
        elif label_text == "gS":
            label_text = r"$\Sigma$"

        # Y軸の範囲を考慮してテキストの位置を調整 (下から5%の位置)
        text_y_position = ylim[0] + (ylim[1] - ylim[0]) * 0.05

        ax.text(
            x_data_np[current_x_idx],
            text_y_position,
            label_text,
            va="bottom",
            ha="center",
            fontsize="large",
            bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.7),
        )

        current_x_idx += k_point_divisions[i]

    # 最後のk_point_divisionの後に最後の垂直線を引く（もし必要なら）
    # データポイントの最大値の位置に線が引かれているか確認し、引かれていなければ引く
    if len(x_data_np) > 0 and (
        current_x_idx - k_point_divisions[-1] < len(x_data_np) - 1 or len(brilloin_zone_path) == 0
    ):
        # 最後の点に線が引かれていない、またはパスが指定されていないがデータがある場合
        if len(x_data_np) > 1 and x_data_np[-1] not in [
            ax.lines[j].get_xdata()[0] for j in range(len(ax.lines)) if ax.lines[j].get_xdata() is not None
        ]:
            ax.axvline(x_data_np[-1], color="black", linestyle="--", linewidth=0.8)

    # X軸の範囲を設定
    ax.set_xlim(np.min(x_data_np), np.max(x_data_np))
    # Y軸の範囲を設定
    ax.set_ylim(ylim)

    # X軸の目盛りラベルを非表示に（ブリリアンゾーンパスのラベルを使うため）
    ax.tick_params(labelbottom=False, bottom=False)

    # グリッドを追加
    ax.grid(True, linestyle=":", alpha=0.7)

    # 画像の保存ロジック
    if is_save:
        output_dir = os.path.dirname(savefig_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        fig.savefig(savefig_path, dpi=300, bbox_inches="tight")  # figオブジェクトから保存
        print(f"Plot saved to: {savefig_path}")
    if is_plot:
        plt.show()
    plt.close(fig)


def plot_pdos(
    pdos_dir_path,
    highest_occupied=None,
    EFermi=None,
    title=None,
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
    ax = fig.add_subplot(111)  # Axesオブジェクトを明示的に取得

    if not title is None:
        fig.suptitle(title, fontsize=16)

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
    y_max_plot = -np.inf  # プロットされるデータの最大値を追跡

    # ヘルパー関数: ファイルの読み込みとエラーチェックを共通化
    def _read_data_from_file(file_pattern, single_file=True):
        files = glob.glob(file_pattern)
        if len(files) == 0:
            raise FileNotFoundError(f"Error: No files found matching pattern '{file_pattern}'")
        if single_file and len(files) != 1:
            raise ValueError(f"Error: Expected exactly one file for pattern '{file_pattern}', but found {len(files)}.")
        return files[0] if single_file else files  # 単一ファイルならそのパス、複数ならパスのリストを返す

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
                for line_num in range(1, len(data_lines)):  # 1行目から解析を開始（ヘッダーをスキップ）
                    parts = data_lines[line_num].split()
                    if len(parts) >= 3:  # x, y, integral_yの3列以上を期待
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

            if np.any(x_filtered_indices):  # フィルタリング結果が空でないことを確認
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

            if not pdos_files:  # _wfc_ファイルが見つからない場合を処理
                print(f"Skipping PDOS plot: No _wfc_ files found in {pdos_dir_path}.")
            else:
                # 最初のファイルからxデータを読み込む
                with open(pdos_files[0]) as f:
                    first_file_data = f.readlines()
                    if len(first_file_data) < 2:
                        raise ValueError(f"Error: PDOS file '{pdos_files[0]}' has insufficient data for x-axis.")
                    x_data_pdos = np.array(
                        [float(first_file_data[i].split()[0]) for i in range(1, len(first_file_data))]
                    )

                x_filtered_indices = (xlim[0] < x_data_pdos - energy_border) & (x_data_pdos - energy_border < xlim[1])

                if not np.any(x_filtered_indices):
                    print(f"Warning: No PDOS data found within xlim range {xlim}.")
                else:
                    element_y_data = {}  # 元素ごとのYデータを集計
                    for file in pdos_files:
                        element_match = re.search(r"\((\D*)\)", os.path.basename(file))
                        if not element_match:
                            print(f"Warning: Could not extract element from filename: {file}. Skipping.")
                            continue
                        element = element_match.group(1)

                        with open(file, "r") as f:
                            data = f.readlines()
                            if len(data) < 2:  # 十分な行があるかチェック
                                print(f"Warning: PDOS file '{file}' has insufficient data. Skipping.")
                                continue

                            y_vals = np.array([float(data[t].split()[1]) for t in range(1, len(data))])

                            if element in element_y_data:
                                element_y_data[element] += y_vals
                            else:
                                element_y_data[element] = y_vals

                    if isinstance(color_dict, dict):
                        for el, y_data in element_y_data.items():
                            ax.plot(x_data_pdos - energy_border, y_data, label=el, c=color_dict.get(el, "black"))
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
                    if len(parts) >= 3:  # x, tot_dos, tot_pdosの3列以上を期待
                        try:
                            parsed_data.append((float(parts[0]), float(parts[1]), float(parts[2])))
                        except ValueError:
                            print(
                                f"Warning: Skipping malformed line in Total PDOS/DOS file: '{data_lines[line_num].strip()}'"
                            )
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

    ax.set_ylabel("Density of States")  # Y軸ラベル
    ax.set_xlim(xlim)

    if not ylim is None:
        ax.set_ylim(ylim)
    elif y_max_plot > -np.inf:  # 何らかのデータがプロットされた場合のみ自動Y軸範囲を設定
        # Y軸の下限を0より少し下にする（y_max_plotが正の場合）
        ax.set_ylim((-0.05 * y_max_plot if y_max_plot > 0 else -1, y_max_plot * 1.1))
    else:  # データがプロットされず、ylimも指定されていない場合のデフォルト
        ax.set_ylim((-1, 10))  # 妥当なデフォルト範囲

    # フェルミ準位/最高被占準位を示す垂直線を追加
    ax.axvline(x=0, color="gray", linestyle="--", linewidth=0.8, label="Fermi/Highest Occupied Level")
    ax.legend()
    ax.grid(True, linestyle=":", alpha=0.7)  # グリッド線を追加

    # 画像の保存ロジック
    if is_save:
        output_dir = os.path.dirname(savefig_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        fig.savefig(savefig_path, dpi=300, bbox_inches="tight")  # Figureオブジェクトから保存
        print(f"Plot saved to: {savefig_path}")
    if is_plot:
        plt.show()
    plt.close(fig)


def plot_relax_out(import_out_path, title=None, figsize=(12, 4), save_path=None):  # save_path引数を追加
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
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    if not title is None:
        fig.suptitle(title, fontsize=16)
    data_sets = [
        (list_total_energy, "Total Energy (Ry)", "r"),
        (list_Total_force, "Total Force", "g"),
        (list_P, "Pressure (kbar)", "b"),
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
            ax.text(
                0.5, 0.5, "No data found", transform=ax.transAxes, ha="center", va="center", fontsize=12, color="gray"
            )
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
            if first_val >= data_mean:  # 平均以上なら下にオフセット
                text_first_y_offset = -current_y_offset
                first_va = "top"  # テキストの上端をデータポイントに合わせる
            else:  # 平均より下なら上にオフセット
                text_first_y_offset = current_y_offset
                first_va = "bottom"  # テキストの下端をデータポイントに合わせる

            # X軸方向のオフセットは右に
            ax.text(
                first_idx + current_x_offset,
                first_val + text_first_y_offset,
                f"{first_val:.6f}",
                ha="left",
                va=first_va,
                fontsize=8,
                color="blue",
                bbox=dict(boxstyle="round,pad=0.2", fc="yellow", ec="none", alpha=0.7),
            )
            ax.scatter(first_idx, first_val, color=color, s=20, zorder=5)

            # 最後の点
            last_idx, last_val = len(data_list) - 1, data_list[-1]

            # Y値が平均より上か下かでテキストの垂直位置を決定
            if last_val >= data_mean:  # 平均以上なら下にオフセット
                text_last_y_offset = -current_y_offset
                last_va = "top"
            else:  # 平均より下なら上にオフセット
                text_last_y_offset = current_y_offset
                last_va = "bottom"

            # X軸方向のオフセットは左に
            ax.text(
                last_idx - current_x_offset,
                last_val + text_last_y_offset,
                f"{last_val:.6f}",
                ha="right",
                va=last_va,
                fontsize=8,
                color="blue",
                bbox=dict(boxstyle="round,pad=0.2", fc="yellow", ec="none", alpha=0.7),
            )
            ax.scatter(last_idx, last_val, color=color, s=20, zorder=5)

            # データの差をグラフ内に表示
            data_range = last_val - first_val

            # Y軸のデータ範囲に基づいてテキスト位置を決定
            y_min, y_max = ax.get_ylim()
            y_span = y_max - y_min

            # 最後のデータポイントがサブプロットのY軸のどこにあるか確認
            # 例えば、下1/3にいるならテキストを上に、上1/3にいるならテキストを下に配置
            relative_last_y = (last_val - y_min) / y_span

            if relative_last_y < 0.33:  # グラフの下1/3に最後のデータがある場合
                # テキストを上に配置 (va='top')
                text_va = "top"
                text_y_pos = 0.98  # Axes座標のY軸の上端近く
            elif relative_last_y > 0.66:  # グラフの上1/3に最後のデータがある場合
                # テキストを下に配置 (va='bottom')
                text_va = "bottom"
                text_y_pos = 0.02  # Axes座標のY軸の下端近く
            else:  # 中間の場合 (デフォルトの右下)
                text_va = "bottom"
                text_y_pos = 0.02

            ax.text(
                0.98,
                text_y_pos,  # Xは右端に固定、Yは動的に決定
                f"Diff: {data_range:.6f}",
                transform=ax.transAxes,
                ha="right",
                va=text_va,
                fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", fc="wheat", ec="k", lw=0.5, alpha=0.7),
            )

        ax.grid(True)
        ax.set_xlabel("Step")  # 共通のX軸ラベル

    plt.tight_layout(rect=[0, 0, 1, 0.96])  # 全体タイトルと重ならないように調整

    # 画像の保存ロジック
    if save_path:  # save_pathがNoneでない場合（つまり、保存したい場合）
        if not isinstance(save_path, str) or not save_path:  # 文字列でない、または空文字列の場合
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

    plt.show()  # 保存する場合もしない場合も、画面には表示する


def plot_scf_out(
    kpoints, abc, total_energies, Total_forces, Ps, title=None, figsize=(12, 4), save_path=None
):  # save_path引数を追加

    kpoints = np.array(kpoints)
    index = np.argsort(np.sum((kpoints * abc) ** 2, axis=1))
    kpoints = kpoints[index]
    txt_kpoints = [f"{i[0]} {i[1]} {i[2]}" for i in kpoints]
    total_energies = np.array(total_energies)[index]
    Total_forces = np.array(Total_forces)[index]
    Ps = np.array(Ps)[index]

    fig, axes = plt.subplots(1, 3, figsize=figsize)
    if not title is None:
        fig.suptitle(title, fontsize=16)
    data_sets = [
        (total_energies, "Total Energy (Ry)", "r"),
        (Total_forces, "Total Force", "g"),
        (Ps, "Pressure (kbar)", "b"),
    ]

    for idx, (data_list, title, color) in enumerate(data_sets):
        ax = axes[idx]
        ax.set_title(title)

        # データがない場合はプロットせずにメッセージを表示
        if data_list is None:
            ax.text(
                0.5, 0.5, "No data found", transform=ax.transAxes, ha="center", va="center", fontsize=12, color="gray"
            )
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        ax.plot(range(len(data_list)), data_list, color=color)

        x_min, x_max = ax.get_xlim()
        x_min, x_max = x_min - 0.05 * (x_max - x_min), x_max + 0.05 * (x_max - x_min)
        ax.set_xlim(x_min, x_max)

        y_min, y_max = ax.get_ylim()
        y_max = y_max + 0.05 * (y_max - y_min)
        ax.set_ylim(y_min, y_max)

        # X軸の目盛りを整数にする設定
        ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

        for i in range(len(data_list)):
            ax.text(
                i,
                data_list[i],
                txt_kpoints[i],
                ha="center",
                va="bottom",
                fontsize=8,
                color="blue",
                bbox=dict(boxstyle="round,pad=0.2", fc="yellow", ec="none", alpha=0.7),
            )
            ax.scatter(i, data_list[i], color=color, s=20, zorder=5)

        if len(data_list) > 0:

            # データの差をグラフ内に表示
            data_range = data_list[-1] - data_list[0]

            # Y軸のデータ範囲に基づいてテキスト位置を決定
            y_min, y_max = ax.get_ylim()
            y_span = y_max - y_min

            # 最後のデータポイントがサブプロットのY軸のどこにあるか確認
            # 例えば、下1/3にいるならテキストを上に、上1/3にいるならテキストを下に配置
            relative_last_y = (data_list[-1] - y_min) / y_span

            if relative_last_y < 0.33:  # グラフの下1/3に最後のデータがある場合
                # テキストを上に配置 (va='top')
                text_va = "top"
                text_y_pos = 0.98  # Axes座標のY軸の上端近く
            elif relative_last_y > 0.66:  # グラフの上1/3に最後のデータがある場合
                # テキストを下に配置 (va='bottom')
                text_va = "bottom"
                text_y_pos = 0.02  # Axes座標のY軸の下端近く
            else:  # 中間の場合 (デフォルトの右下)
                text_va = "bottom"
                text_y_pos = 0.02

            ax.text(
                0.98,
                text_y_pos,  # Xは右端に固定、Yは動的に決定
                f"Diff: {data_range:.6f}",
                transform=ax.transAxes,
                ha="right",
                va=text_va,
                fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", fc="wheat", ec="k", lw=0.5, alpha=0.7),
            )

        ax.grid(True)
        ax.set_xlabel("Step")  # 共通のX軸ラベル

    plt.tight_layout(rect=[0, 0, 1, 0.96])  # 全体タイトルと重ならないように調整

    # 画像の保存ロジック
    if save_path:  # save_pathがNoneでない場合（つまり、保存したい場合）
        if not isinstance(save_path, str) or not save_path:  # 文字列でない、または空文字列の場合
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

    plt.show()  # 保存する場合もしない場合も、画面には表示する
