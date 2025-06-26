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
    plt.figure(figsize=figsize)

    # グラフ保存のチェックを最初に行う
    if is_save:
        if savefig_path is None or not isinstance(savefig_path, str) or not savefig_path.strip():
            raise ValueError("エラー: 'is_save' が True の場合、'savefig_path' に有効なパス文字列を指定してください。")

    # gnuファイルが見つかるかチェック
    if not os.path.exists(gnu_path):
        raise FileNotFoundError(f"エラー: gnuファイルが見つかりません: {gnu_path}")

    # フェルミ準位または最高被占準位の基準を設定
    if not highest_occupied is None:
        border = highest_occupied
        plt.ylabel("E - highest_occupied_level (eV)", fontsize="x-large") # フォントサイズ調整
    elif not EFermi is None:
        border = EFermi
        plt.ylabel("E - EFermi (eV)", fontsize="x-large") # フォントサイズ調整
    else:
        raise ValueError("エラー: 'highest_occupied' または 'EFermi' のいずれかを設定する必要があります。")

    try:
        with open(gnu_path, "r") as bands_gnu:
            data = bands_gnu.readlines()
            
            # 空行を区切りとしてバンドデータを抽出
            separate_indices = [-1] + [i for i, j in enumerate(data) if not j.strip()] # strip()で空白のみの行も処理
            
            # データが全くない場合や、空行だけの場合を考慮
            if len(separate_indices) < 2:
                raise ValueError(f"gnuファイル '{gnu_path}' に有効なバンドデータがありません。")

            x_data = []
            y_data = [] # 各バンドのY値を格納するリストのリスト

            # 最初のバンドのX軸データを抽出
            # separate_indices[0] + 1 は最初のバンドデータの開始行
            # separate_indices[1] は最初のバンドデータの終了行（次の区切り）
            if separate_indices[0] + 1 >= separate_indices[1]:
                raise ValueError(f"gnuファイル '{gnu_path}' の最初のバンドデータが空です。")
            
            for line_idx in range(separate_indices[0] + 1, separate_indices[1]):
                parts = data[line_idx].split()
                if len(parts) >= 2: # XとYのデータがあることを確認
                    try:
                        x_data.append(float(parts[0]))
                    except ValueError:
                        continue # 数値でない場合はスキップ
                
            x = np.array(x_data)
            
            # X軸データが空の場合はエラー
            if len(x) == 0:
                raise ValueError(f"gnuファイル '{gnu_path}' からX軸データを解析できませんでした。")

            # 各バンドのY軸データを抽出
            for t in range(len(separate_indices) - 1): # バンドの数だけループ
                band_y_vals = []
                start_line = separate_indices[t] + 1
                end_line = separate_indices[t+1]
                
                for line_idx in range(start_line, end_line):
                    parts = data[line_idx].split()
                    if len(parts) >= 2: # XとYのデータがあることを確認
                        try:
                            band_y_vals.append(float(parts[1]))
                        except ValueError:
                            # 数値でない場合はスキップ。データの整合性が重要ならraiseでも良い
                            continue
                
                # バンドデータがX軸データと長さが一致するかチェック
                if len(band_y_vals) == len(x):
                    y_data.append(np.array(band_y_vals))
                else:
                    print(f"警告: バンド {t+1} のデータ長がX軸データと一致しません。スキップします。")
            
            if not y_data:
                raise ValueError(f"gnuファイル '{gnu_path}' から有効なY軸バンドデータが解析できませんでした。")
            
            y = np.array(y_data) # YデータをNumpy配列に変換

    except (FileNotFoundError, ValueError) as e:
        print(f"バンドプロットをスキップします: {e}")
        # プロット処理を続行できないので、ここで関数を終了
        plt.close() # エラーで終了する前にFigureを閉じる
        return


    # バンド構造をプロット
    plt.plot(x, y.T - border, color='blue', linewidth=0.8) # 線色と太さを指定

    # ブリリアンゾーンパスの垂直線とラベル
    current_x_index = 0
    for i in range(len(brilloin_zone_path)):
        if current_x_index >= len(x): # インデックスが範囲外にならないようにガード
            print(f"警告: Brillouin Zone Path の長さがデータポイントの数を超えています。")
            break

        # 垂直線
        plt.axvline(x[current_x_index], color="black", linestyle='--', linewidth=0.8) # plt.vlinesよりaxvlineが推奨される

        # ラベルの整形
        text_label = brilloin_zone_path[i]
        if text_label == "gG":
            text_label = r"$\Gamma$"
        elif text_label == "gS":
            text_label = r"$\Sigma$"
        # textのY位置をylimの下限より少し上にする、またはグラフの上端に配置
        # Y軸の範囲を考慮してテキストの位置を調整
        text_y_position = ylim[0] + (ylim[1] - ylim[0]) * 0.05 # 下から5%の位置
        
        plt.text(x[current_x_index], text_y_position, text_label,
                 va="bottom", ha="center", fontsize="large", # va='bottom'で下限より少し上に
                 bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.7)) # 背景追加

        current_x_index += k_point_divisions[i]
        
        # 最後のk_point_divisionの後に最後の垂直線を引く
        if i == len(brilloin_zone_path) - 1 and current_x_index - k_point_divisions[i] < len(x) -1 :
             plt.axvline(x[-1], color="black", linestyle='--', linewidth=0.8) # 最後のX値に線を引く

    # X軸の範囲を設定
    plt.xlim(np.min(x), np.max(x))
    # Y軸の範囲を設定
    plt.ylim(ylim)
    
    # X軸の目盛りラベルを非表示に（ブリリアンゾーンパスのラベルを使うため）
    plt.tick_params(labelbottom=False, bottom=False)
    
    # グリッドを追加
    plt.grid(True, linestyle=':', alpha=0.7)

    # 画像の保存ロジック
    if is_save:
        output_dir = os.path.dirname(savefig_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plt.savefig(savefig_path, dpi=300, bbox_inches="tight") # DPIを300に、bbox_inches="tight"で余白自動調整
        print(f"プロットを保存しました: {savefig_path}")
    if is_plot:
        plt.show()
    plt.close()


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

    plt.figure(figsize=figsize)

    # グラフ保存のチェックを最初に行う
    if is_save:
        # savefig_pathがNoneか、文字列でないか、または空文字列の場合にエラー
        if savefig_path is None or not isinstance(savefig_path, str) or not savefig_path.strip():
            raise ValueError("エラー: 'is_save' が True の場合、'savefig_path' に有効なパス文字列を指定してください。")

    # フェルミ準位または最高被占準位の基準を設定
    if not highest_occupied is None:
        border = highest_occupied
        plt.xlabel("E - highest_occupied_level (eV)")
    elif not EFermi is None:
        border = EFermi
        plt.xlabel("E - EFermi (eV)")
    else:
        # どちらの基準も設定されていない場合にエラー
        raise ValueError("エラー: 'highest_occupied' または 'EFermi' のいずれかを設定する必要があります。")

    y_max = -1

    # ヘルパー関数: ファイルの読み込みとエラーチェックを共通化
    def _read_data_from_file(file_pattern, single_file=True):
        files = glob.glob(file_pattern)
        if len(files) == 0:
            raise FileNotFoundError(f"エラー: '{file_pattern}' に一致するファイルが見つかりません。")
        if single_file and len(files) != 1:
            raise ValueError(f"エラー: '{file_pattern}' のファイルは1つだけを期待しましたが、{len(files)}個見つかりました。")
        return files[0] if single_file else files # 単一ファイルならそのパス、複数ならパスのリストを返す

    # dosのプロット
    if "dos" in plot_list:
        try:
            dos_file = _read_data_from_file(f"{pdos_dir_path}/*.dos")
            with open(dos_file, "r") as dos:
                data = dos.readlines()
                # データ行が十分にあるかチェック
                if len(data) < 2:
                    raise ValueError(f"DOSファイル '{dos_file}' に十分なデータがありません。")
                
                parsed_data = []
                for line_num in range(1, len(data)): # 1行目から解析を開始（ヘッダーをスキップ）
                    parts = data[line_num].split()
                    if len(parts) >= 3: # x, y, integral_yの3列以上を期待
                        try:
                            parsed_data.append((float(parts[0]), float(parts[1]), float(parts[2])))
                        except ValueError:
                            # 解析できない行はスキップ
                            continue
                
                if not parsed_data:
                    raise ValueError(f"DOSファイル '{dos_file}' から有効なデータが解析できませんでした。")

                x = np.array([d[0] for d in parsed_data])
                y = np.array([d[1] for d in parsed_data])
                integral_y = np.array([d[2] for d in parsed_data])

            TF = (xlim[0] < x - border) & (x - border < xlim[1])
            
            # TF（フィルタリング結果）が空でないことを確認してから最大値を取る
            if np.any(TF):
                plt.plot(x - border, y, label="dos")
                y_max = max(y_max, np.max(y[TF])) # y[TF]を直接使用
                plt.plot(x - border, integral_y, label="integral dos")
                y_max = max(y_max, np.max(integral_y[TF]))
            else:
                print(f"警告: xlim範囲 {xlim} 内にDOSデータが見つかりませんでした。")

        except (FileNotFoundError, ValueError) as e:
            print(f"DOSプロットをスキップします: {e}")

    # pdosのプロット
    if "pdos" in plot_list:
        try:
            pdos_files = _read_data_from_file(f"{pdos_dir_path}/*_wfc*", single_file=False)
            
            if not pdos_files: # _wfc_ファイルが見つからない場合を処理
                print(f"PDOSプロットをスキップします: {pdos_dir_path} に _wfc_ ファイルが見つかりませんでした。")
            else:
                # 最初のファイルからxデータを読み込む
                with open(pdos_files[0]) as pdos_f:
                    first_file_data = pdos_f.readlines()
                    if len(first_file_data) < 2:
                        raise ValueError(f"PDOSファイル '{pdos_files[0]}' にx軸用の十分なデータがありません。")
                    x = np.array([float(first_file_data[i].split()[0]) for i in range(1, len(first_file_data))])

                TF = (xlim[0] < x - border) & (x - border < xlim[1])
                
                if not np.any(TF):
                    print(f"警告: xlim範囲 {xlim} 内にPDOSデータが見つかりませんでした。")
                else:
                    dict_y = {}
                    for file in pdos_files:
                        # ファイル名から元素名を抽出
                        element_match = re.search(r"\((\D*)\)", os.path.basename(file))
                        if not element_match:
                            print(f"警告: ファイル名 '{file}' から元素名を抽出できませんでした。スキップします。")
                            continue
                        element = element_match.group(1)

                        with open(file, "r") as pdos_f:
                            data = pdos_f.readlines()
                            if len(data) < 2: # 十分な行があるかチェック
                                print(f"警告: PDOSファイル '{file}' に十分なデータがありません。スキップします。")
                                continue
                            
                            y_vals = np.array([float(data[t].split()[1]) for t in range(1, len(data))])
                            
                            if element in dict_y:
                                dict_y[element] += y_vals
                            else:
                                dict_y[element] = y_vals
                    
                    if isinstance(color_dict, dict):
                        for el, y_data in dict_y.items():
                            # color_dictにキーがない場合も安全に処理 (getメソッドを使用)
                            plt.plot(x - border, y_data, label=el, c=color_dict.get(el, 'black'))
                            y_max = max(y_max, np.max(y_data[TF]))
                    else:
                        for el, y_data in dict_y.items():
                            plt.plot(x - border, y_data, label=el)
                            y_max = max(y_max, np.max(y_data[TF]))
        except (FileNotFoundError, ValueError) as e:
            print(f"PDOSプロットをスキップします: {e}")

    # tot_pdos および tot_dos のプロット（ファイル読み込みを共通化）
    if "tot_pdos" in plot_list or "tot_dos" in plot_list:
        try:
            tot_pdos_file = _read_data_from_file(f"{pdos_dir_path}/*.pdos_tot")
            with open(tot_pdos_file, "r") as dos:
                data = dos.readlines()
                if len(data) < 2:
                    raise ValueError(f"Total PDOS/DOSファイル '{tot_pdos_file}' に十分なデータがありません。")
                
                parsed_data = []
                for line_num in range(1, len(data)):
                    parts = data[line_num].split()
                    if len(parts) >= 3: # x, tot_dos, tot_pdosの3列以上を期待
                        try:
                            parsed_data.append((float(parts[0]), float(parts[1]), float(parts[2])))
                        except ValueError:
                            continue
                
                if not parsed_data:
                    raise ValueError(f"Total PDOS/DOSファイル '{tot_pdos_file}' から有効なデータが解析できませんでした。")

                x = np.array([d[0] for d in parsed_data])
                y_tot_dos = np.array([d[1] for d in parsed_data])
                y_tot_pdos = np.array([d[2] for d in parsed_data])

                TF = (xlim[0] < x - border) & (x - border < xlim[1])
                
                if np.any(TF):
                    if "tot_pdos" in plot_list:
                        plt.plot(x - border, y_tot_pdos, label="tot pdos")
                        y_max = max(y_max, np.max(y_tot_pdos[TF]))
                    if "tot_dos" in plot_list:
                        plt.plot(x - border, y_tot_dos, label="tot dos")
                        y_max = max(y_max, np.max(y_tot_dos[TF]))
                else:
                    print(f"警告: xlim範囲 {xlim} 内にTotal PDOS/DOSデータが見つかりませんでした。")

        except (FileNotFoundError, ValueError) as e:
            print(f"Total PDOS/DOSプロットをスキップします: {e}")

    plt.ylabel("状態密度 (Density of States)") # Y軸ラベルをより一般的なものに変更
    plt.xlim(xlim)

    if not ylim is None:
        plt.ylim(ylim)
    elif y_max > -1: # 何らかのデータがプロットされた場合のみ自動Y軸範囲を設定
        # Y軸の下限を0より少し下にする（y_maxが正の場合）
        plt.ylim((-0.05 * y_max if y_max > 0 else -1, y_max * 1.1))
    else: # データがプロットされず、ylimも指定されていない場合のデフォルト
        plt.ylim((-1, 10)) # 妥当なデフォルト範囲

    # フェルミ準位/最高被占準位を示す垂直線を追加
    plt.axvline(x=0, color='gray', linestyle='--', linewidth=0.8, label='フェルミ準位/最高被占準位')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.7) # グリッド線を追加して可読性向上

    # 画像の保存ロジック
    if is_save:
        # 保存先のディレクトリが存在しない場合は作成
        output_dir = os.path.dirname(savefig_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        # DPIを300に上げて高画質化、bbox_inches="tight"で余白を自動調整
        plt.savefig(savefig_path, dpi=300, bbox_inches="tight")
        print(f"プロットを保存しました: {savefig_path}")
    if is_plot:
        plt.show()
    plt.close()
