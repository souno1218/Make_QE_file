# Make_QE_file 使い方
## はじめに
**(最初の方はReadmeの日本語原文です)**   
Quantum Espressoで計算するためにはinputファイルが必要で、これをある程度自動的に作成するpythonスクリプト群です。   
以下のような機能があります   
- cifからinputファイルを作成   
- relax, vc-relaxのoutputファイルから構造を読み取り、次のinputファイルを作成   
- relax, vc-relaxのoutputファイルから構造を読み取り、cifに変換   
- band, pdos などのplot   

想定しているユースケースとしては、「一つの物質について計算が収束していて、類似物質がたくさんあり、同条件で計算を回したい場合」です。   
最初の一つの物質でも使えなくはないですが、どちらかというと上手くいった条件を使い回すように作成されたスクリプトです。   
そのため、templateとして収束したinputファイルを入れて、物質固有の部分を変更するといった実装です。   

## インストール
基本的にpython標準モジュールと、有名なnumpy, pandas, matplotlibに加えmoyopyというものを使うだけなので、   
環境を汚すことはないと思うので、やらなくてもいいです。   
今回はMake_QE_fileを使用する用の環境(ここでは仮に`Env_Make_QE_file`とします)を用意します。   
特にどう作成しても問題ないですが   
```zsh
conda create -n "Env_Make_QE_file"
conda activate Env_Make_QE_file
```

とかでしょうか。   
以下でmoyopyなど必要なものも含めてインストールします。
```zsh
pip install git+https://github.com/souno1218/Make_QE_file.git
```

## 流れ
想定している使用の流れは   
1. relax -> vc-relax   
2. 1 -> scf -> nscf -> projwfc -> dos   
3. 1 -> scf -> bands -> band_x  
4. plot   

です(自分は数字ごとにディレクトリを分けてました)。   
`how_to_use/auto.py`では上記の流れでの各関数の使用方法を書いています。   

## 使用方法と機能

---
### cif_to_params(import_cif_path)
#### 概要
cifから結晶情報のdictを作成します。cifには以下の情報が必要です。   
- _cell_length_a
- _cell_length_b
- _cell_length_c
- _cell_angle_alpha
- _cell_angle_beta
- _cell_angle_gamma
- _space_group_IT_number(or _symmetry_Int_Tables_number)
- _atom_site_label
- _atom_site_type_symbol
- _atom_site_fract_x
- _atom_site_fract_y
- _atom_site_fract_z
#### Parameters:
- `import_cif_path` (str)   
  参照するcifのパス   
#### Returns:
- `params_cif` (dict)   
  結晶情報、これからinputを作成   

---
### make_input(args)
#### 概要
cif_to_paramsで作成したような結晶情報のdictからinputを作成します。   
構造は{crystal_sg}でinputファイルに記入されるため、`params_structure`にはspace_group_numberが必要になります。   
#### Parameters:
- `calc` (str)   
  対応している`calc`の値はrelax, vc-relax, scf, nscf, projwfc, dos, bands, band_xです(必須)   
  これ以外は自身がやったことがないためわかりません   
  calcによって必要な引数が違うので、それは下を参考に   
- `save_path` (str)   
  作成するinputファイルのパス(必須)
- `prefix` (str)   
  物質名など、判別可能な名前(必須)   
- `template_path=None` (str)   
  テンプレートのパス   
  何も指定しなければこちらが用意するテンプレが入ります   
- `pseudo_dir=None` (str)   
  calc = scf,nscf,relax,vc-relax,bandsで使用(必須)   
  pseudoファイルが入っているディレクトリのパス   
  ここから探してきてファイルパスとcutoffを指定します   
- `params_structure=None` (dict)   
  calc = scf,nscf,relax,vc-relax,bandsで使用(必須)   
  cif_to_paramsなどで作成した構造データを入れる   
- `k_fineness_Magnification=None`　(int)   
  calc = scf,nscf,relax,vc-relaxで使用   
  参考にしたサイト曰く「格子定数 x k点メッシュ数 = 10 ~ 12 Å 」という決め方に則って、この「10 ~ 12 Å」を指定します   
  何も指定しない場合大きめにとって20Åが入ります   
- `min_kpoint=None` (int)   
  calc = scf,nscf,relax,vc-relaxで使用   
  上記の際にあまりにもk点メッシュ数が少なく(1とか)になる時があるので、最低値を決めます   
  何も指定しない場合calc=nscfで4,その他で2になります   
- `fix=False` (bool)
  calc = scf,nscf,relax,vc-relaxで使用   
  原子を固定するかどうか。   
- `nbnd=None` (int)   
  calc = nscf, bandsで使用(必須)   
  バンド数で、scf計算で使用したnbndをoutput_to_nbndで読み取り、上記calcで使用   
- `brilloin_zone_path=None` (list(str))   
  calc = bandsで使用(必須)   
  Gamma点として"gG",Sigma点として"gS"が使用できる   
- `k_point_divisions=None` (list(int))   
  calc = bandsで使用(必須)   
  brilloin_zone_pathと同じだけ指定してください   
  つまり len(brilloin_zone_path) = len(k_point_divisions)
- `ibrav=None` (int)   
  calc = bandsで使用(必須?)   
  なんでかわからないですが、自分のQEのバージョンでは`calc=bands`では`ibrav`が必要でした   
- `emax=None` (float)   
  calc = projwfc, dosで使用   
  状態密度データのエネルギー範囲の上限らしい   
- `emin=None` (float)   
  calc = projwfc, dosで使用   
  状態密度データのエネルギー範囲の下限値らしい   
- `deltae=None` (float)   
  calc = projwfc, dosで使用   
  状態密度データのエネルギー刻み幅らしい   
   
##### emax, emin, deltaeについて
関数引数として与えられた場合それを使用。   
与えられなかった場合、テンプレートファイルの値を使用。   
テンプレートファイルにも値がない場合、以下のデフォルト値を使用。   
`emax, emin, deltae = 50, -50, 0.01`   

---
### output_to_params(calc, import_out_path, base_params)
#### 概要
relax,vc-relax計算のoutputで結晶情報を更新します。   
その際、計算前のparams(つまりcif_to_paramsの返り値)を元にして作るので`base_params`にはそれを入れます。   
#### Parameters:
- `calc` (str)   
  対応している`calc`の値はrelax, vc-relaxです   
- `import_out_path` (str)   
  outputファイルパス   
- `base_params` (dict)   
  元となった構造、最初のcif_to_paramsの返り値など   
#### Returns:
- `return_params` (dict)   
  結晶情報、これから次のinputを作成したり、write_new_cifに入れてcifを作ったりする   

---
### write_new_cif(out_path, material_name, params_structure)
#### 概要
output_to_paramsなどで更新した構造から新しいcifを作成します。   
#### Parameters:
- `out_path` (str)   
  出力するcifのoutputファイルパス   
- `material_name` (str)   
  物質名など
- `params_structure` (dict)   
  `output_to_params`などで更新した構造   

---
### output_to_nbnd(import_out_path)
#### 概要
scf計算のoutputから計算に使用したnbndを探して出力します。   
`number of Kohn-Sham states`が該当する部分です。   
この関数の返り値はcalc = nscf, bandsでinputファイルを作成するときに使用します。   
#### Parameters:
- `import_out_path` (str)   
  参照するoutputファイルパス   
#### Returns:
- `nbnd` (int)   
  scf計算で計算に使用したnbnd   

---
### get_EFermi(output_path)
#### 概要
計算結果のファイルからフェルミエネルギーを探して出力します。   
`the Fermi energy is`と`EFermi`が該当する部分です。   
この関数の返り値はplot_bandとplot_pdosで使用します。   
#### Parameters:
- `output_path` (str)   
  参照するoutputファイルパス   
#### Returns:
- `EFermi` (float)   
  EFermiの値   

---
### get_highest_occupied(output_path)
#### 概要
計算結果のファイルからhighest occupied levelを探して出力します。   
`highest occupied, lowest unoccupied level`と`highest occupied level`が該当する部分です。   
この関数の返り値はplot_bandとplot_pdosで使用します。   
#### Parameters:
- `output_path` (str)   
  参照するoutputファイルパス   
#### Returns:
- `highest_occupied` (float)   
  highest occupied levelの値   

---
### plot_band(args)
#### 概要
band計算によって作成されるgnuファイルを入れることで、band計算結果をプロットします。   
基準はhighest occupied levelかEFermiで、どちらかの入力は必須です。   
#### Parameters:
- `gnu_path` (str)   
  gnuファイルのパス
- `k_point_divisions` (list(int))   
  band計算のinput作成時に使用したものと同じものを使用   
- `brilloin_zone_path` (list(str))   
  band計算のinput作成時に使用したものと同じものを使用   
- `EFermi` (float)   
  基準として、グラフの0を定める   
  read_EFermiで探した値を入れる   
  EFermiかhighest_occupiedのどちらかは必要   
- `highest_occupied` (float)   
  基準として、グラフの0を定める   
  read_highest_occupiedで探した値を入れる   
  EFermiかhighest_occupiedのどちらかは必要   
- `is_save=False` (bool)   
  保存を行うかどうか   
- `is_plot=False` (bool)   
  plotするかどうか、対話型で使用時に指定   
- `savefig_path` (bool)   
  is_save=Trueのとき、保存するパス   
  拡張子はjpegやpngなど   
- `ylim=[-5, 5]` (list(num))   
  y軸のプロットする範囲   

---
### plot_pdos(args)
#### 概要
band計算を行ったディレクトリで実行することで、band計算結果をプロットします。   
基準はhighest occupied levelかEFermiで、どちらかの入力は必須です。   
#### Parameters:
- `pdos_dir_path` (str)   
  `{pdos_dir_path}/*.dos`や`{pdos_dir_path}/*_wfc*`となっているディレクトリのパス
- `plot_list=["pdos"]` (list(str))   
  何をプロットするか、   
  dos, pdos, tot_pdos, tot_pdos, tot_dosがある。
- `EFermi` (float)   
  基準として、グラフの0を定める   
  read_EFermiで探した値を入れる   
  EFermiかhighest_occupiedのどちらかは必要   
- `highest_occupied` (float)   
  基準として、グラフの0を定める   
  read_highest_occupiedで探した値を入れる   
  EFermiかhighest_occupiedのどちらかは必要   
- `is_save=False` (bool)   
  保存を行うかどうか   
- `is_plot=False` (bool)   
  plotするかどうか、対話型で使用時に指定   
- `savefig_path` (bool)   
  is_save=Trueのとき、保存するパス   
  拡張子はjpegやpngなど   
- `xlim=[-10, 10]` (list(num))   
  x軸のプロットする範囲、highest_occupiedかEFermiが0になる   
- `ylim=None` (list(num))   
  y軸のプロットする範囲   
- `color_dict` (dict(str:str))   
  元素ごとにまとまってプロットするため、元素ごとに色を変えることが可能。   
  設定しない場合黒になる。   


## 参考サイト
[公式Doc](https://www.quantum-espresso.org/Doc/INPUT_PW.html)   

#### (日本語)
[quantum ESPRESSO tutorial(東北大学物性理論研究室)](http://www.cmpt.phys.tohoku.ac.jp/~koretsune/SATL_qe_tutorial/)   

[雑多な記録(QuantumESPRESSO)](https://www2.yukawa.kyoto-u.ac.jp/~koudai.sugimoto/dokuwiki/doku.php?id=quantumespresso)   

[Quantum ESPRESSO入力ファイル作成手順](https://qiita.com/xa_member/items/727c1a62930611babaf7)   
