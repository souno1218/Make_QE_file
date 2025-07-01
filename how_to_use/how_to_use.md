# Make_QE_file Usage Guide
## Introduction
**(The beginning is the original Japanese text of the Readme)**   
This is a collection of Python scripts designed to semi-automatically create input files for Quantum Espresso calculations.   
It provides the following functionalities:   
- Create input files from CIF files.   
- Read structures from relax and vc-relax output files and create subsequent input files.   
- Read structures from relax and vc-relax output files and convert them to CIF.   
- Plot band structures, PDOS, etc.   

The anticipated use case is for scenarios where "calculations for one material have converged, there are many similar materials, and you want to run calculations under the same conditions."   
While it's not impossible to use it for the initial material, these scripts are primarily designed to reuse successful conditions.   
Therefore, the implementation involves providing a converged input file as a template and then modifying the material-specific sections.   

## Installation
Since it primarily uses standard Python modules, along with popular libraries like NumPy, Pandas, Matplotlib, and Moyopy,   
it's unlikely to clutter your environment, so you might not strictly need to create a dedicated one.   
However, we will prepare an environment specifically for using Make_QE_file (here, we'll tentatively call it `Env_Make_QE_file`).   
You can create it in any way you prefer, for example:   
```zsh
conda create -n "Env_Make_QE_file"
conda activate Env_Make_QE_file
```
Then, install Moyopy and other necessary components as follows:   
```zsh
pip install git+[https://github.com/souno1218/Make_QE_file.git](https://github.com/souno1218/Make_QE_file.git)
```

## Workflow
The anticipated workflow is:   
1.  relax -> vc-relax   
2.  1 -> scf -> nscf -> projwfc -> dos   
3.  1 -> scf -> bands -> band_x   
4.  plot   
   
(I used to separate directories by number for these steps).   
`how_to_use/auto.py` contains instructions on how to use each function in the workflow described above.   

## About crystal_sg
For specifying atomic positions, we utilize ATOMIC_POSITIONS {crystal_sg}.   
According to the official documentation:   
> crystal_sg :
> atomic positions are in crystal coordinates, i.e. in relative coordinates of the primitive lattice.   
> This option differs from the previous one because in this case only the symmetry inequivalent atoms are given.   
> The variable space_group must indicate the space group number used to find the symmetry equivalent atoms.   
> The other variables that control this option are uniqueb, origin_choice, and rhombohedral.   

Therefore, for this option, the _space_group_IT_number (or _symmetry_Int_Tables_number) in the CIF file is required.   

## Encountered Issues
A list of issues encountered that are somewhat related, though not directly tied to this module.   

#### Memory Error with IntelMPI 2021.10 and QE 7.2
```text
= BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES
= RANK 31 PID 1078176 RUNNING AT pc~~~~~
= KILLED BY SIGNAL: 9 (Killed)
```
When a calculation terminates with this error, it often indicates a memory issue.   
Specifically, we observed continuous memory consumption growth when using the combination of **Intel MPI 2021.10.0** and **Quantum ESPRESSO 7.2**.   
This issue was not observed when using **Intel MPI 2021.10.0** with **Quantum ESPRESSO 6.8**.   

#### crystal_sg Fails to Expand Correctly
While performing calculations for a material with space_group = 129,    
we noticed that the output expanded into an unusually large cell.   
Upon inspection with VESTA, it was clearly incorrect.   
After some investigation, changing `origin_choice` from its default value of 1 to 2 resolved the issue.   
Therefore, `origin_choice = 2` has been included in the template.   


## Usage and Functions
-----
### cif_to_params(import_cif_path)
#### Overview
Creates a dictionary of crystal information from a CIF file. The CIF file must contain the following information:
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
    Path to the CIF file to be referenced.   
#### Returns:
  - `params_cif` (dict)   
    Crystal information, used to create the input file from.   

-----
### make_input(args)
#### Overview
Creates an input file from a crystal information dictionary, similar to one created by `cif_to_params`.   
The structure will be written into the input file as `{crystal_sg}`, so `params_structure` must include the `space_group_number`.   
#### Parameters:
  - `calc` (str)   
    Supported `calc` values are relax, vc-relax, scf, nscf, projwfc, dos, bands, band_x (required).   
    I'm unfamiliar with others as I haven't used them.   
    Required arguments vary by `calc` type; refer to the sections below.   
  - `save_path` (str)   
    Path to save the created input file (required).   
  - `prefix` (str)   
    Material name or other discernible identifier (required).   
  - `template_path=None` (str)   
    Path to the template.   
    If not specified, a default template provided by this tool will be used.   
  - `pseudo_dir=None` (str)   
    Used with `calc` = scf, nscf, relax, vc-relax, bands (required).   
    Path to the directory containing pseudo files.   
    The tool will search here for file paths and specify cutoff values.   
  - `params_structure=None` (dict)   
    Used with `calc` = scf, nscf, relax, vc-relax, bands (required).   
    Input the structural data created by `cif_to_params` or similar.   
  - `k_fineness_Magnification=None` (int)   
    Used with `calc` = scf, nscf, relax, vc-relax.   
    According to a referenced site, it follows the rule "Lattice constant x K-point mesh = 10 \~ 12 Å". This parameter specifies that "10 ~ 12 Å" value.   
    If not specified, a larger value of 20 Å will be used.   
  - `min_kpoint=None` (int)   
    Used with `calc` = scf, nscf, relax, vc-relax.   
    Sets a minimum value when the K-point mesh count becomes too small (e.g., 1).   
    If not specified, it defaults to 4 for `calc=nscf` and 2 for others.   
  - `fix=False` (bool)
    sed with `calc` = scf, nscf, relax, vc-relax.   
    Whether to fix atoms.   
  - `nbnd=None` (int)   
    Used with `calc` = nscf, bands (required).   
    Number of bands. The `nbnd` used in the SCF calculation is read via `output_to_nbnd` and then used for the above `calc` types.   
  - `brilloin_zone_path=None` (list(str))   
    Used with `calc` = bands (required).   
    "gG" can be used for Gamma point, "gS" for Sigma point.   
  - `k_point_divisions=None` (list(int))   
    Used with `calc` = bands (required).   
    Must specify the same number as `brilloin_zone_path`.   
    That is, `len(brilloin_zone_path) = len(k_point_divisions)`.   
  - `ibrav=None` (int)   
    Used with `calc` = bands (required?).   
    For some reason, `ibrav` was required for `calc=bands` in my QE version.   
  - `emax=None` (float)   
    Used with `calc` = projwfc, dos.   
    Presumably the upper limit of the energy range for density of states data.   
  - `emin=None` (float)   
    Used with `calc` = projwfc, dos.   
    Presumably the lower limit of the energy range for density of states data.   
  - `deltae=None` (float)   
    Used with `calc` = projwfc, dos.   
    Presumably the energy step size for density of states data.   

##### About emax, emin, deltae
If provided as function arguments, those values will be used.   
If not provided, values from the template file will be used.   
If values are also absent from the template file, the following default values will be used:   
`emax, emin, deltae = 50, -50, 0.01`   

-----
### output_to_params(calc, import_out_path, base_params)
#### Overview
Updates crystal information using the output of relax and vc-relax calculations.   
It uses the `params` from before the calculation (i.e., the return value of `cif_to_params`) as a base, so `base_params` should contain that.   
#### Parameters:
  - `calc` (str)   
    Supported `calc` values are relax, vc-relax.   
  - `import_out_path` (str)   
    Path to the output file.   
  - `base_params` (dict)   
    The original structure, such as the return value of the initial `cif_to_params`.   
#### Returns:
  - `return_params` (dict)   
    Crystal information, used to create the next input file or to be passed to `write_new_cif` to create a CIF.   

-----
### write_new_cif(out_path, material_name, params_structure)
#### Overview
Creates a new CIF from the structure updated by output_to_params or similar functions.   

#### Parameters:   
  - `out_path` (str)   
    Path to the output CIF file.   
  - `material_name` (str)   
    Material name, etc.   
  - `params_structure` (dict)   
    Structure updated by `output_to_params` or similar.   

-----
### output_to_nbnd(import_out_path)
#### Overview
Searches for and outputs the `nbnd` used in an SCF calculation from its output file.   
The relevant section is "number of Kohn-Sham states".   
The return value of this function is used when creating input files for `calc` = nscf, bands.   
#### Parameters:
  - `import_out_path` (str)   
    Path to the output file to be referenced.   
#### Returns:
  - `nbnd` (int)   
    The `nbnd` used in the SCF calculation.   

-----
### get_EFermi(output_path)
#### Overview
Searches for and outputs the Fermi energy from a calculation result file.   
The relevant sections are "the Fermi energy is" and "EFermi".   
The return value of this function is used in `plot_band` and `plot_pdos`.   
#### Parameters:
  - `output_path` (str)   
    Path to the output file to be referenced.   
#### Returns:
  - `EFermi` (float)   
    The value of EFermi.   

-----
### get_highest_occupied(output_path)
#### Overview
Searches for and outputs the highest occupied level from a calculation result file.   
The relevant sections are "highest occupied, lowest unoccupied level" and "highest occupied level".   
The return value of this function is used in `plot_band` and `plot_pdos`.   
#### Parameters:
  - `output_path` (str)   
    Path to the output file to be referenced.   
#### Returns:
  - `highest_occupied` (float)   
    The value of the highest occupied level.   

-----
### plot_band(args)
#### Overview
Plots band calculation results by taking gnu files created by band calculations.   
The reference point is either `highest_occupied_level` or `EFermi`; one of them must be provided.   
#### Parameters:
  - `gnu_path` (str)   
    Path to the gnu file.   
  - `k_point_divisions` (list(int))   
    Use the same values as those used when creating the band calculation input.   
  - `brilloin_zone_path` (list(str))   
    Use the same values as those used when creating the band calculation input.   
  - `EFermi` (float)   
    Sets the graph's 0 as a reference.   
    Input the value found by `read_EFermi`.   
    Either `EFermi` or `highest_occupied` is required.   
  - `highest_occupied` (float)   
    Sets the graph's 0 as a reference.   
    Input the value found by `read_highest_occupied`.   
    Either `EFermi` or `highest_occupied` is required.   
  - `is_save=False` (bool)   
    Whether to save the plot.   
  - `is_plot=False` (bool)   
    Whether to display the plot (specify when using interactively).   
  - `savefig_path` (str)   
    If `is_save=True`, the path to save the plot.   
    Extensions include jpeg, png, etc.   
  - `ylim=[-5, 5]` (list(num))   
    Plotting range for the y-axis.   

-----
### plot_pdos(args)
#### Overview
Plots band calculation results by executing in the directory where the band calculation was performed.   
The reference point is either `highest_occupied_level` or `EFermi`; one of them must be provided.   
#### Parameters:
  - `pdos_dir_path` (str)   
    Path to the directory containing `{pdos_dir_path}/*.dos` and `{pdos_dir_path}/*_wfc*` files.   
  - `plot_list=["pdos"]` (list(str))   
    What to plot,   
    Options include: dos, pdos, tot_pdos, tot_dos.   
  - `EFermi` (float)   
    Sets the graph's 0 as a reference.   
    Input the value found by `read_EFermi`.   
    Either `EFermi` or `highest_occupied` is required.   
  - `highest_occupied` (float)   
    Sets the graph's 0 as a reference.   
    Input the value found by `read_highest_occupied`.   
    Either `EFermi` or `highest_occupied` is required.   
  - `is_save=False` (bool)   
    Whether to save the plot.   
  - `is_plot=False` (bool)   
    Whether to display the plot (specify when using interactively).   
  - `savefig_path` (str)   
    If `is_save=True`, the path to save the plot.   
    Extensions include jpeg, png, etc.   
  - `xlim=[-10, 10]` (list(num))   
    Plotting range for the x-axis, where `highest_occupied` or `EFermi` becomes 0.   
  - `ylim=None` (list(num))   
    Plotting range for the y-axis.   
  - `color_dict` (dict(str:str))   
    Allows changing colors for each element, as elements are plotted together.   
    If not set, it defaults to black.   

## References
[Official Doc](https://www.quantum-espresso.org/Doc/INPUT_PW.html)
#### (Japanese)
[quantum ESPRESSO tutorial(東北大学物性理論研究室)](http://www.cmpt.phys.tohoku.ac.jp/~koretsune/SATL_qe_tutorial/)   

[雑多な記録(QuantumESPRESSO)](https://www2.yukawa.kyoto-u.ac.jp/~koudai.sugimoto/dokuwiki/doku.php?id=quantumespresso)   

[Quantum ESPRESSO入力ファイル作成手順](https://qiita.com/xa_member/items/727c1a62930611babaf7)   
