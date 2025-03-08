from .check_output import check_output
from .write_new_cif import write_new_cif
from .cif_to_params import cif_to_params
from .make_input import make_input
from .output_to_params import output_to_params
from .write_new_cif import write_new_cif
from .utils import check_JOB_DONE, check_output

import os
import importlib.resources as pkg_resources

templates_path = os.fspath(pkg_resources.path("Make_QE_file", "templates"))
print(f"templates_path : {templates_path}")

# https://zenn.dev/h_waka/articles/67d0dd45d8dc50
