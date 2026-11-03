# zsoil-output

Python SDK for reading ZSoil FEM result files: mesh (`C_Mesh`), solver history
(`C_HistoryOfExecution`, `C_Rcf_info`), and per-element-type results readers
(`C_NodalResults`, `C_BeamResults`, `C_TrussResults`, `C_ShellResults`,
`C_Contact_2D_Results`, `C_Contact_3D_Results`, `C_ContinuumResults`, ...).

Formerly shipped only as a folder inside the ZSoil installation
(`Tools v<year>/ZSoilPy3`) and added to `sys.path` by hand; now also
installable as a regular Python package.

## Install

```bash
pip install zsoil-output
```
or
```bash
pip install git+https://github.com/ZSoil-FEM/zsoil-output.git
```

## Usage

`zsoil_output` re-exports its main entry-point classes at the top level:

`Mesh`, `Element`, `HistoryOfExecution`, `RCF_info`, `RCF_layers_info`, `Diagram`,
`Diagram3D`, `EleResults`, `EleLayersResults`, `Beam_EleResults`, `Beam_S_EleResults`,
`Continuum_EleResults`, `Shell_EleResults`, `Truss_EleResults`, `NodalResults`,
`Contact_2D_Results`, `Contact_3D_Results`, `Heat_Exch_Results`, `ShellLayers_Results`,
`SuperConvPatchRecovery`.

`Element`'s `GROUP_*` constants (e.g. `Element.GROUP_BEAM`) are how you select which
elements to read results for.

A few supporting pieces aren't re-exported at the top level and are imported from
their own submodule instead — e.g. the `.xls` export helpers via `zsoil_output.ExcelUtl`,
or the fictitious-beam integration helpers (`postprocess_MNQ_time_history_for_fict_uniform_beam`
and friends) via `zsoil_output.C_ContinuumResults`.

```python
import zsoil_output
res_path = r"C:\path\to\MyZSoilResultsFiles"
his  = zsoil_output.HistoryOfExecution(res_path)
mesh = zsoil_output.Mesh(res_path)
rcf  = zsoil_output.RCF_info(res_path)
beam_rsl = zsoil_output.Beam_EleResults(mesh, his, rcf)
```

or, importing the names directly:

```python
from zsoil_output import HistoryOfExecution, Mesh, RCF_info, Beam_EleResults
res_path = r"C:\path\to\MyZSoilResultsFiles"
his, mesh, rcf = HistoryOfExecution(res_path), Mesh(res_path), RCF_info(res_path)
beam_rsl = Beam_EleResults(mesh, his, rcf)
```

### Using a local clone instead of the pip-installed package

To use a local clone (e.g. while developing the SDK itself) instead of, or
alongside, the pip-installed package, add the clone's `src` folder to `sys.path`
before importing — no install needed:

```python
import sys
sys.path.insert(0, r"C:\path\to\MyLocalPackageClone\src")   # the folder that CONTAINS zsoil_output\
import zsoil_output
```

See [zsoil-output-examples](https://github.com/ZSoil-FEM/zsoil-output-examples)
for worked examples reading real ZSoil projects with this SDK.
