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

## Contributing

Contributions go through a fork and a pull request — no write access to this
repository is needed.

1. **Fork** — on [github.com/ZSoil-FEM/zsoil-output](https://github.com/ZSoil-FEM/zsoil-output),
   click **Fork** (top right). This creates `github.com/<your-user>/zsoil-output`.

2. **Clone your fork** and register this repository as `upstream`:
   ```bash
   git clone https://github.com/<your-user>/zsoil-output.git
   cd zsoil-output
   git remote add upstream https://github.com/ZSoil-FEM/zsoil-output.git
   ```

3. **Install in editable mode** — Python then imports `zsoil_output` straight from
   your clone's `src` folder, so your edits take effect without reinstalling:
   ```bash
   pip uninstall zsoil-output     # only if a regular install is already present
   pip install -e .
   ```
   Check that the clone is the one being used:
   ```bash
   python -c "import zsoil_output; print(zsoil_output.__file__)"
   ```
   The printed path should point into your clone's `src\zsoil_output` folder.
   Use the same Python / conda environment you run your ZSoil scripts with.

4. **Create a branch** for each change, starting from the latest upstream `main`:
   ```bash
   git fetch upstream
   git checkout -b my-fix upstream/main
   ```

5. **Make your changes, test, commit and push** to your fork:
   ```bash
   git add -A
   git commit -m "Short description of the change"
   git push -u origin my-fix
   ```

6. **Open a pull request** — GitHub shows a *Compare & pull request* button on
   your fork after the push. Target `ZSoil-FEM/zsoil-output`, branch `main`, and
   describe what changed and how you tested it.

To keep your fork up to date later:
```bash
git fetch upstream
git checkout main
git merge upstream/main
git push origin main
```

Please use `git` (not the GitHub web *Upload files* button) to submit changes:
uploads replace whole files, including their line endings, which makes the
pull request very hard to review.

See [zsoil-output-examples](https://github.com/ZSoil-FEM/zsoil-output-examples)
for worked examples reading real ZSoil projects with this SDK.
