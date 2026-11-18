# Lattice VinaPilotX — details

Lattice VinaPilotX runs **AutoDock Vina** from inside UCSF ChimeraX. It takes a prepared receptor and prepared ligands (PDBQT), lets you define the search box on the structure, runs single or batch docking, and collects the results. This page explains every feature and how it works. (Working name; the tool is not yet final-named.)

> Inputs for docking are `.pdbqt` files. The **Prepare** tab can make them from PDB/mmCIF (receptor) and SDF/MOL/MOL2/SMILES (ligands) with Meeko.

## Highlights

- One window, six tabs: **Settings → Prepare → Inputs → Box → Run → Results**.
- **Meeko preparation** in a private environment (micromamba, with a venv + pip fallback) set up from the Settings tab; receptor cleaning is done in ChimeraX.
- Finds your Vina, or **downloads the official Vina** for your OS/CPU and verifies it.
- Receptor opens in the scene when selected; other models are hidden to avoid clashes.
- **Ligand table** with size and chemistry info per ligand; the **largest ligand is highlighted**; export to CSV.
- Single ligand, a folder of ligands, or a **multi-ligand library file** (split automatically).
- Search-box editor: numeric entry, from selection, mouse move/resize, undo/redo.
- **Batch docking** with parallel jobs, per-ligand progress bars, cancel, and **resume**.
- Reproducible: Vina's reported **seed**, exact command, versions and parameters are saved with every run.
- Expandable per-pose results, ratings, per-pose H-bonds/clashes, instant pose display with a capped cache, saved analysis, CSV export.
- Optional switch that makes ChimeraX open `.pdbqt` files natively (fixes the "PDBT / Vinardock" error).
- Everything is also scriptable with `dock ...` commands.

---

## 1. Settings tab

### AutoDock Vina program
Lattice VinaPilotX does not contain Vina and installs nothing into ChimeraX; it starts Vina as a separate process.

- **Browse / Verify**: you can point at `vina.exe` (or `vina`) or at the folder containing it. If you give a folder, files named `vina*` are considered (never `vina_split*`), newest version first. The file is test-run with `--version`; the reported version is shown. The choice is remembered.
- **Download latest Vina**: looks up the newest release through the GitHub API (falls back to the `/releases/latest` redirect, then to a built-in version number). It picks the asset for your platform — `vina_<ver>_win.exe`, `vina_<ver>_mac_x86_64`, `vina_<ver>_mac_aarch64`, `vina_<ver>_linux_x86_64`, `vina_<ver>_linux_aarch64` — downloads it to a `.part` file, checks the SHA-256 when GitHub supplies one, makes it executable, test-runs `--version`, and only then keeps it. Failed or cancelled downloads leave no partial file. Downloads are stored in ChimeraX's user-data folder under `VinaPilotX/vina/`.
- On startup, a remembered path is used; otherwise Lattice VinaPilotX looks for `vina` on `PATH` and for earlier downloads.

### Working directory
All outputs of a run go here (layout in section 6). The folder is created if missing.

### Opening PDBQT files automatically
Some ChimeraX installs associate the `.pdbqt` suffix with a format called **PDBT**. Plain `open file.pdbqt` (menu, drag-and-drop, double-click) then fails with *"Install the Vinardock bundle to open PDBT format files"*, and the file can only be opened through ViewDock or with `open file.pdbqt format pdbqt`.

The cause: ChimeraX's data-format manager consults each format's `default_for` suffix list **before** the user's saved suffix preference, so changing the preference is not enough. When this box is ticked, Lattice VinaPilotX:

1. adds `.pdbqt` to the `default_for` list of ChimeraX's built-in **AutoDock PDBQT** format,
2. removes only the `.pdbqt` claim from the PDBT format (PDBT keeps its own `.pdbt` suffix),
3. saves the suffix preference and clears any cached format choice.

Multi-pose Vina output then opens in ViewDock and single structures open as ordinary models. The original state is remembered while the override is active, so **unticking restores it exactly**. The setting is stored and re-applied every time ChimeraX starts (and once more just after startup, in case another bundle registers formats late). Uncheck it if it ever conflicts with another bundle. The command `vinax pdbqt native true|false` does the same as the checkbox.

---

### Preparation environment (Meeko)
Meeko and RDKit are installed in a private environment, never into ChimeraX. Every Meeko command runs as a separate process.

- **Set up Meeko** with *Automatic*: downloads the official **micromamba** program (GitHub releases, then micro.mamba.pm), creates a conda-forge environment (`python 3.11, meeko, rdkit, scipy, numpy, gemmi, prody`) and checks it. If that fails, it falls back to a **venv + pip** environment (system Python first, then the Python shipped with ChimeraX; ProDy is optional there). You can also force one of the two methods.
- **Use existing Python…** uses a Python that already has Meeko and RDKit (a conda or venv Python; nothing is installed).
- **Verify** re-checks the environment (Meeko version, RDKit, ProDy). **Remove** deletes an environment this tool created.
- The environment lives in ChimeraX's user-data folder under `VinaPilotX/prep/` and is remembered between sessions. Progress and pip/micromamba output are shown in the log box.
- Meeko's command-line options changed between versions, so the tool reads `--help` of the installed Meeko and adapts (e.g. `--read_with_prody` / `--read_pdb` / `-i`, `--delete_bad_res` / `--allow_bad_res`).

---

## 2. Prepare tab (receptor and ligands)

### Receptor (PDB / mmCIF → PDBQT): three steps
Pick a file or an already open model (a copy is used; your model is never changed).

1. **Clean receptor** opens a copy as `<name>_prep` and keeps it in the scene so you can inspect it. It uses ChimeraX **Dock Prep** (`dockprep`) with your choices: delete water, delete ions, delete alternate locations, incomplete side chains (complete with Dunbrack / Dynameomics / Richardson rotamers, truncate to Ala or Gly, or leave), and hydrogens; chains and "other ligands / cofactors" (everything that is not protein, nucleic acid, water or ions) are removed with `delete`. Add charges and mol2 output are not used, because Meeko computes its own charges.
   - **Hydrogens: Add if missing (recommended)** adds hydrogens only when the structure has none, so protonation you set earlier (hydrogens, `HID`/`HIE`/`HIP`, `ASH`, `GLH`, `LYN`, `CYM` names) is never changed. *Always add* runs ChimeraX `addh` on every run; *Don't add* leaves the structure as it is.
   - Meeko does **not** need hydrogens in its input: it adds the polar hydrogens from the residue names itself (histidines get its default state unless named `HID`/`HIE`/`HIP`). ChimeraX's hydrogens are optional but pick histidine states from the surroundings.
2. **Minimize (optional)** runs ChimeraX `minimize … dockPrep false maxSteps N` on the copy (default 100 steps; "until convergence" is slow and not recommended). It runs in the background: wait until the Log shows it has finished. Dock Prep is not repeated, so hydrogens and protonation stay as they are. Standard residues use the force field's own charges, so no extra charges are needed unless non-standard residues or cofactors remain.
3. **Convert with Meeko** saves the copy as `prepared/receptor/<name>_clean.pdb`, runs `mk_prepare_receptor -o <name> -p` and writes `prepared/receptor/<name>.pdbqt`. If a residue matches no Meeko template, it is retried **without its hydrogens** (Meeko then adds standard ones); residues that still fail are skipped (*Skip residues Meeko cannot handle*) and listed, because a skipped residue is missing from the receptor. The result is placed in the receptor field of the Inputs tab.

### Ligands (SDF / MOL / MOL2 / SMILES → PDBQT)
1. Add files or folders (sub-folders included). `.smi` files are lines of `SMILES name`.
2. **RDKit staging**: molecules without hydrogens, or with flat 2D coordinates, get hydrogens and a 3D conformer (ETKDG + MMFF); molecules that already have 3D coordinates and hydrogens are used as they are. Protonation states are **not** changed. Names come from the SDF title/`_Name` or the file name; duplicate names get `_2`, `_3` …; unreadable molecules are listed as skipped.
3. **Meeko**: `mk_prepare_ligand` writes one PDBQT per molecule to `prepared/ligands/` (Gasteiger charges, Meeko's defaults).
4. The PDBQT files are added to the Inputs tab automatically.

---

## 3. Inputs tab

### Receptor
Select a prepared receptor `.pdbqt`. It is checked and then opened in the scene:

- **Checks**: has ATOM/HETATM records; atom types are valid AutoDock types (error otherwise); no `ROOT/TORSDOF` records (warning: looks like a ligand or flexible receptor); a single model; presence of polar hydrogens (`HD`; warning if none).
- **Hiding**: other top-level models (an unrelated protein, earlier results) are hidden so they do not clash with the docking view. The search-box model is left alone. Use the ChimeraX Models panel to show them again.
- **Box presets**: the Box tab's center/size fields are pre-filled with the receptor's centre and extent (a blind-docking starting point). **No box is drawn** until you create one.
- Reselecting the same file does not reopen it; selecting another receptor closes the previous one.

### Ligands
Add files, add a folder (optionally with sub-folders), or add library files that contain many ligands.

**How a ligand is recognised.** A ligand block ends at its `TORSDOF` record. `MODEL/ENDMDL` wrapper lines are dropped. A file with one block is a single ligand; a file with several blocks is a library and every block becomes one ligand. Names come from `REMARK Name =` when present (library entries) or the file name (plus `_n` for unnamed library entries); names are made filesystem-safe and unique. Adding the same file twice is ignored. A file with no ligand records, or a block missing `ROOT`, `TORSDOF`, atoms or valid atom types, is shown in red and skipped at run time.

**Table columns**

| Column | Meaning and how it is computed |
|---|---|
| Ligand | Name as above (SMILES shown on hover if the file has a `REMARK SMILES` line). |
| Heavy atoms | Atoms that are not `H`, `HD` or `HS`. |
| Torsions | The `TORSDOF` value (rotatable bonds). |
| Dimensions (Å) | Extents of the heavy-atom centres along the molecule's **principal axes** (from an SVD of the centred coordinates), largest first. Independent of how the file is oriented. These describe the conformation stored in the file; flexible ligands can extend further. Atom radii are not added. |
| Longest (Å) | The largest distance between any two heavy atoms. |
| Charge | Sum of the PDBQT partial-charge column. |
| Acceptors/Donors | Count of acceptor atoms (`OA`, `NA`, `SA`) / count of polar hydrogens (`HD`). |
| Folder, File | Where the ligand came from; library entries show `[#n]`. |

**Largest ligand.** The row with the largest *longest distance* is highlighted, and a note lists its name, dimensions and longest distance. The note ends with a *suggestion*, not a rule: the box should exceed these sizes in all directions so the ligand can move freely. With one very large ligand in a mixed set, or a very asymmetric ligand, sizing the box from the largest ligand may not be what you want.

**Save table (CSV)** writes the table (in the order currently shown) with extra columns: separate dimensions, SMILES, index in file, validity and problems.

---

## 4. Box tab

The search space is Vina's `center_x/y/z` and `size_x/y/z` in Å. Center and size are two columns (X, Y, Z rows); if the panel is too narrow it scrolls sideways.

- **New box** draws a box from the numbers. Nothing is drawn until you ask.
- **Create box from selection** uses the current atom selection:
  - *Selection bounds + padding* (recommended): per axis, `min − p` to `max + p`; centre is the midpoint, size is the difference. No rounding.
  - *Vina ligand-centroid (legacy)*: centre = mean of the selected coordinates; per axis `size = ceil(2 × (max |q − centre| + p))`. This mirrors how Vina derives a box from a ligand. It is not recommended for receptor pockets.
  - Changing padding or method recalculates the box from the saved selection until you edit the numbers by hand.
- **Mouse interaction**: left-drag moves the box, right-drag resizes the face under the cursor (opposite face fixed). It temporarily replaces the standard left/right mouse modes and restores them when unticked or when the tool closes. Each gesture is one undo step.
- **Rounded grid span**: Vina discretises the box into voxels; `n = ceil(size / spacing)`, span `= n × spacing` (default spacing 0.375 Å). The optional outline shows this; it does not change the exported size.
- **Check selection** reports how many selected atoms are inside the box and the maximum overhang. **Copy Vina config** copies the six values. **View box**, **Clear**, colour and "draw on top" are display helpers.
- Everything (create, edit, mouse gestures, colour, clear) is undoable with ChimeraX's undo/redo.
- Coordinates are ChimeraX **scene** coordinates. If you moved a structure in the scene without changing its saved coordinates, the box and the files Vina reads can disagree.

The box editor is derived from the standalone ChimeraX Vina Grid Box bundle (same author); the mouse code adapts LGPL code from ChimeraX-Vinardock (see `THIRD_PARTY_NOTICES.md`).

---

## 5. Run tab

### Settings
| Setting | Default | Notes |
|---|---|---|
| Exhaustiveness | 8 | Search effort. Higher is slower and more thorough. |
| Number of poses | 9 | Maximum poses written per ligand. |
| Energy range | 3 kcal/mol | Poses worse than the best by more than this are dropped. |
| Scoring function | vina | `vina` or `vinardo`. |
| CPU cores per job | automatic | 0 = Vina detects and uses all available cores. |
| Parallel jobs | 1 | Number of docking jobs running at once. |
| Random seed | random | An integer makes runs reproducible. |

**Reset to defaults** restores all of these.

**CPU cores and parallel jobs.** Each ligand is an independent Vina process. *Parallel jobs* is how many such processes run at the same time. If CPU cores per job is "automatic" and parallel jobs is greater than 1, the machine's cores are divided between the jobs (`cores ÷ jobs`, at least 1 each); if you set cores per job explicitly, that number is used for each job. Parallel jobs are **not** multi-ligand co-docking: Vina can dock several ligands simultaneously into one site in a single run, but Lattice VinaPilotX does not do that yet; here every ligand is docked alone.

**Seeds.** By default Vina chooses its own random seed for every ligand. Lattice VinaPilotX reads the seed Vina reports in its output (`random seed: N`) and records it. If you enter an integer, `--seed N` is passed for every ligand. Seeds are shown in the Results table and saved in `summary.csv`, `manifest.json` and each ligand's log.

### What happens when you press Run
1. **Pre-flight checks**: Vina runs; receptor is valid; a box exists; parameters are in range; working folder set; warnings for boxes above 27,000 Å³ (Vina's own warning threshold) and for non-ASCII paths on Windows.
2. **Staging**: the receptor is copied to `receptor/`; each chosen ligand is written to `ligands/<name>.pdbqt` (libraries are split; individual ligands you removed from the table are not written); the box and search parameters are written to `vina_config.txt`.
3. **Jobs**: one Vina process per ligand, launched as
   `vina --config vina_config.txt --receptor <receptor> --ligand <ligand> --out results/<name>_out.pdbqt [--cpu N] [--seed N]`.
   File paths are passed as separate arguments (never written into the config), so folders with spaces are safe.
4. **Progress**: Vina prints a bar of 51 `*` characters; Lattice VinaPilotX counts them to fill that ligand's progress bar (one bar per running ligand). A counter shows *Ligands completed: x / y*. Whether the bar moves smoothly depends on Vina flushing its output.
5. **Console output** of every ligand is saved to `logs/<name>.log`; the first line records the exact command.
6. **Cancel** stops running Vina processes and marks unfinished ligands as cancelled.
7. A failing ligand (Vina error, no poses) is marked *failed* with the last useful error lines; it does not stop the batch.

### Resume
If a ligand already has a result and a log in the working folder, it is reused (status "done (existing)") instead of re-docked, so an interrupted batch continues where it stopped. Failed ligands are retried. Tick *Re-dock ligands that already have results* to force everything to run again. Note that resume does not check whether the receptor, box or parameters changed; use a new working folder when they do.

---

## 6. Results tab

- **Expandable tree.** One row per ligand (best affinity, number of poses, seed, status, best rating). Expand it to list every pose: affinity, RMSD l.b., RMSD u.b., rating, H-bonds, clashes. Poses stay in pose order when you sort; ligand rows sort by any column.
- **Instant display.** A ligand's whole Vina output (all poses) is opened once as one model, so its poses appear as sub-models (e.g. #3.1, #3.2 …). Selecting a pose shows that sub-model and hides the others, so clicking through poses is instant. At most *Max ligand files kept open* (default 5) stay loaded; opening another closes the least recently used one. The receptor is opened once and reused.
- **Rating.** A 1–5 star box per pose, stored in `analysis.json`. **Save rated poses…** copies poses at or above a chosen rating, plus a CSV and the receptor, into `selected_poses/`.
- **H-bonds / Clashes.** Select pose(s) and click (up to 3 per click). H-bonds come from one ChimeraX `hbonds` run with the "color H-bonds not meeting precise criteria differently" option (`twoColors`/`slopColor`); bonds that only meet the relaxed criteria get their own colour and show as "3 (+1)" (3 strict, 1 relaxed). Clashes use `clashes`. All bonds of all poses live in just two models, "Lattice H-bonds" and "Lattice clashes", and are shown only for the displayed pose. The receptor residues revealed around a pose's bonds are hidden again when you switch away from that pose.
- **Saved interactions.** Each pose's bonds are stored in `pose_cache/<ligand>_pose<N>_interactions.json` (atom indexes) and the counts in `analysis.json`. Showing the pose again, also after restarting ChimeraX and using *Load folder…*, redraws them without recalculating. Saved bonds are skipped if the receptor no longer matches.
- **Colours.** Settings tab → Interaction colours: H-bonds (default #00bfff), H-bonds with relaxed criteria (default #2f6fd6) and clashes (default #ff6a00). Changes apply immediately to drawn bonds and are remembered.
- **Open in ViewDock** runs `viewdock` on the already-open model of the selected ligand(s).
- **Export CSV** writes every pose with rating, H-bonds, clashes and partners. **Load folder…** rebuilds the tree from a working folder, including ratings and saved H-bond/clash results, so analysis continues where you stopped.
- The Run-tab log is saved to `run_log.txt` in the working folder.

---

## 7. Files written to the working directory

```
prepared/receptor/         cleaned PDB and Meeko receptor PDBQT (Prepare tab)
prepared/ligands/          Meeko ligand PDBQT files; prepared/staging/ holds the RDKit SDFs
receptor/<name>.pdbqt      copy of the receptor used
ligands/<name>.pdbqt       one file per ligand (libraries split)
results/<name>_out.pdbqt   Vina poses (REMARK VINA RESULT lines give the affinities)
logs/<name>.log            command line + Vina console output
vina_config.txt            box and search parameters shared by all ligands
summary.csv                one row per POSE: ligand, pose, affinity, rmsd_lb, rmsd_ub, rating,
                           hbonds, hbonds_relaxed, clashes, hbond_partners, seed, status, error
best_per_ligand.csv        one row per ligand (best affinity, poses, seed, status)
run_log.txt                timestamped copy of the Run-tab log
analysis.json              ratings, H-bond and clash results
pose_cache/                saved H-bond/clash atom pairs per pose (deleting it loses only the drawing; counts stay in analysis.json)
selected_poses/            written by "Save rated poses…"
manifest.json              tool/Vina versions, receptor, box, parameters, seed policy and seeds,
                           an example command, status counts, Python/OS
```

---

## Known limits
- Docking inputs must be PDBQT; the Prepare tab creates them. Ligand protonation states are not assigned; flexible/covalent receptor options of Meeko are not exposed yet.
- Resume does not detect changed receptor/box/parameters.
- Ligand dimensions describe the stored conformation only.
- Progress-bar smoothness depends on Vina's output buffering.
- Tested here with a stand-in Vina and a simulated ChimeraX; test with your own Vina and ChimeraX before relying on results.

## Cite
Cite AutoDock Vina when you publish results: Trott & Olson, *J. Comput. Chem.* 2010, 31, 455; Eberhardt et al., *J. Chem. Inf. Model.* 2021, 61, 3891. Cite ChimeraX as requested by UCSF.

---

## Commands

```
vinax gui [settings|inputs|box|run|results]     open the window
vinax vina [path FILE|FOLDER] [download true]   show, set or download Vina
vinax box center x,y,z size x,y,z               create/update the box from numbers
vinax box selection true [padding 4] [method bounds|centroid]
vinax run RECEPTOR ligands PATH workdir DIR
         [center x,y,z size x,y,z]              (default: the box in the scene)
         [exhaustiveness 8] [num_modes 9] [energy_range 3] [cpu 0] [parallel 1]
         [seed 0] [scoring vina|vinardo] [overwrite true] [recursive true]
         PATH = a file, a folder, a library file, or a wildcard such as C:/ligs/*.pdbqt
vinax cancel                                    cancel the running job
vinax pdbqt [native true|false]                 show/set the native .pdbqt opening override
```
