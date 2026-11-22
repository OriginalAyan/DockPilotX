"""The ``vinapilotx`` command family."""

import os

from chimerax.core.commands import (
    BoolArg, CmdDesc, EnumOf, FloatArg, Float3Arg, IntArg, StringArg,
)
from chimerax.core.errors import UserError

from .controller import get_controller
from .docking import DockingParams, SCORING_FUNCTIONS


TABS = ["settings", "protein", "ligands", "inputs", "box", "run", "results"]


def _tool(session, show=True):
    from chimerax.core import tools
    from .tool import VinaPilotXTool
    tool = tools.get_singleton(session, VinaPilotXTool, "Lattice VinaPilotX", create=True)
    if tool is not None and show:
        tool.tool_window.shown = True
    return tool


def cmd_gui(session, tab=None):
    tool = _tool(session)
    if tab is not None:
        tool.show_tab(TABS.index(tab))


gui_desc = CmdDesc(optional=[("tab", EnumOf(TABS))],
                   synopsis="open the Lattice VinaPilotX window")


def cmd_open(session, path):
    """Open a project.vinapx file, or a folder that contains one."""
    from .project import Project
    path = os.path.expanduser(path)
    if Project.locate(path) is None:
        raise UserError("No project.vinapx found in %s" % path)
    tool = _tool(session)
    try:
        proj = tool.open_project(path)
    except Exception as err:
        raise UserError("Could not open the project: %s" % err)
    session.logger.info("Opened the Lattice VinaPilotX project in %s" % proj.workdir)


open_desc = CmdDesc(required=[("path", StringArg)], synopsis="open a Lattice VinaPilotX project (.vinapx)")


def cmd_vina(session, path=None, download=False):
    ctrl = get_controller(session)
    if download:
        ctrl.download_vina()
        session.logger.info("Downloading the latest AutoDock Vina in the background ...")
        return
    if path:
        ok, info = ctrl.set_vina(os.path.expanduser(path))
        if not ok:
            raise UserError(info)
        return
    if ctrl.vina_exe:
        session.logger.info("AutoDock Vina %s: %s" % (ctrl.vina_version, ctrl.vina_exe))
    else:
        session.logger.info("AutoDock Vina is not set up. Use 'vinax vina path <file>' "
                            "or 'vinax vina download true'.")


vina_desc = CmdDesc(optional=[("path", StringArg)], keyword=[("download", BoolArg)],
                    synopsis="set, verify or download the AutoDock Vina program")


def cmd_box(session, center=None, size=None, selection=False, padding=4.0,
            method="bounds"):
    tool = _tool(session, show=False)
    panel = tool.box_panel
    if selection:
        message = panel.from_selection(method="bounds" if method == "bounds" else "vina", padding=padding)
        if message:
            session.logger.info(message)
    elif center is not None and size is not None:
        panel.set_box_from_cli(center, size, panel.spacing.value())
    elif center is not None or size is not None:
        raise UserError("Give both center and size, or use 'selection true'.")
    box = panel.current_box()
    if box is None:
        session.logger.info("No search box.")
    else:
        c, s = box
        session.logger.info("Search box: center %.3f,%.3f,%.3f  size %.3f,%.3f,%.3f" % (*c, *s))


box_desc = CmdDesc(
    keyword=[("center", Float3Arg), ("size", Float3Arg), ("selection", BoolArg), ("padding", FloatArg),
             ("method", EnumOf(["bounds", "centroid"]))],
    synopsis="create/update the search box from numbers or the current selection")


def cmd_dock(session, receptor, ligands, workdir, center=None, size=None, exhaustiveness=8,
             num_modes=9, energy_range=3.0, cpu=0, parallel=1, seed=0, scoring="vina",
             overwrite=False, recursive=False):
    ctrl = get_controller(session)
    if (center is None) != (size is None):
        raise UserError("Give both center and size, or neither to use the box in the scene.")
    if center is None:
        from .undo_support import find_box_model
        box = find_box_model(session)
        if box is None:
            raise UserError("No search box in the scene. Create one with 'vinax box ...' "
                            "or give center and size.")
        center, size = box.center, box.size
    params = DockingParams(exhaustiveness=exhaustiveness, num_modes=num_modes, energy_range=energy_range,
                           cpu=cpu, parallel_jobs=parallel, seed=seed, scoring=scoring)
    import glob
    lig_path = os.path.expanduser(ligands)
    ligand_inputs = sorted(glob.glob(lig_path)) if any(c in lig_path for c in "*?[") else [lig_path]
    if not ligand_inputs:
        raise UserError("No ligand files match: %s" % ligands)
    ok, result = ctrl.start_docking(
        os.path.expanduser(receptor), ligand_inputs,
        os.path.expanduser(workdir), tuple(center), tuple(size), params,
        recursive=recursive, overwrite=overwrite)
    if not ok:
        raise UserError(str(result))
    session.logger.info("Docking started for %d ligand(s). Progress appears in the log and the "
                        "Lattice VinaPilotX window." % result)


dock_desc = CmdDesc(
    required=[("receptor", StringArg)],
    keyword=[("ligands", StringArg), ("workdir", StringArg), ("center", Float3Arg), ("size", Float3Arg),
             ("exhaustiveness", IntArg), ("num_modes", IntArg), ("energy_range", FloatArg), ("cpu", IntArg),
             ("parallel", IntArg), ("seed", IntArg), ("scoring", EnumOf(list(SCORING_FUNCTIONS))),
             ("overwrite", BoolArg), ("recursive", BoolArg)],
    required_arguments=["ligands", "workdir"],
    synopsis="dock a ligand file, a folder of ligands, a wildcard set, or a ligand library with AutoDock Vina")


def cmd_cancel(session):
    get_controller(session).cancel_docking()


cancel_desc = CmdDesc(synopsis="cancel the running docking job")


def cmd_pdbqt(session, native=None):
    ctrl = get_controller(session)
    if native is None:
        session.logger.info("Native .pdbqt opening is %s." % ("ON" if ctrl.native_pdbqt_wanted else "OFF"))
        return
    if not ctrl.set_native_pdbqt(native):
        raise UserError("Could not change the .pdbqt setting; see the log.")


pdbqt_desc = CmdDesc(keyword=[("native", BoolArg)],
                     synopsis="make ChimeraX's native AutoDock PDBQT reader the default for .pdbqt files (or undo that)")


COMMANDS = {
    "vinax gui": (cmd_gui, gui_desc),
    "vinax open": (cmd_open, open_desc),
    "vinax vina": (cmd_vina, vina_desc),
    "vinax box": (cmd_box, box_desc),
    "vinax run": (cmd_dock, dock_desc),
    "vinax cancel": (cmd_cancel, cancel_desc),
    "vinax pdbqt": (cmd_pdbqt, pdbqt_desc),
}
