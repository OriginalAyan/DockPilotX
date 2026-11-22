from chimerax.core.toolshed import BundleAPI

TOOL_NAME = "Lattice VinaPilotX"


class _VinaPilotXAPI(BundleAPI):
    api_version = 1

    @staticmethod
    def initialize(session, bundle_info):
        """Apply the (user-switchable) native .pdbqt reader override at startup."""
        from .controller import Settings
        try:
            wanted = bool(Settings().get("native_pdbqt", True))
        except Exception:
            wanted = True
        if not wanted:
            return
        from . import format_preference as fp

        def apply(log_missing=True):
            try:
                changed, _prev = fp.prefer_native_pdbqt(session.data_formats, force=True)
            except LookupError as err:
                if log_missing:
                    session.logger.warning("Lattice VinaPilotX could not set the native .pdbqt default: %s" % err)
                return False
            if changed:
                session.logger.info("Lattice VinaPilotX: .pdbqt files now open with ChimeraX's native AutoDock "
                                    "PDBQT reader (switch off in Lattice VinaPilotX > Settings).")
            return True

        apply(log_missing=False)
        try:
            from Qt.QtCore import QTimer
            QTimer.singleShot(0, apply)  # once more after other bundles have registered their formats
        except Exception:
            pass

    @staticmethod
    def start_tool(session, bi, ti):
        if ti.name != TOOL_NAME:
            raise ValueError("Unknown tool name: %s" % ti.name)
        from chimerax.core import tools
        from .tool import VinaPilotXTool
        return tools.get_singleton(session, VinaPilotXTool, ti.name, create=True)

    @staticmethod
    def register_command(bi, ci, logger):
        from chimerax.core.commands import register
        from .commands import COMMANDS
        try:
            func, desc = COMMANDS[ci.name]
        except KeyError:
            raise ValueError("Unknown command name: %s" % ci.name)
        if desc.synopsis is None:
            desc.synopsis = ci.synopsis
        register(ci.name, desc, func, logger=logger)

    @staticmethod
    def run_provider(session, name, mgr, **kw):
        """Opens ``.vinapx`` project files (File > Open, drag and drop, or ``open project.vinapx``)."""
        if mgr is session.open_command:
            from chimerax.open_command import OpenerInfo

            class _ProjectOpener(OpenerInfo):
                def open(self, session, data, file_name, **kw):
                    from .commands import cmd_open
                    cmd_open(session, str(data))
                    return [], "Opened the Lattice VinaPilotX project %s" % file_name
            return _ProjectOpener()
        raise ValueError("Unknown provider: %s" % name)

    @staticmethod
    def get_class(class_name):
        if class_name == "VinaPilotXTool":
            from .tool import VinaPilotXTool
            return VinaPilotXTool
        raise ValueError("Unknown class name '%s'" % class_name)


bundle_api = _VinaPilotXAPI()
