"""Undo/redo for receptor preparation steps (hydrogens): the model is swapped between two checkpoint files."""

try:
    from chimerax.core.undo import UndoAction
except Exception:                                   # keeps the module importable outside ChimeraX
    class UndoAction:
        def __init__(self, name, can_redo=True):
            self.name, self.can_redo = name, can_redo


class ModelSwapUndo(UndoAction):
    """Undo restores the model saved *before* a step, redo the one saved *after* it.

    *swap(path, which)* is called with the checkpoint file and 'before' or 'after'; it must be idempotent."""

    def __init__(self, name, before_file, after_file, swap):
        super().__init__(name)
        self.before_file, self.after_file, self.swap = before_file, after_file, swap

    def undo(self):
        self.swap(self.before_file, "before")

    def redo(self):
        self.swap(self.after_file, "after")


def register(session, action):
    mgr = getattr(session, "undo", None)
    if mgr is not None and hasattr(mgr, "register"):
        try:
            mgr.register(action)
            return True
        except Exception:
            pass
    return False
