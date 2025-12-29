from collections.abc import Callable
from functools import partial, wraps
from typing import TYPE_CHECKING, Union

import pandas as pd

if TYPE_CHECKING:
    from mewpy.germ.models import Model
    from mewpy.germ.variables import Variable


class HistoryManager:

    def __init__(self):

        self._history = []
        self._undo_able_commands = []
        self._temp_stack = []
        self._redo_able_commands = []

    def __str__(self):
        return f"History: {len(self._undo_able_commands)} undos and {len(self._redo_able_commands)} redos"

    def __repr__(self):
        """Rich representation showing history state."""
        lines = []
        lines.append("=" * 60)
        lines.append("History Manager")
        lines.append("=" * 60)

        # Undo/redo availability
        try:
            undo_count = len(self._undo_able_commands)
            redo_count = len(self._redo_able_commands)

            lines.append(f"{'Undo available:':<20} {undo_count}")
            lines.append(f"{'Redo available:':<20} {redo_count}")

            # Total history size
            history_size = len(self._history)
            lines.append(f"{'Total history:':<20} {history_size}")
        except:
            pass

        # Current position indicator
        try:
            if undo_count > 0 and redo_count == 0:
                lines.append(f"{'Position:':<20} At end (can undo)")
            elif undo_count == 0 and redo_count > 0:
                lines.append(f"{'Position:':<20} At start (can redo)")
            elif undo_count > 0 and redo_count > 0:
                lines.append(f"{'Position:':<20} Middle (can undo/redo)")
            else:
                lines.append(f"{'Position:':<20} Empty")
        except:
            pass

        # Show recent history entries
        try:
            if len(self._history) > 0:
                lines.append(f"{'Recent actions:':<20}")
                recent = self._history[-3:] if len(self._history) > 3 else self._history
                for entry in recent:
                    method_name = entry[0] if len(entry) > 0 else "unknown"
                    obj_str = entry[3] if len(entry) > 3 else ""
                    # Truncate object string if too long
                    if len(obj_str) > 25:
                        obj_str = obj_str[:22] + "..."
                    if obj_str:
                        lines.append(f"{'  -':<20} {method_name} ({obj_str})")
                    else:
                        lines.append(f"{'  -':<20} {method_name}")
        except:
            pass

        lines.append("=" * 60)
        return "\n".join(lines)

    def _repr_html_(self):
        """Pandas-like HTML representation for Jupyter notebooks."""
        from mewpy.util.html_repr import render_html_table

        rows = []

        # Undo/redo availability
        try:
            undo_count = len(self._undo_able_commands)
            redo_count = len(self._redo_able_commands)

            rows.append(("Undo available", str(undo_count)))
            rows.append(("Redo available", str(redo_count)))

            # Total history size
            history_size = len(self._history)
            rows.append(("Total history", str(history_size)))
        except:
            pass

        # Current position indicator
        try:
            if undo_count > 0 and redo_count == 0:
                rows.append(("Position", "At end (can undo)"))
            elif undo_count == 0 and redo_count > 0:
                rows.append(("Position", "At start (can redo)"))
            elif undo_count > 0 and redo_count > 0:
                rows.append(("Position", "Middle (can undo/redo)"))
            else:
                rows.append(("Position", "Empty"))
        except:
            pass

        # Show recent history entries
        try:
            if len(self._history) > 0:
                rows.append(("Recent actions", ""))
                recent = self._history[-3:] if len(self._history) > 3 else self._history
                for entry in recent:
                    method_name = entry[0] if len(entry) > 0 else "unknown"
                    obj_str = entry[3] if len(entry) > 3 else ""
                    # Truncate object string if too long
                    if len(obj_str) > 25:
                        obj_str = obj_str[:22] + "..."
                    if obj_str:
                        rows.append(("  -", f"{method_name} ({obj_str})"))
                    else:
                        rows.append(("  -", method_name))
        except:
            pass

        return render_html_table("History Manager", rows)

    @property
    def history(self):

        return pd.DataFrame(data=self._history, columns=["method", "args", "kwargs", "object"])

    @property
    def undo_able_commands(self):
        return self._undo_able_commands

    @property
    def redo_able_commands(self):
        return self._redo_able_commands

    def _do(self, undo=True):

        if undo:
            method = self.undo_able_commands.pop()
            redo_command = self._temp_stack.pop()
            self.redo_able_commands.append(redo_command)

        else:
            method = self.redo_able_commands.pop()

        method()

    def undo(self) -> None:
        self._do(undo=True)

    def redo(self) -> None:
        self._do(undo=False)

    def reset(self) -> None:

        while len(self.undo_able_commands) > 0:
            self.undo()

    def restore(self) -> None:

        while len(self.redo_able_commands) > 0:
            self.redo()

    def __call__(self, *args, **kwargs) -> None:

        self.queue_command(*args, **kwargs)

    def queue_command(
        self,
        undo_func: Callable,
        func: Callable,
        undo_args: tuple = None,
        undo_kwargs: dict = None,
        args: tuple = None,
        kwargs: dict = None,
        obj: "Model" = None,
    ) -> None:

        if not undo_args:
            undo_args = ()

        if not undo_kwargs:
            undo_kwargs = {}

        self.undo_able_commands.append(partial(undo_func, *undo_args, **undo_kwargs))

        if not args:
            args = ()

        if not kwargs:
            kwargs = {}

        self._temp_stack.append(partial(func, *args, **kwargs))

        self._history.append((func.__name__, str(args), str(kwargs), str(obj)))


def recorder(func: Callable):

    @wraps(func)
    def wrapper(self: Union["Model", "Variable"], value):

        history = self.history

        old_value = getattr(self, func.__name__)

        if old_value != value:

            history.queue_command(undo_func=func, undo_args=(self, old_value), func=func, args=(self, value))

        return func(self, value)

    return wrapper
