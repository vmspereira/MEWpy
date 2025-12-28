from typing import TYPE_CHECKING, Dict, Iterable, List, Union

import pandas as pd

if TYPE_CHECKING:
    from mewpy.solvers.solution import Solution


class MultiSolution:
    """
    A MultiSolution object is a collection of Solution objects.
    It can be used to compare different simulation methods or to compare the same method with different parameters.
    A MultiSolution object can be created by passing a list of Solution objects to the constructor.
    This object can be exported into a pandas DataFrame or Summary-like object
    using the to_frame(), to_summary() methods, respectively.
    """

    def __init__(self, *solutions: "Solution"):
        """
        A MultiSolution object is a collection of Solution objects.
        It can be used to compare different simulation methods or to compare the same method with different parameters.
        :param solutions: a list of Solution objects
        """
        if not solutions:
            _solutions = {}
        else:
            _solutions = {}

            for solution in solutions:
                setattr(self, f"{solution.method}", solution)
                _solutions[solution.method] = solution

        self._solutions = _solutions

    @property
    def solutions(self) -> Dict[str, "Solution"]:
        """
        Returns a dict of Solution objects by the method name
        :return: a dict of Solution objects by the method name
        """
        return self._solutions

    def to_frame(self) -> pd.DataFrame:
        """
        Returns a pandas DataFrame with the results of the MultiSolution object
        :return: pandas DataFrame by the method name
        """
        frames = []
        columns = []

        for method, solution in self._solutions.items():
            frames.append(solution.to_frame().frame)
            columns.append(method)

        df = pd.concat(frames, axis=1, join="outer", keys=columns)

        return df

    def to_summary(self) -> pd.DataFrame:
        """
        Returns a pandas DataFrame with the summary of the MultiSolution object
        :return: pandas DataFrame by the method name
        """
        frames = []
        columns = []

        for method, solution in self._solutions.items():
            frames.append(solution.to_summary().frame)
            columns.append(method)

        df = pd.concat(frames, axis=1, join="outer", keys=columns)

        return df

    def __repr__(self):
        """Rich representation showing all methods."""
        lines = []
        lines.append("=" * 60)
        lines.append("MultiSolution")
        lines.append("=" * 60)
        lines.append(f"Number of solutions: {len(self._solutions)}")

        if self._solutions:
            lines.append("\nMethods:")
            for method in self._solutions:
                lines.append(f"  - {method}")
            lines.append("\nUse .solutions to access individual solution objects")
            lines.append("Use .to_frame() to get results as DataFrame")

        lines.append("=" * 60)
        return "\n".join(lines)

    def __str__(self):
        return "MultiSolution:" ",".join([method for method in self._solutions])


# TODO: methods stubs and type hinting
class DynamicSolution:
    """
    A DynamicSolution object is a collection of Solution objects.
    It is similar to the MultiSolution object, but it is used to store the results of a dynamic simulation using the
    time point rather than the method name.
    """

    def __init__(self, *solutions: "Solution", time: Iterable = None):
        """
        A DynamicSolution object is a collection of Solution objects.
        It is similar to the MultiSolution object, but it is used to store the results of a dynamic simulation using the
        time point rather than the method name.
        :param solutions: a list of Solution objects
        :param time: a linear space of time points
        """
        if not solutions:
            _solutions = {}
            time = []

        else:
            _solutions = {}

            if time is None:
                time = [i for i in range(len(solutions))]

            else:
                time = list(time)

            for t, solution in zip(time, solutions):
                setattr(self, f"t_{t}", solution)
                _solutions[f"t_{t}"] = solution

        self._solutions = _solutions
        self._time = time

    @property
    def solutions(self) -> Dict[str, "Solution"]:
        """
        Returns a dict of Solution objects by the time point
        :return: a dict of Solution objects by the time point
        """
        return self._solutions

    def to_frame(self) -> pd.DataFrame:
        """
        Returns a pandas DataFrame with the results of the DynamicSolution object
        :return: pandas DataFrame by the time point
        """
        frames = []
        columns = []

        for time, solution in self._solutions.items():
            frames.append(solution.to_frame().frame)
            columns.append(time)

        df = pd.concat(frames, axis=1, join="outer", keys=columns)

        return df

    def to_summary(self) -> pd.DataFrame:
        """
        Returns a pandas DataFrame with the summary of the DynamicSolution object
        :return: pandas DataFrame by the time point
        """
        frames = []
        columns = []

        for time, solution in self._solutions.items():
            frames.append(solution.to_summary().frame)
            columns.append(time)

        df = pd.concat(frames, axis=1, join="outer", keys=columns)

        return df

    def __repr__(self):
        """Rich representation showing time points."""
        lines = []
        lines.append("=" * 60)
        lines.append("DynamicSolution")
        lines.append("=" * 60)
        lines.append(f"Number of time points: {len(self._time)}")

        if self._time:
            lines.append(f"Time range: {min(self._time)} to {max(self._time)}")
            lines.append(f"\nTime points ({len(self._time)} total):")
            # Show first few and last few if many time points
            if len(self._time) <= 10:
                for t in self._time:
                    lines.append(f"  - t_{t}")
            else:
                for t in self._time[:5]:
                    lines.append(f"  - t_{t}")
                lines.append(f"  ... ({len(self._time) - 10} more time points)")
                for t in self._time[-5:]:
                    lines.append(f"  - t_{t}")

            lines.append("\nUse .solutions to access individual time point solutions")
            lines.append("Use .to_frame() to get results as DataFrame")

        lines.append("=" * 60)
        return "\n".join(lines)

    def __str__(self):
        return "DynamicSolution:" ",".join([time for time in self._solutions])


# TODO: methods stubs and type hinting
class KOSolution:
    """
    A KOSolution object is a collection of Solution objects.
    It is similar to the MultiSolution object, but it is used to store the results of a KO simulations.
    """

    def __init__(self, solutions: Union[List["Solution"], Dict[str, "Solution"]], kos: List[str] = None):
        """
        A KOSolution object is a collection of Solution objects.
        It is similar to the MultiSolution object, but it is used to store the results of a KO simulations.
        :param solutions: a list of Solution objects or a dict of Solution objects by the KO name
        :param kos: a list of KO names
        """
        if not solutions:
            _solutions = {}
            kos = []

        else:
            _solutions = {}

            if not kos:

                if isinstance(solutions, dict):

                    kos = list(solutions.keys())
                    solutions = list(solutions.values())

                else:

                    kos = [i for i, _ in enumerate(solutions)]

            else:

                if isinstance(solutions, dict):
                    solutions = list(solutions.values())

            for ko, solution in zip(kos, solutions):
                setattr(self, f"ko_{ko}", solution)
                _solutions[f"ko_{ko}"] = solution

        self._solutions = _solutions
        self._time = kos

    @property
    def solutions(self) -> Dict[str, "Solution"]:
        """
        Returns a dict of Solution objects by the KO name
        :return: a dict of Solution objects by the KO name
        """
        return self._solutions

    def to_frame(self) -> pd.DataFrame:
        """
        Returns a pandas DataFrame with the results of the KOSolution object
        :return: pandas DataFrame by the KO name
        """
        frames = []
        columns = []

        for ko, solution in self._solutions.items():
            frames.append(solution.to_frame().frame)
            columns.append(ko)

        df = pd.concat(frames, axis=1, join="outer", keys=columns)

        return df

    def to_summary(self) -> pd.DataFrame:
        """
        Returns a pandas DataFrame with the summary of the KOSolution object
        :return: pandas DataFrame by the KO name
        """
        frames = []
        columns = []

        for ko, solution in self._solutions.items():
            frames.append(solution.to_summary().frame)
            columns.append(ko)

        df = pd.concat(frames, axis=1, join="outer", keys=columns)

        return df

    def __repr__(self):
        """Rich representation showing knocked out targets."""
        lines = []
        lines.append("=" * 60)
        lines.append("KOSolution")
        lines.append("=" * 60)
        lines.append(f"Number of knockouts: {len(self._time)}")

        if self._time:
            lines.append(f"\nKnockouts ({len(self._time)} total):")
            # Show first few and last few if many KOs
            if len(self._time) <= 10:
                for ko in self._time:
                    lines.append(f"  - {ko}")
            else:
                for ko in self._time[:5]:
                    lines.append(f"  - {ko}")
                lines.append(f"  ... ({len(self._time) - 10} more knockouts)")
                for ko in self._time[-5:]:
                    lines.append(f"  - {ko}")

            lines.append("\nUse .solutions to access individual KO solutions")
            lines.append("Use .to_frame() to get results as DataFrame")

        lines.append("=" * 60)
        return "\n".join(lines)

    def __str__(self):
        return "KOSolution:" ",".join([time for time in self._solutions])
