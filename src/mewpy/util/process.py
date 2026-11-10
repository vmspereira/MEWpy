# Copyright (C) 2019-2026 Vitor Pereira and Centre of Biological Engineering,
#     University of Minho, Portugal

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
##############################################################################
Multiprocessing module.

Authors: Vitor Pereira
##############################################################################
"""
import copy
import importlib.util
import logging
from abc import ABC, abstractmethod

from .constants import EAConstants, ModelConstants

logger = logging.getLogger(__name__)

MP_Evaluators = []

import multiprocessing
from multiprocessing import Process
from multiprocessing.pool import Pool as MPPool

# MEWpy pools use a local 'fork' context: 'spawn' (the default on macOS since Python 3.8) requires
# pickling all objects, which fails with CPLEX/REFRAMED, while 'fork' copies the process memory.
# The user's global start method is left untouched.
_FORK_CONTEXT = None
if "fork" in multiprocessing.get_all_start_methods():
    _FORK_CONTEXT = multiprocessing.get_context("fork")

MP_Evaluators.append("nodaemon")

# pathos
try:
    import pathos.multiprocessing as multiprocessing
    from pathos.multiprocessing import Pool

    MP_Evaluators.append("mp")

except ImportError:
    import multiprocessing
    from multiprocessing.pool import Pool

    MP_Evaluators.append("mp")

# dask
try:
    import dask

    MP_Evaluators.append("dask")
except ImportError:
    pass
# pyspark
try:
    from pyspark import SparkConf, SparkContext

    MP_Evaluators.append("spark")
except ImportError:
    pass


class NoDaemonProcess(Process):

    def _get_daemon(self):
        return False

    def _set_daemon(self, value):
        pass

    daemon = property(_get_daemon, _set_daemon)


if _FORK_CONTEXT is not None:

    class NoDaemonForkProcess(_FORK_CONTEXT.Process):
        daemon = property(NoDaemonProcess._get_daemon, NoDaemonProcess._set_daemon)


class NoDaemonProcessPool(MPPool):
    """A pool whose workers are not daemonic, so that they can start their own processes."""

    @staticmethod
    def Process(ctx, *args, **kwds):
        # Pool creates its workers with Process(ctx, ...) since Python 3.8
        if _FORK_CONTEXT is not None and ctx is _FORK_CONTEXT:
            return NoDaemonForkProcess(*args, **kwds)
        return NoDaemonProcess(*args, **kwds)


def _new_pool(n_processes, pool_class=None):
    """A process pool using the local 'fork' context when the standard library pool is used."""
    if pool_class is None:
        pool_class = Pool
    if isinstance(pool_class, type) and issubclass(pool_class, MPPool):
        return pool_class(n_processes, context=_FORK_CONTEXT)
    return pool_class(n_processes)


def cpu_count():
    """The number of cpus
    Return EAConstants.NUM_CPUS if it is set to a positive number
    otherwise return half of the available cpus.
    """
    if EAConstants.NUM_CPUS > 0:
        return EAConstants.NUM_CPUS
    else:
        try:
            return max(1, multiprocessing.cpu_count() // 2)
        except (ImportError, NotImplementedError):
            return 1


class Evaluable(ABC):

    @abstractmethod
    def evaluator(self, candidates, *args):
        raise NotImplementedError


class Evaluator(ABC):
    """An interface for multiprocessing evaluators

    Raises:
        NotImplementedError: Requires an evaluated method to
        be implemented.
    """

    @abstractmethod
    def evaluate(self, candidates, args):
        raise NotImplementedError


class MultiProcessorEvaluator(Evaluator):

    def __init__(self, evaluator, mp_num_cpus):
        """A multiprocessing evaluator

        Args:
            evaluator(function): Evaluation function.
            mp_num_cpus(int): Number of CPUs

        When using COBRApy, memory resources are not released after each
        pool map. As such, the pool needs to be instantiated and closed at
        each iteration.
        """
        self.mp_num_cpus = mp_num_cpus
        self.evaluator = evaluator
        self.__name__ = self.__class__.__name__

    def evaluate(self, candidates, args):
        """
        Values in args will be ignored and not passed to the evaluator to avoid unnecessary pickling in inspyred.
        """
        pool = _new_pool(self.mp_num_cpus)
        try:
            return pool.map(self.evaluator, candidates)
        finally:
            # close, not terminate: terminate sends SIGTERM to workers, which may carry signal
            # handlers inherited from the parent (e.g. Ray's)
            pool.close()

    def __call__(self, candidates, args):
        return self.evaluate(candidates, args)


class NoDaemonMultiProcessorEvaluator(Evaluator):

    def __init__(self, evaluator, mp_num_cpus):
        """A multiprocessing evaluator

        Args:
            evaluator(function): Evaluation function.
            mp_num_cpus(int): Number of CPUs
        """
        self.mp_num_cpus = mp_num_cpus
        self.evaluator = evaluator
        self.__name__ = self.__class__.__name__
        logger.debug("nodaemon")

    def evaluate(self, candidates, args):
        """
        Values in args will be ignored and not passed to the evaluator
        to avoid unnecessary pickling in inspyred.
        """
        pool = _new_pool(self.mp_num_cpus, NoDaemonProcessPool)
        try:
            return pool.map(self.evaluator, candidates)
        finally:
            pool.close()

    def __call__(self, candidates, args):
        return self.evaluate(candidates, args)


class DaskEvaluator(Evaluator):

    def __init__(self, evaluator, mp_num_cpus, scheduler="processes"):
        """A Dask multiprocessing evaluator

        Args:
            evaluator (function): Evaluation function.
            mp_num_cpus (int): Number of CPUs.
        """
        self.evaluator = evaluator
        self.scheduler = scheduler
        self.__name__ = self.__class__.__name__

    def evaluate(self, candidates, args):
        with dask.config.set(scheduler=self.scheduler):
            return list(dask.compute(*[dask.delayed(self.evaluator)(c) for c in candidates]))

    def __call__(self, candidates, args):
        return self.evaluate(candidates, args)


class SparkEvaluator(Evaluator):
    def __init__(self, evaluator, mp_num_cpus):
        """A Spark multiprocessing evaluator

        Args:
            evaluator (function): Evaluation function.
            mp_num_cpus (int): Number of CPUs.
        """
        self.evaluator = evaluator
        self.spark_conf = SparkConf().setAppName("mewpy").setMaster(f"local[{mp_num_cpus}]")
        self.spark_context = SparkContext(conf=self.spark_conf)
        self.__name__ = self.__class__.__name__

    def evaluate(self, candidates, args):
        solutions_to_evaluate = self.spark_context.parallelize(candidates)
        return solutions_to_evaluate.map(lambda s: self.evaluator(s)).collect()

    def __call__(self, candidates, args):
        return self.evaluate(candidates, args)


# ray: imported when a Ray evaluator is created, not with mewpy (it is slow and writes to stderr)
if importlib.util.find_spec("ray") is not None:
    MP_Evaluators.append("ray")


class RayActor:
    """
    Each actor (worker) has a solver instance to overcome the need
    to serialize solvers which may not be pickable.
    The solver is not reset before each evaluation.
    """

    def __init__(self, problem):
        self.problem = copy.deepcopy(problem)

    def evaluate_candidates(self, candidates):
        """Evaluates a sublist of candidates"""
        if getattr(self.problem, "evaluator", False):
            return self.problem.evaluator(candidates, None)
        else:
            res = []
            for candidate in candidates:
                res.append(self.problem.evaluate(candidate))
            return res


class RayActorF:
    """
    Each actor (worker) has a solver instance to overcome the need to
    serialize solvers which may not be pickable.
    The solver is not reset before each evaluation.
    """

    def __init__(self, func):
        self.func = copy.deepcopy(func)

    def evaluate_candidates(self, candidates):
        """Evaluates a sublist of candidates"""
        res = []
        for candidate in candidates:
            res.append(self.func(candidate))
        return res


_RAY_ACTOR_CLASSES = {}


def _ray_actor_class(cls):
    """The Ray actor class of RayActor/RayActorF, made on first use."""
    import ray

    if cls not in _RAY_ACTOR_CLASSES:
        _RAY_ACTOR_CLASSES[cls] = ray.remote(cls)
    return _RAY_ACTOR_CLASSES[cls]


class RayEvaluator(Evaluator):
    def __init__(self, problem, number_of_actors, isfunc=False):
        """A ray actor responsible for performing evaluations.

        Args:
            problem: A class implementing an evaluator(list_of_candidates,**kwargs)
            number_of_actors (int): Number of workers
        """
        import ray

        ray.init(ignore_reinit_error=True)
        actor_class = _ray_actor_class(RayActorF if isfunc else RayActor)
        self.actors = [actor_class.remote(problem) for _ in range(number_of_actors)]
        self.number_of_actors = len(self.actors)
        self.__name__ = self.__class__.__name__
        logger.info("Using %s workers.", self.number_of_actors)

    def evaluate(self, candidates, args):
        """
        Divides the candidates into sublists to be evaluated by each actor
        """
        import ray

        size = len(candidates) // self.number_of_actors
        n = len(candidates) % self.number_of_actors
        sub_lists = []
        p = 0
        for _ in range(self.number_of_actors):
            d = size + 1 if n > 0 else size
            sub_lists.append(candidates[p : p + d])
            p = p + d
            n -= 1
        values = []
        for i in range(self.number_of_actors):
            actor = self.actors[i]
            values.append(actor.evaluate_candidates.remote(sub_lists[i]))
        r = ray.get(values)
        result = []
        for x in r:
            result += x
        return result

    def __call__(self, candidates, args):
        return self.evaluate(candidates, args)


def get_mp_evaluators():
    """ "Returns the list of available multiprocessing evaluators."""
    return MP_Evaluators


def get_evaluator(problem, n_mp=None, evaluator=None):
    """Returns a multiprocessing evaluator

    Args:
        problem: a class implementing an evaluate(candidate) function
        n_mp (int, optional): The number of cpus. Defaults to cpu_count().
        evaluator (str, optional): The evaluator name: options 'ray','dask','spark'.
            Defaults to ModelConstants.MP_EVALUATOR, read when the function is called.

    Returns:
        Evaluator: A multiprocessing evaluator instance based on the specified evaluator type
    """
    n_mp = cpu_count() if n_mp is None else n_mp
    evaluator = ModelConstants.MP_EVALUATOR if evaluator is None else evaluator
    if evaluator == "ray" and "ray" in MP_Evaluators:
        return RayEvaluator(problem, n_mp)
    elif evaluator == "nodaemon" and "nodaemon" in MP_Evaluators:
        return NoDaemonMultiProcessorEvaluator(problem.evaluate, n_mp)
    elif evaluator == "dask" and "dask" in MP_Evaluators:
        return DaskEvaluator(problem.evaluate, n_mp)
    elif evaluator == "spark" and "spark" in MP_Evaluators:
        return SparkEvaluator(problem.evaluate, n_mp)
    else:
        return MultiProcessorEvaluator(problem.evaluate, n_mp)


def get_fevaluator(func, n_mp=None, evaluator=None):
    """Returns a multiprocessing evaluator

    Args:
        func: an evaluation function
        n_mp (int, optional): The number of cpus. Defaults to cpu_count().
        evaluator (str, optional): The evaluator name: options 'ray','dask','spark'.
            Defaults to ModelConstants.MP_EVALUATOR, read when the function is called.

    Returns:
        Evaluator: A multiprocessing evaluator instance based on the specified evaluator type
    """
    n_mp = cpu_count() if n_mp is None else n_mp
    evaluator = ModelConstants.MP_EVALUATOR if evaluator is None else evaluator
    if evaluator == "ray" and "ray" in MP_Evaluators:
        return RayEvaluator(func, n_mp, isfunc=True)
    elif evaluator == "nodaemon" and "nodaemon" in MP_Evaluators:
        return NoDaemonMultiProcessorEvaluator(func, n_mp)
    elif evaluator == "dask" and "dask" in MP_Evaluators:
        return DaskEvaluator(func, n_mp)
    elif evaluator == "spark" and "spark" in MP_Evaluators:
        return SparkEvaluator(func, n_mp)
    else:
        return MultiProcessorEvaluator(func, n_mp)
