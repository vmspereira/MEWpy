"""
Unified factory for creating MetabolicModel instances from external simulators.

This factory provides a consistent interface for creating GERM-compatible models
from external simulators (COBRApy, reframed) using either:
1. Original MetabolicModel with simulator backends
2. SimulatorBasedMetabolicModel (pure simulator wrapper)

The factory automatically selects the best approach and provides a unified interface.
"""
from typing import Union, Any, TYPE_CHECKING

if TYPE_CHECKING:
    from mewpy.simulation.simulation import Simulator
    from .metabolic import MetabolicModel
    from .simulator_model import SimulatorBasedMetabolicModel

def create_model_from_simulator(simulator: 'Simulator', 
                               approach: str = 'wrapper',
                               identifier: Any = None,
                               **kwargs) -> Union['MetabolicModel', 'SimulatorBasedMetabolicModel']:
    """
    Create a GERM-compatible MetabolicModel from a simulator.
    
    :param simulator: External simulator instance (COBRApy, reframed, etc.)
    :param approach: 'wrapper' (SimulatorBasedMetabolicModel) or 'backend' (MetabolicModel with backend)
    :param identifier: Model identifier (defaults to simulator's model ID)
    :param kwargs: Additional arguments for model creation
    :return: GERM-compatible MetabolicModel
    """
    if identifier is None:
        identifier = getattr(simulator, 'id', 'simulator_model')
    
    if approach == 'wrapper':
        # Use SimulatorBasedMetabolicModel (pure simulator wrapper)
        from .simulator_model import SimulatorBasedMetabolicModel
        return SimulatorBasedMetabolicModel(simulator, identifier=identifier, **kwargs)
    
    elif approach == 'backend':
        # Use original MetabolicModel with simulator backend
        # This would require converting simulator data to GERM variables first
        # and then setting up the backend - more complex but preserves full GERM functionality
        raise NotImplementedError("Backend approach not yet implemented")
    
    else:
        raise ValueError(f"Unknown approach: {approach}. Use 'wrapper' or 'backend'")


def load_cobra_model(model_path: str, 
                    approach: str = 'wrapper',
                    identifier: Any = None,
                    **sim_kwargs) -> Union['MetabolicModel', 'SimulatorBasedMetabolicModel']:
    """
    Load a COBRApy model and create a GERM-compatible MetabolicModel.
    
    :param model_path: Path to the model file (SBML, JSON, etc.)
    :param approach: 'wrapper' or 'backend' approach
    :param identifier: Model identifier
    :param sim_kwargs: Additional arguments for simulator creation
    :return: GERM-compatible MetabolicModel
    """
    try:
        import cobra
        from mewpy.simulation.cobra import Simulation
        
                # Load COBRApy model
        cobra_model = cobra.io.read_sbml_model(model_path)
        
        # Create simulator
        simulator = Simulation(cobra_model, **sim_kwargs)
        
        # Create GERM model using specified approach
        return create_model_from_simulator(
            simulator, 
            approach=approach, 
            identifier=identifier or cobra_model.id
        )
        
    except ImportError as e:
        raise ImportError(f"COBRApy is required to load cobra models: {e}")


def load_reframed_model(model_path: str,
                       approach: str = 'wrapper', 
                       identifier: Any = None,
                       **sim_kwargs) -> Union['MetabolicModel', 'SimulatorBasedMetabolicModel']:
    """
    Load a reframed model and create a GERM-compatible MetabolicModel.
    
    :param model_path: Path to the model file (SBML, JSON, etc.)
    :param approach: 'wrapper' or 'backend' approach
    :param identifier: Model identifier
    :param sim_kwargs: Additional arguments for simulator creation
    :return: GERM-compatible MetabolicModel
    """
    try:
        from reframed.io.sbml import load_cbmodel
        from mewpy.simulation.reframed import Simulation
        
        # Load reframed model
        reframed_model = load_cbmodel(model_path)
        
        # Create simulator
        simulator = Simulation(reframed_model, **sim_kwargs)
        
        # Create GERM model using specified approach
        return create_model_from_simulator(
            simulator,
            approach=approach,
            identifier=identifier or reframed_model.id
        )
        
    except ImportError as e:
        raise ImportError(f"reframed is required to load reframed models: {e}")


def from_cobra_model(cobra_model, 
                    approach: str = 'wrapper',
                    identifier: Any = None,
                    **sim_kwargs) -> Union['MetabolicModel', 'SimulatorBasedMetabolicModel']:
    """
    Create a GERM-compatible MetabolicModel from a COBRApy model object.
    
    :param cobra_model: COBRApy Model instance
    :param approach: 'wrapper' or 'backend' approach
    :param identifier: Model identifier
    :param sim_kwargs: Additional arguments for simulator creation
    :return: GERM-compatible MetabolicModel
    """
    try:
        from mewpy.simulation.cobra import Simulation
        
        # Create simulator
        simulator = Simulation(cobra_model, **sim_kwargs)
        
        # Create GERM model using specified approach
        return create_model_from_simulator(
            simulator,
            approach=approach, 
            identifier=identifier or cobra_model.id
        )
        
    except ImportError as e:
        raise ImportError(f"COBRApy simulation support is required: {e}")


def from_reframed_model(reframed_model,
                       approach: str = 'wrapper',
                       identifier: Any = None, 
                       **sim_kwargs) -> Union['MetabolicModel', 'SimulatorBasedMetabolicModel']:
    """
    Create a GERM-compatible MetabolicModel from a reframed CBModel object.
    
    :param reframed_model: reframed CBModel instance
    :param approach: 'wrapper' or 'backend' approach
    :param identifier: Model identifier
    :param sim_kwargs: Additional arguments for simulator creation
    :return: GERM-compatible MetabolicModel
    """
    try:
        from mewpy.simulation.reframed import Simulation
        
        # Create simulator
        simulator = Simulation(reframed_model, **sim_kwargs)
        
        # Create GERM model using specified approach
        return create_model_from_simulator(
            simulator,
            approach=approach,
            identifier=identifier or reframed_model.id
        )
        
    except ImportError as e:
        raise ImportError(f"reframed simulation support is required: {e}")


def from_simulator(simulator: 'Simulator',
                  approach: str = 'wrapper',
                  identifier: Any = None) -> Union['MetabolicModel', 'SimulatorBasedMetabolicModel']:
    """
    Create a GERM-compatible MetabolicModel from any simulator instance.
    
    :param simulator: Simulator instance (COBRApy, reframed, etc.)
    :param approach: 'wrapper' or 'backend' approach
    :param identifier: Optional model identifier
    :return: GERM-compatible MetabolicModel
    """
    return create_model_from_simulator(simulator, approach=approach, identifier=identifier)


# Convenience aliases for different loading methods
load_external_model = {
    'cobra': load_cobra_model,
    'reframed': load_reframed_model
}

from_external_model = {
    'cobra': from_cobra_model,
    'reframed': from_reframed_model
}

# Default functions (using wrapper approach)
def load_model(model_path: str, backend: str = 'cobra', **kwargs):
    """
    Load an external model using the specified backend.
    
    :param model_path: Path to model file
    :param backend: 'cobra' or 'reframed'
    :param kwargs: Additional arguments
    :return: GERM-compatible MetabolicModel
    """
    if backend not in load_external_model:
        raise ValueError(f"Unknown backend: {backend}. Use 'cobra' or 'reframed'")
    
    return load_external_model[backend](model_path, **kwargs)


def from_model(external_model, backend: str = None, **kwargs):
    """
    Create GERM model from external model object.
    
    :param external_model: External model object
    :param backend: 'cobra' or 'reframed' (auto-detected if None)
    :param kwargs: Additional arguments
    :return: GERM-compatible MetabolicModel
    """
    if backend is None:
        # Auto-detect backend
        model_type = type(external_model).__name__
        if 'cobra' in model_type.lower() or hasattr(external_model, 'reactions'):
            try:
                import cobra
                if isinstance(external_model, cobra.Model):
                    backend = 'cobra'
            except ImportError:
                pass
        
        if backend is None:
            try:
                from reframed.core.cbmodel import CBModel
                if isinstance(external_model, CBModel):
                    backend = 'reframed'
            except ImportError:
                pass
        
        if backend is None:
            raise ValueError("Could not auto-detect backend. Please specify 'cobra' or 'reframed'")
    
    if backend not in from_external_model:
        raise ValueError(f"Unknown backend: {backend}. Use 'cobra' or 'reframed'")
    
    return from_external_model[backend](external_model, **kwargs)


# Main unified factory function
def unified_factory(source, **kwargs):
    """
    Unified factory for creating GERM-compatible models from various sources.
    
    :param source: Can be:
        - External model object (COBRApy Model, reframed CBModel)
        - Simulator instance
        - File path to model
        - Model identifier string
    :param kwargs: Additional arguments
    :return: GERM-compatible MetabolicModel
    """
    # Handle external model objects
    try:
        import cobra
        if isinstance(source, cobra.Model):
            return from_cobra_model(source, **kwargs)
    except ImportError:
        pass
    
    try:
        from reframed.core.cbmodel import CBModel
        if isinstance(source, CBModel):
            return from_reframed_model(source, **kwargs)
    except ImportError:
        pass
    
    # Handle simulator objects
    try:
        from mewpy.simulation.simulation import Simulator
        if isinstance(source, Simulator):
            return from_simulator(source, **kwargs)
    except ImportError:
        pass
    
    # Handle string inputs (file paths or identifiers)
    if isinstance(source, str):
        import os
        if os.path.exists(source):
            # It's a file path - auto-detect backend and load
            if source.endswith('.xml') or source.endswith('.sbml'):
                return load_cobra_model(source, **kwargs)
            else:
                # Default to cobra for other formats
                return load_cobra_model(source, **kwargs)
        else:
            # It's an identifier - create an empty model
            # For backwards compatibility with IO engines, just issue deprecation warning
            import warnings
            warnings.warn(
                "Creating MetabolicModel from identifier is deprecated. "
                "Use unified_factory with external model objects instead.",
                DeprecationWarning,
                stacklevel=2
            )
            from .metabolic import MetabolicModel
            return MetabolicModel(identifier=source)
    
    raise TypeError(f"Cannot create model from {type(source)}. "
                   f"Expected external model object, simulator, file path, or identifier string.")
