"""
DEPRECATED: Legacy factory functions.

This module is kept for backwards compatibility but is now just a wrapper
around unified_factory. New code should use unified_factory directly.

ALL INTERACTIONS WITH EXTERNAL MODELS GO THROUGH THE SIMULATOR INTERFACE.
"""

# Import and re-export unified factory functions for backwards compatibility
from .unified_factory import (
    from_cobra_model,
    from_reframed_model,
    from_simulator,
    load_cobra_model,
    load_reframed_model,
)

# Convenience aliases for backwards compatibility
load_model = {"cobra": load_cobra_model, "reframed": load_reframed_model}

from_model = {"cobra": from_cobra_model, "reframed": from_reframed_model}
