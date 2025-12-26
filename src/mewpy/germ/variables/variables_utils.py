from typing import TYPE_CHECKING, Tuple, Union

if TYPE_CHECKING:
    from mewpy.germ.variables import Variable


def initialize_coefficients(coefficients: Union[Tuple[float, ...], None]) -> Tuple[float, ...]:
    """
    Initialize coefficients for Gene, Target, or Regulator variables.

    Helper function to standardize coefficient initialization across variable types.
    If no coefficients provided, defaults to (0.0, 1.0) for inactive/active states.

    :param coefficients: Optional tuple of coefficients, or None for defaults
    :return: Tuple of coefficients (minimum 0.0, 1.0 if None)
    """
    if not coefficients:
        return (0.0, 1.0)
    else:
        return tuple(coefficients)


def coefficients_setter(instance: "Variable", value: Union[Tuple[float, float], Tuple[float], float]):
    """
    Setter for the coefficients attribute of a variable.
    :param instance: the variable instance
    :param value: the value to be set
    :return:
    """
    if value is None:
        value = (0, 1)

    elif isinstance(value, (int, float)):
        value = (value, value)

    elif len(value) == 1:
        value = (value[0], value[0])

    elif len(value) > 1:
        value = tuple(value)

    else:
        raise ValueError("Invalid value for coefficients")

    # if it is a reaction, bounds must be set
    if hasattr(instance, "_bounds"):
        instance._bounds = value

    # if it is a metabolite, the bounds coefficient of the exchange reaction must be returned
    elif hasattr(instance, "exchange_reaction"):
        if hasattr(instance.exchange_reaction, "_bounds"):
            instance.exchange_reaction._bounds = value

    else:
        instance._coefficients = value

    if instance.model:
        instance.model.notify()
