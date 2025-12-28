# Copyright (C) 2019- Centre of Biological Engineering,
#     University of Minho, Portugal
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""
##############################################################################
Logging configuration for MEWpy

Author: Vitor Pereira
##############################################################################
"""
import logging
import logging.config


def configure_logging(level="INFO"):
    """
    Configure logging for MEWpy.

    Parameters
    ----------
    level : str, optional
        Logging level. One of: DEBUG, INFO, WARNING, ERROR, CRITICAL.
        Default is INFO.
    """
    DEFAULT_LOGGING_CONFIG = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "basic": {
                "format": "[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s",
                "datefmt": "%Y-%m-%d %H:%M:%S",
            }
        },
        "handlers": {
            "console": {
                "formatter": "basic",
                "class": "logging.StreamHandler",
                "stream": "ext://sys.stderr",
                "level": level,
            }
        },
        "loggers": {
            "mewpy": {
                "handlers": ["console"],
                "level": level,
                "propagate": False,
            },
        },
    }

    logging.config.dictConfig(DEFAULT_LOGGING_CONFIG)


def get_logger(module):
    """
    Get a logger for the given module.

    Parameters
    ----------
    module : str
        Module name, typically __name__

    Returns
    -------
    logging.Logger
        Logger instance
    """
    return logging.getLogger(module)
