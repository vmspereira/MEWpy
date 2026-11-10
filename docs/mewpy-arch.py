# Copyright (C) 2026 Vitor Pereira

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
"""Draws the architecture figure of the documentation, docs/mewpy-arch.png.

Run from the repository root: python docs/mewpy-arch.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

EDGE = "#44546a"
PROBLEM_BG, PROBLEM_BOX = "#f2f2f2", "#dbe5f1"
METHOD_EDGE = "#f0b400"
SIMULATION_BG = "#ebf1de"
EA_BG, EA_BOX = "#f2dcdb", "#f9ecec"
BASE_BG = "#e7e6e6"

METHODS = [
    ("Stoichiometric\nConstraints", ["Reactions", "Genes", "Medium", "Cofactor swap"]),
    ("Enzymatic\nConstraints", ["GECKO", "sMOMENT", "ETFL"]),
    ("Regulatory\nConstraints", ["OptRAM", "OptORF"]),
    ("Kinetic\nModels", ["Vmax"]),
    ("Microbial\nCommunities", ["Organisms"]),
]
SIMULATION = [
    "FBA",
    "pFBA",
    "MOMA",
    "lMOMA",
    "ROOM",
    "FVA",
    "RFBA",
    "SRFBA",
    "PROM",
    "CoRegFlux",
    "dynamic FBA",
    "dynamic RFBA",
    "SteadyCom",
    "kinetic (ODE)",
]


def box(ax, x, y, w, h, color, edge=EDGE, lw=1.6, radius=0.6):
    ax.add_patch(
        FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0,rounding_size=%s" % radius, facecolor=color, edgecolor=edge, lw=lw
        )
    )


def arrow(ax, start, end, both=False):
    style = "<|-|>" if both else "-|>"
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle=style, mutation_scale=16, lw=2.2, color="black"))


def bullets(ax, x, y, items, size=10.5, step=2.4):
    for i, item in enumerate(items):
        ax.text(x, y - i * step, "• " + item, fontsize=size, va="top")


def draw(path):
    fig = plt.figure(figsize=(11, 9.2), dpi=150)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 84)
    ax.axis("off")

    # the optimization problem
    box(ax, 1, 49, 98, 34, PROBLEM_BG, lw=2)
    ax.text(50, 51, "Computational Strain Optimization Problem", ha="center", fontsize=13, weight="bold")
    box(ax, 2.5, 75, 21, 6.5, PROBLEM_BOX)
    ax.text(3.5, 78.2, "Targeted Product", fontsize=12, weight="bold", va="center")
    box(ax, 2.5, 63, 21, 10.5, PROBLEM_BOX)
    ax.text(3.5, 71.8, "Strategy:", fontsize=12, weight="bold", va="center")
    bullets(ax, 3.8, 69.3, ["Over/Under expression", "Knockout"])
    box(ax, 2.5, 54.5, 21, 7, PROBLEM_BOX)
    ax.text(3.5, 58, "Environmental\nconditions", fontsize=12, weight="bold", va="center")

    box(ax, 25.5, 54.5, 72.5, 27, PROBLEM_BG)
    ax.text(26.8, 79.6, "Method:", fontsize=12, weight="bold", va="center")
    width, gap, left = 13.2, 0.9, 26.8
    for i, (title, items) in enumerate(METHODS):
        x = left + i * (width + gap)
        box(ax, x, 56, width, 21.5, "white", edge=METHOD_EDGE, lw=1.8, radius=0.3)
        ax.text(x + width / 2, 74.4, title, ha="center", va="center", fontsize=10.5, weight="bold")
        bullets(ax, x + 1, 70.3, items, size=10)

    # solving the problem
    arrow(ax, (71, 48.5), (71, 44.5))
    ax.text(69.5, 46.5, "SOLVE", ha="right", va="center", fontsize=11, weight="bold")
    arrow(ax, (75, 44.5), (75, 48.5))
    ax.text(76.5, 46.5, "SOLUTIONS", ha="left", va="center", fontsize=11, weight="bold")

    # phenotype simulation
    box(ax, 1, 7.5, 21, 36, SIMULATION_BG)
    ax.text(2.5, 40.3, "Phenotype\nSimulation", fontsize=12, weight="bold", va="center")
    half = (len(SIMULATION) + 1) // 2
    bullets(ax, 2.5, 36.3, SIMULATION[:half], size=9, step=3.9)
    bullets(ax, 11.5, 36.3, SIMULATION[half:], size=9, step=3.9)

    # evaluation
    box(ax, 25.5, 33.5, 33, 10, EA_BG)
    ax.text(42, 40.6, "Optimization Objectives:", ha="center", va="center", fontsize=12, weight="bold")
    ax.text(
        42,
        36.5,
        "(BPCY; WYIELD; target flux; batch titer,\nyield and productivity; etc.)",
        ha="center",
        va="center",
        fontsize=9.5,
    )
    box(ax, 30, 21, 24, 6.5, EA_BG)
    ax.text(42, 24.25, "Evaluation", ha="center", va="center", fontsize=12, weight="bold")
    box(ax, 30, 7.5, 24, 8.5, EA_BG)
    ax.text(42, 11.75, "Model\nConstraints", ha="center", va="center", fontsize=12, weight="bold")
    arrow(ax, (42, 27.8), (42, 33.2), both=True)

    # the evolutionary algorithm
    box(ax, 61, 7.5, 38, 36, EA_BG)
    ax.text(
        80,
        40.3,
        "Single- or Multi-Objective\nEvolutionary Algorithm",
        ha="center",
        va="center",
        fontsize=12,
        weight="bold",
    )
    ax.text(
        80, 34.6, "(GA; SA; NSGA-II; SPEA2; NSGA-III)\ninspyred or jMetalPy", ha="center", va="center", fontsize=9.5
    )
    box(ax, 64, 22, 32, 6.5, EA_BOX)
    ax.text(80, 25.25, "Evolutionary Process", ha="center", va="center", fontsize=12, weight="bold")
    box(ax, 64, 9.5, 32, 6.5, EA_BOX)
    ax.text(70.5, 12.75, "Decode", ha="center", va="center", fontsize=12, weight="bold")
    ax.text(84, 12.75, "[0.5; 2.0; 0.0; 1.5]", ha="center", va="center", fontsize=10.5, color="#c00000")
    arrow(ax, (80, 16.3), (80, 21.7), both=True)

    # data flow between simulation, evaluation and the algorithm
    arrow(ax, (22.4, 24.25), (29.6, 24.25))
    arrow(ax, (54.4, 24.25), (60.6, 24.25))
    arrow(ax, (60.6, 11.75), (54.4, 11.75))
    arrow(ax, (29.6, 11.75), (22.4, 11.75))

    # models and solvers
    box(ax, 1, 0.8, 98, 5, BASE_BG)
    ax.text(
        50,
        3.3,
        "Models: COBRApy, REFRAMED, GERM and kinetic (ODE)        Solvers: CPLEX, Gurobi, SCIP, GLPK",
        ha="center",
        va="center",
        fontsize=11,
        weight="bold",
    )

    fig.savefig(path, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    draw(Path(__file__).with_name("mewpy-arch.png"))
