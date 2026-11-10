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

import copy

from reframed.core.cbmodel import CBModel


class SMomentModel(CBModel):

    def __init__(self, model, enzyme_reaction_prefix="R_ENZYME_DELIVERY_"):
        """SMomentModel an extension of REFRAMED CBModel for sMOMENT models.

        :param model: a path to a CBM xml file or an instance of REFRAMED SBModel
        :param enzyme_reaction_prefix: str enzyme prefix, with or without the SBML reaction prefix (R_)
                   GECKO_ANALOGON: "R_ENZYME_DELIVERY_"
                   GECKO: "R__TG_ER_"
        """
        # TODO: Define as an extension of simulator to be compatible with both
        # COBRApy and REFRAMED.

        if isinstance(model, str):
            from reframed.io.sbml import load_cbmodel

            model = load_cbmodel(model)
        elif not isinstance(model, CBModel):
            raise ValueError("The model should be a path or a CBModel")

        super(SMomentModel, self).__init__(model.id)
        # REFRAMED reaction ids carry the SBML prefix (R_): a prefix given without it, e.g. "ENZYME_DELIVERY_",
        # is completed by the simulator, as Simulator.resolve does for ids. The default ids are unchanged.
        if not any(rx.startswith(enzyme_reaction_prefix) for rx in model.reactions):
            from mewpy.simulation import get_simulator

            native_prefix = get_simulator(model).to_native(enzyme_reaction_prefix, "reaction")
            if any(rx.startswith(native_prefix) for rx in model.reactions):
                enzyme_reaction_prefix = native_prefix
        self.enzyme_prefix = enzyme_reaction_prefix
        # import CBModel's data
        self.compartments = copy.deepcopy(model.compartments)
        self.metabolites = copy.deepcopy(model.metabolites)
        self.reactions = copy.deepcopy(model.reactions)
        self.genes = copy.deepcopy(model.genes)
        self.enzymes = []
        for rx in self.reactions:
            if rx.startswith(self.enzyme_prefix):
                self.enzymes.append(rx[len(self.enzyme_prefix) :])

    @property
    def protein_rev_reactions(self):
        return {}

    @property
    def proteins(self):
        return self.enzymes

    def set_pool_bounds(self, lb, up):
        pass
