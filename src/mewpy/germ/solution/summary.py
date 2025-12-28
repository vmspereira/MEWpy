import pandas as pd


class Summary:
    """
    A summary of a Solution. This object contains the main information of the solution.
    """

    def __init__(
        self,
        inputs: pd.DataFrame = None,
        outputs: pd.DataFrame = None,
        objective: pd.DataFrame = None,
        df: pd.DataFrame = None,
        metabolic: pd.DataFrame = None,
        regulatory: pd.DataFrame = None,
    ):
        """
        A summary of a Solution

        :param inputs: the inputs of the model
        :param outputs: the outputs of the model
        :param objective: the objective of the model
        :param df: the data frame of the summary
        :param metabolic: the metabolic summary
        :param regulatory: the regulatory summary
        :return:
        """
        if inputs is None:
            inputs = pd.DataFrame()

        if outputs is None:
            outputs = pd.DataFrame()

        if objective is None:
            objective = pd.DataFrame()

        if df is None:
            df = pd.DataFrame()

        if metabolic is None:
            metabolic = pd.DataFrame()

        if regulatory is None:
            regulatory = pd.DataFrame()

        self.inputs = inputs
        self.outputs = outputs
        self.objective = objective
        self.df = df
        self.metabolic = metabolic
        self.regulatory = regulatory

    def __repr__(self):
        """Text representation of the summary."""
        lines = []
        lines.append("=" * 60)
        lines.append("Solution Summary")
        lines.append("=" * 60)

        # Objective section
        if not self.objective.empty:
            lines.append("\nObjective:")
            lines.append(str(self.objective))

        # Inputs section
        if not self.inputs.empty:
            lines.append(f"\nInputs: {len(self.inputs)} entries")
            if len(self.inputs) <= 5:
                lines.append(str(self.inputs))
            else:
                lines.append("(Use .inputs to view full table)")

        # Outputs section
        if not self.outputs.empty:
            lines.append(f"\nOutputs: {len(self.outputs)} entries")
            if len(self.outputs) <= 5:
                lines.append(str(self.outputs))
            else:
                lines.append("(Use .outputs to view full table)")

        # Metabolic section
        if not self.metabolic.empty:
            lines.append(f"\nMetabolic: {len(self.metabolic)} reactions")
            if len(self.metabolic) <= 5:
                lines.append(str(self.metabolic))
            else:
                lines.append("(Use .metabolic to view full table)")

        # Regulatory section
        if not self.regulatory.empty:
            lines.append(f"\nRegulatory: {len(self.regulatory)} targets")
            if len(self.regulatory) <= 5:
                lines.append(str(self.regulatory))
            else:
                lines.append("(Use .regulatory to view full table)")

        # Full dataframe info
        if not self.df.empty:
            lines.append(f"\nFull Summary: {self.df.shape[0]} rows × {self.df.shape[1]} columns")
            lines.append("(Use .df to view full table)")

        lines.append("=" * 60)
        return "\n".join(lines)

    def _repr_html_(self):
        """
        It returns a html representation of the linear problem
        :return:
        """
        return self.df.to_html()
