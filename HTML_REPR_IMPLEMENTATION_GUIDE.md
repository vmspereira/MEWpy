# HTML Representation Implementation Guide

This guide shows how to add pandas-like HTML representations (`_repr_html_()` methods) to MEWpy classes for better Jupyter notebook display.

## Completed Classes

The following classes already have `_repr_html_()` methods implemented:

1. ✅ **SimulationResult** - `src/mewpy/simulation/simulation.py`
2. ✅ **ODEModel** - `src/mewpy/model/kinetic.py`
3. ✅ **CommunityModel** - `src/mewpy/com/com.py`

These use the new `html_repr.render_html_table()` helper function.

## Classes That Need `_repr_html_()` Methods

### Phase 1 Classes (Core)
- **Model** - `src/mewpy/germ/models/model.py` (has old HTML, needs update)
- **Reaction** - `src/mewpy/germ/variables/reaction.py`
- **Gene** - `src/mewpy/germ/variables/gene.py`
- **Metabolite** - `src/mewpy/germ/variables/metabolite.py`
- **Simulator** - `src/mewpy/simulation/simulation.py`

### Phase 2 Classes (Optimization & Regulatory)
- **Problem** - `src/mewpy/problems/problem.py`
- **Regulator** - `src/mewpy/germ/variables/regulator.py`
- **Interaction** - `src/mewpy/germ/variables/interaction.py`
- **EvaluationFunction** - `src/mewpy/optimization/evaluation/evaluator.py`

### Phase 3 Classes (Advanced)
- **Expression** - `src/mewpy/germ/algebra/expression.py` (has old HTML, needs update)
- **VariableContainer** - `src/mewpy/germ/lp/linear_containers.py`
- **ConstraintContainer** - `src/mewpy/germ/lp/linear_containers.py`
- **Environment** - `src/mewpy/simulation/environment.py` (has old HTML, needs update)
- **HistoryManager** - `src/mewpy/util/history.py`

## Implementation Pattern

### Step 1: Import the Helper

```python
def _repr_html_(self):
    """Pandas-like HTML representation for Jupyter notebooks."""
    from mewpy.util.html_repr import render_html_table
```

### Step 2: Build Rows List

Extract data from your `__repr__()` method logic and create a list of (label, value) tuples:

```python
    rows = []

    # Example: Simple attribute
    if self.name:
        rows.append(("Name", self.name))

    # Example: Computed value
    if hasattr(self, "reactions"):
        rows.append(("Reactions", str(len(self.reactions))))

    # Example: Formatted value
    if self.objective_value is not None:
        rows.append(("Objective", f"{self.objective_value:.6g}"))

    # Example: Conditional with direction
    if self.maximize:
        rows.append(("Direction", "Maximize"))

    # Example: Indented sub-items (use "  " prefix for label)
    if self.constraints:
        rows.append(("Constraints", str(len(self.constraints))))
        rows.append(("  Environment", str(env_count)))
        rows.append(("  Model", str(model_count)))
```

### Step 3: Return Rendered HTML

```python
    return render_html_table("Class Name: {self.id}", rows)
```

## Complete Example: Reaction Class

Here's a complete example for the Reaction class:

```python
def _repr_html_(self):
    """Pandas-like HTML representation for Jupyter notebooks."""
    from mewpy.util.html_repr import render_html_table

    rows = []

    # Name
    if hasattr(self, "name") and self.name and self.name != self.id:
        rows.append(("Name", self.name))

    # Equation
    try:
        equation = self.equation
        if len(equation) > 80:
            equation = equation[:77] + "..."
        rows.append(("Equation", equation))
    except:
        rows.append(("Equation", "<not available>"))

    # Bounds
    try:
        lb, ub = self.bounds
        rows.append(("Bounds", f"({lb:.4g}, {ub:.4g})"))
    except:
        pass

    # Reversibility
    try:
        reversible = "Yes" if self.reversible else "No"
        rows.append(("Reversible", reversible))
    except:
        pass

    # Boundary
    try:
        boundary = "Yes" if self.boundary else "No"
        rows.append(("Boundary", boundary))
    except:
        pass

    # GPR
    try:
        if self.gpr and not self.gpr.is_none:
            gpr_str = self.gpr.to_string()
            if len(gpr_str) > 50:
                gpr_str = gpr_str[:47] + "..."
            rows.append(("GPR", gpr_str))
    except:
        pass

    # Genes count
    try:
        gene_count = len(self.genes)
        if gene_count > 0:
            rows.append(("Genes", str(gene_count)))
    except:
        pass

    # Metabolites count
    try:
        met_count = len(self.metabolites)
        if met_count > 0:
            rows.append(("Metabolites", str(met_count)))
    except:
        pass

    # Compartments
    try:
        comps = self.compartments
        if comps:
            comp_str = ", ".join(comps) if len(comps) <= 3 else f"{len(comps)} compartments"
            rows.append(("Compartments", comp_str))
    except:
        pass

    return render_html_table(f"Reaction: {self.id}", rows)
```

## Tips for Implementation

1. **Mirror __repr__() Logic**: Your `_repr_html_()` should show the same information as `__repr__()`, just formatted as HTML

2. **Use try/except**: Wrap attribute access in try/except to handle missing attributes gracefully

3. **Format Numbers**: Use Python format strings for consistent number display:
   - `f"{value:.4g}"` for general numbers
   - `f"{value:.6g}"` for high precision

4. **Truncate Long Strings**: Keep values readable:
   ```python
   if len(text) > 80:
       text = text[:77] + "..."
   ```

5. **Handle Indentation**: Use `"  "` prefix for labels to create hierarchy:
   ```python
   rows.append(("Parameters", str(total)))
   rows.append(("  Constant", str(const_count)))
   rows.append(("  Variable", str(var_count)))
   ```

6. **Empty Values**: For indented items without values, use empty string:
   ```python
   rows.append(("  organism_1", ""))
   ```

7. **Test in Jupyter**: The HTML is designed for Jupyter notebooks. Test by displaying the object directly in a notebook cell.

## Testing

Create a test file to verify your implementation:

```python
# Test that HTML is generated
obj = YourClass(...)
html = obj._repr_html_()
assert html is not None
assert "YourClass" in html
assert "mewpy-table" in html
print(f"✓ HTML generated ({len(html)} chars)")
```

## Styling

The `render_html_table()` function provides consistent styling:
- **Header**: Green background (#2e7d32), white text, bold
- **Rows**: Alternating with hover effect (#f5f5f5)
- **Labels**: 30% width, bold, dark gray (#333)
- **Values**: 70% width, regular, medium gray (#666)
- **Indentation**: 24px left padding for sub-items
- **Font**: System fonts (Apple/Segoe UI/Roboto)
- **Size**: 12px body, 14px header

## Priority Order

Recommended implementation order:

1. **High Priority** (most frequently used in notebooks):
   - Reaction, Gene, Metabolite (Phase 1)
   - Problem, EvaluationFunction (Phase 2)

2. **Medium Priority**:
   - Model (update existing), Simulator (Phase 1)
   - Regulator, Interaction (Phase 2)

3. **Lower Priority**:
   - Expression (update existing), Environment (update existing) (Phase 3)
   - VariableContainer, ConstraintContainer, HistoryManager (Phase 3)

## Questions?

Refer to the implemented examples:
- **SimulationResult**: Complex with nested constraints
- **ODEModel**: Hierarchical parameters display
- **CommunityModel**: Variable length organism lists

All use the same pattern and helper function for consistency.
