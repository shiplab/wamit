# HydroScripts

This is a set of several python scripts to pre and post-processing WAMIT outputs.

## Analysis Wamit

```python
import analysis_wamit

r = analysis_wamit.raos(plota=0, dof_plot=[1,2,3,4,5,6], inc_plot=[0,45,90,135,180])
```

```plota```    - flag to plot outputs

```dof_plot``` - degrees of freedom to plot

```inc_plot``` - headings to plot