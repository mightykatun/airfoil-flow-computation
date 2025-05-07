# Airfoil Flow Analysis with Wind Tunnel Validation

Welcome to the Airfoil Flow Analysis project repository! This software computes essential aerodynamic coefficients for airfoil profiles using a robust panel method approach in Python. It accurately predicts Cp (Coefficient of pressure), Cl (lift), and Cd (drag), and provides insightful visualizations of velocity streamlines.

### Features:
- **Panel Method Implementation**: Utilizes a subdivision technique for precise aerodynamic calculations.
- **Geometry Generation**: Includes an embedded NACA airfoil generator for versatile shape configurations.
- **Integration with XFoil**: Seamless integration with XFoil for enhanced accuracy and validation.
- **Streamline Visualization**: Generates graphical representations of flow streamlines for intuitive analysis.

<p align="center">
<img src="https://github.com/user-attachments/assets/15564255-db2c-4705-9331-1a72bbaa2546" width="60%"/>
</p>

### Wind Tunnel Integration:
To complement the software's accuracy, a small-scale wind tunnel has been constructed. This setup includes:
- **Test Section Specifications**: 0.3m x 0.3m with a 9:1 contraction ratio.
- **Measurement Capabilities**: Arduino-based system for real-time data logging of lift, drag, and Cp polars.
- **Experimental Validation**: Validates software predictions against empirical data with a high degree of correlation (within 15% error margin).

For more details on the wind tunnel design and integration, visit [ldak.dev/projects/wind-tunnel](https://ldak.dev/projects/wind-tunnel).

### Usage
*Note that to run on Linux, you need to have Wine installed (see `xfoil.py`)*

All the required configuration is found in `config.txt`:

```yaml
# foil computation
naca_foil: 4412			# 4-digit NACA complient code
angle_of_attack: 0
panel_number: 500		# number of panels used for approximation
v_infinity: 1			# air velocity
grid_size_x: 150
grid_size_y: 150
streamline_comp: True	# set to False to skip expensive and slow streamline computation

# show plot options (they will be saved either way)
foilgen_plot: False
panel_plot: False
cp_plot: False
cp_comparaison: False
streamline_plot: False

# xfoil args
xfoil_run: True
# default values (no need to change them unless you know what you're doing)
panel_bunching: 4
te/le_density: 1
panel_density: 1
top_y/c_lim: 1 1
bottom_y/c_lim: 1 1
```

### Future Developments:
Continued efforts are focused on:
- **Enhancing Computational Efficiency**: Optimizing code for faster calculations and expanded modeling capabilities.
- **Improving Accuracy**: Refining algorithms for more precise aerodynamic predictions across various airfoil shapes.

Feel free to explore and contribute!
