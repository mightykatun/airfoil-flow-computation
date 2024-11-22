# Airfoil Flow Analysis with Wind Tunnel Validation

Welcome to the Airfoil Flow Analysis project repository! This software computes essential aerodynamic coefficients for airfoil profiles using a robust panel method approach in Python. It accurately predicts Cp (Coefficient of pressure), Cl (lift), and Cd (drag), and provides insightful visualizations of velocity streamlines.

### Features:
- **Panel Method Implementation**: Utilizes a subdivision technique for precise aerodynamic calculations.
- **Geometry Generation**: Includes an embedded NACA airfoil generator for versatile shape configurations.
- **Integration with XFoil**: Seamless integration with XFoil for enhanced accuracy and validation.
- **Streamline Visualization**: Generates graphical representations of flow streamlines for intuitive analysis.

<img src="https://github.com/user-attachments/assets/15564255-db2c-4705-9331-1a72bbaa2546" width="60%"/>

### Wind Tunnel Integration:
To complement the software's accuracy, a small-scale wind tunnel has been constructed. This setup includes:
- **Test Section Specifications**: 0.3m x 0.3m with a 9:1 contraction ratio.
- **Measurement Capabilities**: Arduino-based system for real-time data logging of lift, drag, and Cp polars.
- **Experimental Validation**: Validates software predictions against empirical data with a high degree of correlation (within 15% error margin).

For more details on the wind tunnel design and integration, visit [ldak.dev/projects/wind-tunnel](https://ldak.dev/projects/wind-tunnel).

### Future Developments:
Continued efforts are focused on:
- **Enhancing Computational Efficiency**: Optimizing code for faster calculations and expanded modeling capabilities.
- **Improving Accuracy**: Refining algorithms for more precise aerodynamic predictions across various airfoil shapes.

Feel free to explore and contribute!
