import scipy.io as sio
import numpy as np
import plotly.graph_objects as go

# Load the data from the .mat file
data = sio.loadmat('generated_data.mat')

omega_vect_1 = data['omega_vect_1'].flatten()
v_R_vect = data['v_R_vect'].flatten()
disp_3D = data['disp_3D']

# Since disp_3D is a 2D array, we don't need to slice it
Z = disp_3D

# Create a meshgrid for the contour plot
X, Y = np.meshgrid(omega_vect_1, v_R_vect)

# Define the levels for the contour plot
level_vect = np.linspace(0, 1, 31)

# Create the contour plot using Plotly
fig = go.Figure(data=go.Contour(
    z=Z,
    x=omega_vect_1,
    y=v_R_vect,
    contours=dict(
        start=level_vect[0],
        end=level_vect[-1],
        size=(level_vect[-1] - level_vect[0]) / len(level_vect),
        coloring='heatmap'
    ),
    colorbar=dict(
        title=dict(
            text='w(ω)',
            side='right'
        ),
        ticks='outside',
        tickvals=np.linspace(0, 1, 11),  # Fewer tick values to avoid clutter
        ticktext=[f'{val:.2f}' for val in np.linspace(0, 1, 11)],
        x=1.02,  # Position closer to the plot
        y=0.5,
        yanchor='middle',
        lenmode='fraction',
        len=1  # Length of the color bar
    )
))

# Update layout for the plot with smaller size and proper axis titles
fig.update_layout(
    title='Contour plot of FT of deflection',
    xaxis_title='Frequency, ω (rad/s)',
    yaxis_title='Velocity Ratio, α',
    font=dict(
        size=12
    ),
    margin=dict(l=70, r=10, t=70, b=70),  # Adjust the margins to ensure labels are fully visible
    width=600,  # Set the width of the figure
    height=400  # Set the height of the figure
)

# Display the plot
fig.show()
