import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider

# Create a new figure with a black background
fig = plt.figure(figsize=(10, 8), facecolor='black')
ax = fig.add_subplot(projection='3d', facecolor='black')

# Set the limits of the axes
limit = 5
ax.set_xlim([-limit, limit])
ax.set_ylim([-limit, limit])
ax.set_zlim([-limit, limit])

# Draw the positive axes with arrows using quiver
arrow_length = limit
ax.quiver(0, 0, 0, arrow_length, 0, 0, color='cyan', arrow_length_ratio=0.1, linewidth=1) # X-axis
ax.quiver(0, 0, 0, 0, arrow_length, 0, color='cyan', arrow_length_ratio=0.1, linewidth=1) # Y-axis
ax.quiver(0, 0, 0, 0, 0, arrow_length, color='cyan', arrow_length_ratio=0.1, linewidth=1) # Z-axis

# Draw negative axes (optional)
neg_arrow_length = -limit
ax.quiver(0, 0, 0, neg_arrow_length, 0, 0, color='cyan', arrow_length_ratio=0.05, linestyle='--', linewidth=0.5) # -X-axis
ax.quiver(0, 0, 0, 0, neg_arrow_length, 0, color='cyan', arrow_length_ratio=0.05, linestyle='--', linewidth=0.5) # -Y-axis
ax.quiver(0, 0, 0, 0, 0, neg_arrow_length, color='cyan', arrow_length_ratio=0.05, linestyle='--', linewidth=0.5) # -Z-axis

# Remove axis ticks and labels
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
ax.set_xlabel(None)
ax.set_ylabel(None)
ax.set_zlabel(None)
ax.set_axis_off()

# Define the ellipse
a_ellipse = 3.5  # Semi-major axis
b_ellipse = 2  # Semi-minor axis
phi_ellipse_static = np.deg2rad(30)  # Static rotation of the ellipse
theta_ellipse = np.linspace(0, 2 * np.pi, 100)
x_ellipse_data = a_ellipse * np.cos(theta_ellipse) * np.cos(phi_ellipse_static) - b_ellipse * np.sin(theta_ellipse) * np.sin(phi_ellipse_static)
y_ellipse_data = a_ellipse * np.cos(theta_ellipse) * np.sin(phi_ellipse_static) + b_ellipse * np.sin(theta_ellipse) * np.cos(phi_ellipse_static)
z_ellipse_data = np.zeros_like(x_ellipse_data)
ellipse_line, = ax.plot(x_ellipse_data, y_ellipse_data, z_ellipse_data, color='yellow', linewidth=2, label='Ellipse')

# Define the circle
diameter_circle = 2 * a_ellipse
radius_circle = diameter_circle / 2
theta_circle = np.linspace(0, 2 * np.pi, 100)
x_circle_data = radius_circle * np.cos(theta_circle)
y_circle_data = radius_circle * np.sin(theta_circle)
z_circle_data = np.zeros_like(x_circle_data)
circle_line, = ax.plot(x_circle_data, y_circle_data, z_circle_data, color='salmon', linewidth=2, label='Circle')

# Initial angle for point P
initial_angle_deg = 45
initial_angle_rad = np.deg2rad(initial_angle_deg)
x_P_initial = radius_circle * np.cos(initial_angle_rad)
y_P_initial = radius_circle * np.sin(initial_angle_rad)
z_P_initial = 0
point_P = ax.scatter([x_P_initial], [y_P_initial], [z_P_initial], color='blue', s=50, label='Point P')
point_P_text = ax.text(x_P_initial + 0.5, y_P_initial + 1, z_P_initial, 'P', color='white', fontsize=12)

# Initial angle offset for point Q (controlled by the second slider now)
initial_angle_POQ_deg = 60
initial_angle_POQ_rad = np.deg2rad(initial_angle_POQ_deg)

# Initial angle for point Q
initial_angle_Q_rad = initial_angle_rad + initial_angle_POQ_rad
x_Q_initial = radius_circle * np.cos(initial_angle_Q_rad)
y_Q_initial = radius_circle * np.sin(initial_angle_Q_rad)
z_Q_initial = 0
point_Q = ax.scatter([x_Q_initial], [y_Q_initial], [z_Q_initial], color='lime', s=50, label='Point Q')
point_Q_text = ax.text(x_Q_initial + 0.5, y_Q_initial + 1, z_Q_initial, 'Q', color='white', fontsize=12)

# Set the viewing angle
ax.view_init(elev=30, azim=60)

plt.subplots_adjust(left=0.1, top=0.85) # Adjusted top to make more space for sliders

# Create the slider for the angle of P
angleP_slider_ax = fig.add_axes([0.01, 0.90, 0.15, 0.05])
P_slider = Slider(
    ax=angleP_slider_ax,
    label='Angle P (°)',
    valmin=0,
    valmax=359,
    valinit=initial_angle_deg,
    valstep=1
)

# Create the slider for the angle POQ
anglePOQ_slider_ax = fig.add_axes([0.01, 0.84, 0.15, 0.05])
POQ_slider = Slider(
    ax=anglePOQ_slider_ax,
    label='Angle POQ (°)',
    valmin=10,
    valmax=340,
    valinit=initial_angle_POQ_deg,
    valstep=1
)

phi_value_text = fig.text(0.17, 0.96, f"Angle P = {P_slider.val}°", color='white', fontsize=10, verticalalignment='top', horizontalalignment='left')
POQ_value_text = fig.text(0.17, 0.90, f"Angle POQ = {POQ_slider.val}°", color='white', fontsize=10, verticalalignment='top', horizontalalignment='left')

# Function to update the positions of P and Q, and the text
def update(val):
    angle_P_deg = P_slider.val
    angle_P_rad = np.deg2rad(angle_P_deg)
    angle_POQ_deg = POQ_slider.val
    angle_POQ_rad = np.deg2rad(angle_POQ_deg)

    # Update point P
    x_P = radius_circle * np.cos(angle_P_rad)
    y_P = radius_circle * np.sin(angle_P_rad)
    z_P = 0
    point_P._offsets3d = ([x_P], [y_P], [z_P])
    point_P_text.set_position((x_P + 0.5, y_P + 1, z_P))
    phi_value_text.set_text(rf'Angle P = {int(angle_P_deg)}°')

    # Update point Q based on the angle of P and the POQ angle
    angle_Q_rad = angle_P_rad + angle_POQ_rad
    x_Q = radius_circle * np.cos(angle_Q_rad)
    y_Q = radius_circle * np.sin(angle_Q_rad)
    z_Q = 0
    point_Q._offsets3d = ([x_Q], [y_Q], [z_Q])
    point_Q_text.set_position((x_Q + 0.5, y_Q + 1, z_Q))
    POQ_value_text.set_text(rf'Angle POQ = {int(angle_POQ_deg)}°')

    fig.canvas.draw_idle()

# Connect the slider to the update function
P_slider.on_changed(update)
POQ_slider.on_changed(update) # Connect the second slider as well

plt.show()