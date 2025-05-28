import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider

# Create a new figure with a black background
fig = plt.figure(figsize=(10, 8), facecolor='black')
ax = fig.add_subplot(projection='3d', facecolor='black')

# Set the limits of the axes
limit = 5
ax.set_xlim([-limit, limit])
ax.set_ylim([-limit, limit])
ax.set_zlim([-limit, limit + 2])  # Increased z limit to see the second ellipse

# Draw the positive axes with arrows using quiver
arrow_length = limit
ax.quiver(0, 0, 0, arrow_length, 0, 0, color='cyan', arrow_length_ratio=0.1, linewidth=1)  # X-axis
ax.quiver(0, 0, 0, 0, arrow_length, 0, color='cyan', arrow_length_ratio=0.1, linewidth=1)  # Y-axis
ax.quiver(0, 0, 0, 0, 0, arrow_length, color='cyan', arrow_length_ratio=0.1, linewidth=1)  # Z-axis

# Draw negative axes (optional)
neg_arrow_length = -limit
ax.quiver(0, 0, 0, neg_arrow_length, 0, 0, color='cyan', arrow_length_ratio=0.05, linestyle='--', linewidth=0.5)  # -X-axis
ax.quiver(0, 0, 0, 0, neg_arrow_length, 0, color='cyan', arrow_length_ratio=0.05, linestyle='--', linewidth=0.5)  # -Y-axis
ax.quiver(0, 0, 0, 0, neg_arrow_length, 0, color='cyan', arrow_length_ratio=0.05, linestyle='--', linewidth=0.5)  # -Z-axis

# Remove axis ticks and labels
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
ax.set_xlabel(None)
ax.set_ylabel(None)
ax.set_zlabel(None)
ax.set_axis_off()

# Define the first ellipse (in XY-plane, z=0)
a_ellipse = 3.5  # Semi-major axis
b_ellipse = 2  # Semi-minor axis
phi_ellipse_rad = np.deg2rad(30)  # Rotation around Z-axis
theta_ellipse = np.linspace(0, 2 * np.pi, 100)
x_ellipse_data = a_ellipse * np.cos(theta_ellipse) * np.cos(phi_ellipse_rad) - b_ellipse * np.sin(
    theta_ellipse) * np.sin(phi_ellipse_rad)
y_ellipse_data = a_ellipse * np.cos(theta_ellipse) * np.sin(phi_ellipse_rad) + b_ellipse * np.sin(
    theta_ellipse) * np.cos(phi_ellipse_rad)
z_ellipse_data = np.zeros_like(x_ellipse_data)
ellipse_line1, = ax.plot(x_ellipse_data, y_ellipse_data, z_ellipse_data, color='yellow', linewidth=2,
                         label='Ellipse Z=0')

# Define the second ellipse (parallel to XY-plane, z=7)
z_constant_ellipse2 = 7
x_ellipse_data_2 = a_ellipse * np.cos(theta_ellipse) * np.cos(phi_ellipse_rad) - b_ellipse * np.sin(
    theta_ellipse) * np.sin(phi_ellipse_rad)
y_ellipse_data_2 = a_ellipse * np.cos(theta_ellipse) * np.sin(phi_ellipse_rad) + b_ellipse * np.sin(
    theta_ellipse) * np.cos(phi_ellipse_rad)
z_ellipse_data_2 = np.full_like(x_ellipse_data_2, z_constant_ellipse2)
ellipse_line2, = ax.plot(x_ellipse_data_2, y_ellipse_data_2, z_ellipse_data_2, color='lime', linewidth=2,
                         label='Ellipse Z=7')

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

# Initial angle offset for point Q
initial_angle_POQ_deg = 150
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

plt.subplots_adjust(left=0.1, top=0.85)  # Adjusted top to make more space for sliders

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
    valmin=0,
    valmax=359,
    valinit=initial_angle_POQ_deg,
    valstep=1
)

phi_value_text = fig.text(0.17, 0.96, f"Angle P = {P_slider.val}°", color='white', fontsize=10,
                          verticalalignment='top', horizontalalignment='left')
POQ_value_text = fig.text(0.17, 0.90, f"Angle POQ = {POQ_slider.val}°", color='white', fontsize=10,
                          verticalalignment='top', horizontalalignment='left')

# Store the line segments and points for the perpendiculars and intersections
perp_P_line = None
perp_Q_line = None
point_D = None
point_E = None
text_D = None
text_E = None
point_D_prime = None  # To store the projection on the second ellipse
text_D_prime = None  # Text for the projection of D
point_E_prime = None  # To store the projection of E
text_E_prime = None  # Text for the projection of E

# Lists to store the trail of lines
trail_D_prime_E_lines = []
trail_E_prime_D_lines = []

# Function to update the positions of P and Q, and the text
def update(val):
    global perp_P_line, perp_Q_line, point_D, point_E, text_D, text_E, point_D_prime, text_D_prime, point_E_prime, text_E_prime, trail_D_prime_E_lines, trail_E_prime_D_lines

    angle_P_deg = P_slider.val
    angle_P_rad = np.deg2rad(angle_P_deg)
    angle_POQ_deg = POQ_slider.val
    angle_POQ_rad = np.deg2rad(angle_POQ_deg)
    phi_ellipse_rad_current = phi_ellipse_rad  # Unghiul de rotatie al elipsei este constant aici

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

    # Coordonatele a două puncte de pe axa mare a elipsei
    x_major_axis_1 = a_ellipse * np.cos(phi_ellipse_rad_current)
    y_major_axis_1 = a_ellipse * np.sin(phi_ellipse_rad_current)
    z_major_axis_1 = 0
    x_major_axis_2 = -a_ellipse * np.cos(phi_ellipse_rad_current)
    y_major_axis_2 = -a_ellipse * np.sin(phi_ellipse_rad_current)
    z_major_axis_2 = 0

    # Desenarea perpendicularei din P pe axa mare (aproximativ)
    x_perp_P_end = (x_P + y_P * np.tan(phi_ellipse_rad_current)) / (1 + np.tan(phi_ellipse_rad_current) ** 2) * np.cos(
        phi_ellipse_rad_current)
    y_perp_P_end = (x_P + y_P * np.tan(phi_ellipse_rad_current)) / (1 + np.tan(phi_ellipse_rad_current) ** 2) * np.sin(
        phi_ellipse_rad_current)

    if perp_P_line:
        perp_P_line.remove()
    perp_P_line, = ax.plot([x_P, x_perp_P_end],
                            [y_P, y_perp_P_end],
                            [z_P, z_P], linestyle='--', color='gray')

    # Desenarea perpendicularei din Q pe axa mare (aproximativ)
    x_perp_Q_end = (x_Q + y_Q * np.tan(phi_ellipse_rad_current)) / (1 + np.tan(phi_ellipse_rad_current) ** 2) * np.cos(
        phi_ellipse_rad_current)
    y_perp_Q_end = (x_Q + y_Q * np.tan(phi_ellipse_rad_current)) / (1 + np.tan(phi_ellipse_rad_current) ** 2) * np.sin(
        phi_ellipse_rad_current)

    if perp_Q_line:
        perp_Q_line.remove()
    perp_Q_line, = ax.plot([x_Q, x_perp_Q_end],
                            [y_Q, y_perp_Q_end],
                            [z_Q, z_Q], linestyle='--', color='gray')

    # APROXIMARE pentru punctele D și E (intersecția cu elipsa Z=0)
    factor_offset = 1.2  # Un factor pentru a plasa etichetele vizual în apropiere
    x_D_approx = (x_P + x_perp_P_end) / 2 * factor_offset
    y_D_approx = (y_P + y_perp_P_end) / 2 * factor_offset
    x_E_approx = (x_Q + x_perp_Q_end) / 2 * factor_offset
    y_E_approx = (y_Q + y_perp_Q_end) / 2 * factor_offset

    if point_D:
        point_D.remove()
    point_D = ax.scatter([x_D_approx], [y_D_approx], [z_P], color='red', s=30)
    if text_D:
        text_D.remove()
    text_D = ax.text(x_D_approx, y_D_approx, z_P, 'D', color='white')

    if point_E:
        point_E.remove()
    point_E = ax.scatter([x_E_approx], [y_E_approx], [z_Q], color='orange', s=30)
    if text_E:
        text_E.remove()
    text_E = ax.text(x_E_approx, y_E_approx, z_Q, 'E', color='white')

    # Projection of D onto the second ellipse (Z=7)
    x_D_prime = x_D_approx
    y_D_prime = y_D_approx
    z_D_prime = 7

    if point_D_prime:
        point_D_prime.remove()
    point_D_prime = ax.scatter([x_D_prime], [y_D_prime], [z_D_prime], color='cyan', s=30)
    if text_D_prime:
        text_D_prime.remove()
    text_D_prime = ax.text(x_D_prime, y_D_prime, z_D_prime, 'D\'', color='white')

    # Projection of E onto the second ellipse (Z=7)
    x_E_prime = x_E_approx
    y_E_prime = y_E_approx
    z_E_prime = 7

    if point_E_prime:
        point_E_prime.remove()
    point_E_prime = ax.scatter([x_E_prime], [y_E_prime], [z_E_prime], color='magenta', s=30)
    if text_E_prime:
        text_E_prime.remove()
    text_E_prime = ax.text(x_E_prime, y_E_prime, z_E_prime, 'E\'', color='white')

    # Draw new line from D' to E and add to trail
    new_line_D_prime_E, = ax.plot([x_D_prime, x_E_approx], [y_D_prime, y_E_approx], [z_D_prime, z_Q], color='pink', linestyle='--')
    trail_D_prime_E_lines.append(new_line_D_prime_E)

    # Draw new line from E' to D and add to trail
    new_line_E_prime_D, = ax.plot([x_E_prime, x_D_approx], [y_E_prime, y_D_approx], [z_E_prime, z_P], color='pink', linestyle='--')
    trail_E_prime_D_lines.append(new_line_E_prime_D)

    fig.canvas.draw_idle()

def on_poq_slider_change(val):
    global trail_D_prime_E_lines, trail_E_prime_D_lines
    for line in trail_D_prime_E_lines:
        line.remove()
    for line in trail_E_prime_D_lines:
        line.remove()
    trail_D_prime_E_lines = []
    trail_E_prime_D_lines = []
    update(val)


# Connect the slider to the update function
P_slider.on_changed(update)
POQ_slider.on_changed(on_poq_slider_change)

plt.show()