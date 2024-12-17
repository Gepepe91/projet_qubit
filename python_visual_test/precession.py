import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# Création des données pour l'animation
n_frames = 100  # Nombre de frames
theta = np.linspace(0, 2 * np.pi, n_frames)  # Angle azimutal (précession)
phi = np.pi / 4 * np.sin(theta) + np.pi / 4  # Angle polaire (oscillation)

# Conversion sphériques -> cartésiennes pour la flèche
x = np.sin(phi) * np.cos(theta)  # Coordonnée x de l'extrémité de la flèche
y = np.sin(phi) * np.sin(theta)  # Coordonnée y de l'extrémité de la flèche
z = np.cos(phi)                 # Coordonnée z de l'extrémité de la flèche

# Préparation de la figure et des axes 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.set_xlim([-1.1, 1.1])
ax.set_ylim([-1.1, 1.1])
ax.set_zlim([-1.1, 1.1])
ax.set_box_aspect([1, 1, 1])  # Égalisation des proportions

# Ajout de la sphère pour guider l'animation
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x_sphere = np.outer(np.cos(u), np.sin(v))
y_sphere = np.outer(np.sin(u), np.sin(v))
z_sphere = np.outer(np.ones_like(u), np.cos(v))
ax.plot_surface(x_sphere, y_sphere, z_sphere, color='lightblue', alpha=0.3)

# Initialisation de la flèche (vecteur)
arrow = ax.quiver(0, 0, 0, x[0], y[0], z[0], color='red', linewidth=2)

# Fonction de mise à jour de l'animation
def update(frame):
    # Met à jour la direction de la flèche
    arrow.set_segments([[[0, 0, 0], [x[frame], y[frame], z[frame]]]])
    return arrow,

# Création de l'animation
anim = FuncAnimation(fig, update, frames=n_frames, interval=50, blit=False)

# Afficher ou sauvegarder l'animation
plt.show()
anim.save("precession_arrow.gif", writer="ffmpeg", fps=30)
