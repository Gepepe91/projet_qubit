import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Création des données d'exemple
x = np.linspace(0, 2 * np.pi, 100)  # Points x
y = np.sin(x)                       # Exemple de fonction y = f(x)

# Préparer la figure et l'axe
fig, ax = plt.subplots()
line, = ax.plot(x, y, label="y = f(x)")  # Initialisation de la courbe

# Définir les limites des axes
ax.set_xlim(x.min(), x.max())
ax.set_ylim(-1.5, 1.5)
ax.set_title("Animation de y = f(x)")
ax.legend()

# Fonction pour mettre à jour les données de la courbe
def update(frame):
    # Modifier y en fonction du temps (par exemple une onde se déplaçant)
    y = np.sin(x + 0.1 * frame)  # Exemple de mise à jour
    line.set_ydata(y)  # Met à jour les données de la courbe
    return line,

# Créer l'animation
n_frames = 100  # Nombre total de frames
anim = FuncAnimation(fig, update, frames=n_frames, interval=50, blit=True)

# Afficher ou sauvegarder l'animation
plt.show()

anim.save('animation.gif', writer='ffmpeg', fps=30)

