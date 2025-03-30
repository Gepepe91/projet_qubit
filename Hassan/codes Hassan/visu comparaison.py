import numpy as np
import matplotlib.pyplot as plt

# Lire les fichiers statiques (sans nombres complexes)
data_static = np.loadtxt('data.csv', delimiter=' ', skiprows=1)
data_with_noise = np.loadtxt('data_with_noise.csv', delimiter=' ', skiprows=1)

# Lire les fichiers dynamiques (avec nombres complexes)
data_dynamic = np.loadtxt('dynamic_data.csv', delimiter=' ', skiprows=1)
data_oscillating_with_noise = np.loadtxt('oscillating_magnetic_field_with_noise_and_correction.csv', delimiter=' ', skiprows=1)

# Tracer les courbes
plt.figure(figsize=(10, 8))

dif_static = np.abs(data_static[:, -2] - data_with_noise[:, -2])
dif_dynamic = np.abs(data_dynamic[:, -2] - data_oscillating_with_noise[:, 1])

plt.subplot(2, 1, 1)
plt.plot(data_static[:, 0], data_static[:, -2], label='Static', color='b')
plt.plot(data_with_noise[:, 0], data_with_noise[:, -2], label='Static with Noise', color='r')
# plt.plot(data_static[:, 0], dif_static, label='Difference', color='g')
plt.xlabel('Time')
plt.ylabel('abs_alpha2')
plt.title('Static Field vs Static Field with Noise')
plt.legend()
plt.grid(True)

# Deuxième graphique : Dynamic Field vs Oscillating Field with Noise
plt.subplot(2, 1, 2)
plt.plot(data_dynamic[:, 0], data_dynamic[:, -2], label='Dynamic', color='b')
plt.plot(data_oscillating_with_noise[:, 0], data_oscillating_with_noise[:, 1], label='Oscillating with Noise', color='r')
plt.plot(data_oscillating_with_noise[:, 0], data_oscillating_with_noise[:, 3], label='alpha_corrige', color='g')
# plt.plot(data_oscillating_with_noise[:, 0], data_oscillating_with_noise[:, 5], label='alpha_corrige2', color='c')
# plt.plot(data_oscillating_with_noise[:, 0], data_oscillating_with_noise[:, 7], label='alpha_corrige2', color='m')
#plt.plot(data_oscillating_with_noise[:, 0], data_oscillating_with_noise[:, -2], label='bruit_alpha', color='m')
plt.xlabel('Time')
plt.ylabel('abs_alpha2')
plt.title('Dynamic Field vs Oscillating Field with Noise')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

# Nouvelle figure pour l'écart entre la correction et le signal bruité
plt.figure(figsize=(10, 6))

# Calculer l'écart entre les signaux corrigés et bruités
difference_alpha_corrige = np.abs(data_oscillating_with_noise[:, 1] - data_dynamic[:, -2])
difference_alpha_corrige2 = np.abs(data_oscillating_with_noise[:, 3] - data_dynamic[:, -2])

# Tracer les écarts
plt.plot(data_oscillating_with_noise[:, 0], difference_alpha_corrige, label='Ecart alpha_bruit', color='g')
plt.plot(data_oscillating_with_noise[:, 0], difference_alpha_corrige2, label='Ecart alpha_corrige', color='c')
plt.xlabel('Time')
plt.ylabel('Difference')
plt.title('Ecart entre la correction et le signal bruité')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
