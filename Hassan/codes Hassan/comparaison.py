import numpy as np
import matplotlib.pyplot as plt

# Charger les données depuis le fichier CSV
data_oscillating_with_noise = np.loadtxt('simulation_results.csv', delimiter=',', skiprows=1)
data_dynamic = np.loadtxt('dynamic_data.csv', delimiter=' ', skiprows=1)

# Extraire les colonnes pertinentes
noise_levels = np.unique(data_oscillating_with_noise[:, 0])
times = data_oscillating_with_noise[:, 1]

plt.figure(figsize=(12, 6))

# Calculer les différences et tracer les résultats pour chaque niveau de bruit
for noise_level in noise_levels:
    indices = np.where(data_oscillating_with_noise[:, 0] == noise_level)[0]
    abs_alpha_corrige = data_oscillating_with_noise[indices, 4]
    abs_beta_corrige = data_oscillating_with_noise[indices, 5]

    # Assurez-vous que les données théoriques correspondent aux indices
    difference_alpha_corrige = np.abs(abs_alpha_corrige - data_dynamic[:len(abs_alpha_corrige), -2])
    difference_beta_corrige = np.abs(abs_beta_corrige - data_dynamic[:len(abs_beta_corrige), -1])

    # Calculer l'écart moyen toutes les 10 itérations
    mean_difference_alpha = np.convolve(difference_alpha_corrige, np.ones(10)/10, mode='valid')
    mean_difference_beta = np.convolve(difference_beta_corrige, np.ones(10)/10, mode='valid')

    # Tracer les résultats
    plt.plot(times[indices][:len(mean_difference_alpha)], mean_difference_alpha, label=f'Noise Level {noise_level} Alpha')
    plt.plot(times[indices][:len(mean_difference_beta)], mean_difference_beta, label=f'Noise Level {noise_level} Beta')

plt.axhline(y=0.2, color='r', linestyle='--', label='20% Threshold')
plt.xlabel('Time')
plt.ylabel('Mean Difference')
plt.title('Mean Difference from Theoretical Values Over Time')
plt.legend()
plt.grid(True)
plt.show()
