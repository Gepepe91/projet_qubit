import numpy as np
import matplotlib.pyplot as plt

# Charger les données depuis le fichier CSV
data_oscillating_with_noise = np.loadtxt('simulation_results.csv', delimiter=',', skiprows=1)
data_dynamic = np.loadtxt('dynamic_data.csv', delimiter=' ', skiprows=1)

# Extraire les colonnes pertinentes
noise_levels = np.unique(data_oscillating_with_noise[:, 0])
times = data_oscillating_with_noise[:, 1]

# Initialiser les listes pour stocker les différences
differences = []

plt.figure(figsize=(12, 6))

# Calculer les différences et tracer les résultats pour chaque niveau de bruit
for noise_level in noise_levels:
    indices = np.where(data_oscillating_with_noise[:, 0] == noise_level)[0]
    abs_alpha_corrige = data_oscillating_with_noise[indices, 4]
    abs_beta_corrige = data_oscillating_with_noise[indices, 5]

    # Calculer les différences par rapport aux valeurs théoriques
    difference_alpha_corrige = np.abs(abs_alpha_corrige - data_dynamic[:len(abs_alpha_corrige), -2])
    difference_beta_corrige = np.abs(abs_beta_corrige - data_dynamic[:len(abs_beta_corrige), -1])

    # Regrouper les différences alpha et beta
    combined_differences = (difference_alpha_corrige + difference_beta_corrige) / 2
    differences.append(combined_differences)

    # Calculer l'écart moyen toutes les 10 itérations
    window_size = 10
    mean_combined_differences = np.array([np.mean(combined_differences[max(0, i-window_size+1):i+1]) for i in range(len(combined_differences))])

    # Tracer les résultats
    plt.plot(times[indices], mean_combined_differences, label=f'Noise Level {noise_level}')

    # Identifier les moments où les courbes dépassent le seuil de 10%
    threshold = 0.1
    #exceeds_threshold = mean_combined_differences > threshold
    #plt.scatter(times[indices][exceeds_threshold], mean_combined_differences[exceeds_threshold], color='red', zorder=5)

plt.axhline(y=threshold, color='r', linestyle='--', label='10% Threshold')
plt.xlabel('Time')
plt.ylabel('Mean Difference')
plt.title('Mean Difference from Theoretical Values Over Time')
plt.legend()
plt.grid(True)
plt.show()

# Calculer et afficher les statistiques supplémentaires
plt.figure(figsize=(10, 6))

for i, noise_level in enumerate(noise_levels):
    mean_diff = np.mean(differences[i])
    std_diff = np.std(differences[i])
    plt.errorbar(noise_level, mean_diff, yerr=std_diff, fmt='o', label=f'Noise Level {noise_level}')

plt.axhline(y=threshold, color='r', linestyle='--', label='10% Threshold')
plt.xlabel('Noise Level')
plt.ylabel('Mean Difference')
plt.title('Mean and Standard Deviation of Differences by Noise Level')
plt.legend()
plt.grid(True)
plt.show()
