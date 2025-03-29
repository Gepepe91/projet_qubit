import numpy as np
import matplotlib.pyplot as plt

# Charger les données depuis les fichiers CSV
data_oscillating_with_noise = np.loadtxt('simulation_results_repeated.csv', delimiter=',', skiprows=1)
data_dynamic = np.loadtxt('dynamic_data.csv', delimiter=' ', skiprows=1)

# Extraire les colonnes pertinentes
noise_levels = np.unique(data_oscillating_with_noise[:, 0])
times = np.unique(data_oscillating_with_noise[:, 2])

# Initialiser les listes pour stocker les différences et les répétitions
differences = []
repetition_counts = []

plt.figure(figsize=(12, 6))

# Calculer les différences et tracer les résultats pour chaque niveau de bruit
for noise_level in noise_levels:
    noise_data = data_oscillating_with_noise[data_oscillating_with_noise[:, 0] == noise_level]
    
    avg_diff_alpha = []
    avg_diff_beta = []
    repetitions = []
    
    for t in times:
        time_data = noise_data[noise_data[:, 2] == t]
        if time_data.shape[0] == 0:
            continue

        diff_alpha = np.abs(time_data[:, 5] - data_dynamic[:len(time_data), -2])
        diff_beta = np.abs(time_data[:, 6] - data_dynamic[:len(time_data), -1])
        
        avg_diff_alpha.append(np.mean(diff_alpha))
        avg_diff_beta.append(np.mean(diff_beta))
        repetitions.append(len(time_data))
    
    combined_differences = (np.array(avg_diff_alpha) + np.array(avg_diff_beta)) / 2
    differences.append(combined_differences)
    repetition_counts.append(repetitions)

    # Calculer l'écart moyen toutes les 10 itérations
    window_size = 10
    mean_combined_differences = np.array([np.mean(combined_differences[max(0, i-window_size+1):i+1]) for i in range(len(combined_differences))])

    # Tracer les résultats
    plt.plot(times[:len(mean_combined_differences)], mean_combined_differences, label=f'Noise Level {noise_level}')

plt.axhline(y=0.1, color='r', linestyle='--', label='10% Threshold')
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

plt.axhline(y=0.1, color='r', linestyle='--', label='10% Threshold')
plt.xlabel('Noise Level')
plt.ylabel('Mean Difference')
plt.title('Mean and Standard Deviation of Differences by Noise Level')
plt.legend()
plt.grid(True)
plt.show()

# Afficher l'évolution du nombre de répétitions
plt.figure(figsize=(12, 6))
for i, noise_level in enumerate(noise_levels):
    plt.plot(times[:len(repetition_counts[i])], repetition_counts[i], label=f'Noise Level {noise_level}')

plt.xlabel('Time')
plt.ylabel('Repetition Count')
plt.title('Repetition Count Over Time by Noise Level')
plt.legend()
plt.grid(True)
plt.show()
