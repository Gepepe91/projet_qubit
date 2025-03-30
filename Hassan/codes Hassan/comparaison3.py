import numpy as np
import matplotlib.pyplot as plt

# Load the data from the CSV files
data_oscillating_with_noise = np.loadtxt('simulation_results_repeated.csv', delimiter=',', skiprows=1)
data_dynamic = np.loadtxt('dynamic_data.csv', delimiter=' ', skiprows=1)

# Verify the shape of the loaded data
print("Shape of data_oscillating_with_noise:", data_oscillating_with_noise.shape)
print("Shape of data_dynamic:", data_dynamic.shape)

# Extract unique noise levels and times
noise_levels = np.unique(data_oscillating_with_noise[:, 0])
times = np.unique(data_oscillating_with_noise[:, 2])

# Initialize lists to store the averaged differences
averaged_differences = []

plt.figure(figsize=(12, 6))

# Calculate the averaged differences and plot the results for each noise level
for noise_level in noise_levels:
    # Filter data for the current noise level
    noise_data = data_oscillating_with_noise[data_oscillating_with_noise[:, 0] == noise_level]

    # Initialize arrays to store the sum of differences for averaging
    sum_difference_alpha_corrige = np.zeros(len(times))
    sum_difference_beta_corrige = np.zeros(len(times))
    repetition_count = np.zeros(len(times))

    # Iterate over each time step
    for i, t in enumerate(times):
        # Filter data for the current time step
        time_data = noise_data[noise_data[:, 2] == t]

        # Check if time_data has the expected number of columns
        if time_data.shape[1] < 6:
            raise ValueError(f"time_data does not have enough columns at time {t}. Expected at least 6 columns.")

        # Sum the differences for alpha and beta
        sum_difference_alpha_corrige[i] += np.sum(np.abs(time_data[:, 4] - data_dynamic[i, -2]))
        sum_difference_beta_corrige[i] += np.sum(np.abs(time_data[:, 5] - data_dynamic[i, -1]))

        # Count the number of repetitions for averaging
        repetition_count[i] = time_data.shape[0]

    # Calculate the average differences
    print(repetition_count)
    avg_difference_alpha_corrige = sum_difference_alpha_corrige / repetition_count
    avg_difference_beta_corrige = sum_difference_beta_corrige / repetition_count

    # Combine the differences for alpha and beta
    combined_avg_differences = (avg_difference_alpha_corrige + avg_difference_beta_corrige) / 2
    averaged_differences.append(combined_avg_differences)

    # Calculate the moving average of the combined differences
    window_size = 10
    mean_combined_avg_differences = np.array([np.mean(combined_avg_differences[max(0, i-window_size+1):i+1]) for i in range(len(combined_avg_differences))])

    # Plot the results
    plt.plot(times, mean_combined_avg_differences, label=f'Noise Level {noise_level}')

# Plot the 10% threshold line
threshold = 0.1
plt.axhline(y=threshold, color='r', linestyle='--', label='10% Threshold')
plt.xlabel('Time')
plt.ylabel('Mean Difference')
plt.title('Mean Difference from Theoretical Values Over Time')
plt.legend()
plt.grid(True)
plt.show()

# Plot the mean and standard deviation of differences by noise level
plt.figure(figsize=(10, 6))

for i, noise_level in enumerate(noise_levels):
    mean_diff = np.mean(averaged_differences[i])
    std_diff = np.std(averaged_differences[i])
    plt.errorbar(noise_level, mean_diff, yerr=std_diff, fmt='o', label=f'Noise Level {noise_level}')

plt.axhline(y=threshold, color='r', linestyle='--', label='10% Threshold')
plt.xlabel('Noise Level')
plt.ylabel('Mean Difference')
plt.title('Mean and Standard Deviation of Differences by Noise Level')
plt.legend()
plt.grid(True)
plt.show()
