import numpy as np
import matplotlib.pyplot as plt


data_oscillating_with_noise = np.loadtxt('oscillating_magnetic_field_with_noise.csv', delimiter=' ', skiprows=1)


plt.hist(data_oscillating_with_noise[:, -2], bins=100, color='b', alpha=0.7, label='Bx')
plt.show()