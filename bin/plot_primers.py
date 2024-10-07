#!/usr/bin/env python3
import matplotlib.pyplot as plt

# Assuming you have pairs of (start, end) positions of primer sequences along the linear DNA
primer_positions = [(144, 123), (174, 153), (205, 184), (236, 215)]

# Create a figure
plt.figure(figsize=(10, 3))

# Plot the start positions with blue vertical lines
plt.scatter([start for start, end in primer_positions], [1] * len(primer_positions), marker='|', color='blue', s=100, label='Start Position')

# Plot the end positions with red vertical lines
plt.scatter([end for start, end in primer_positions], [1] * len(primer_positions), marker='|', color='red', s=100, label='End Position')

plt.xlabel('Position along Linear DNA')
plt.yticks([])

# Add a legend to distinguish between start and end positions
plt.legend(loc='upper right')

plt.title('Primer Positions in Phi29-Amplified Linear DNA')
plt.show()

