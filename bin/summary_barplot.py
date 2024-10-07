#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import sys

# Access command-line arguments
if len(sys.argv) > 1:
    # The first argument (index 0) is the script name
    # The second argument (index 1) is the batch
    batch = sys.argv[1]
    print("Batch:", batch)
else:
    print("Please provide a valid batch name.")

summary_dir = "/well/bag/users/vbr851/projects/agnostic_diagnostic_v2/output/" + batch + "/summary/"
file_path = summary_dir + "overall_summary.csv"

# Load the data from the CSV file
data = pd.read_csv(file_path)

# Set the Barcode column as the index
data.set_index("barcode", inplace=True)

# Filter the columns with the required prefix
columns_to_plot = [col for col in data.columns if col.startswith("reads_") and col not in ["reads_all", "reads_mapped"]]

# Define a custom color palette
custom_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22']

# Plotting the bar chart
data_to_plot = data[columns_to_plot]
data_to_plot.plot(kind='bar', width=1.2, colormap='Paired')

# Customizing the plot
plt.xlabel('Barcode')
plt.ylabel('Reads Mapped')
plt.title('Bar Plot of Reads Mapped to Viral References')
plt.legend([col.replace("reads_", "") for col in columns_to_plot], bbox_to_anchor=(1, 1), loc='upper left')

# Show the plot
plt.show()

plt.savefig(summary_dir + 'mapped_reads.png', bbox_inches='tight')

