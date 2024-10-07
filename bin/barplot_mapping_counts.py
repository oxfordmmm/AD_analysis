import pandas as pd
import matplotlib.pyplot as plt

# Load the data from the CSV file
data = pd.read_csv("overall_summary.csv")

# Set the Barcode column as the index
data.set_index("Barcode", inplace=True)

# Filter the columns with the required prefix
columns_to_plot = [col for col in data.columns if col.startswith("Reads_") and col not in ["Reads_all", "Reads_Mapped"]]

# Plotting the bar chart with adjusted bar width
data_to_plot = data[columns_to_plot]
ax = data_to_plot.plot(kind='bar', width=0.8, colormap='Paired')  # Adjust width as needed

# Customizing the plot
plt.xlabel('Barcode')
plt.ylabel('Values')
plt.title('Bar Plot from overall_summary.csv')
plt.legend([col.replace("Reads_", "") for col in columns_to_plot], bbox_to_anchor=(1, 1), loc='upper left')

# Save the plot as a PNG file
plt.savefig('bar_plot.png', bbox_inches='tight')

# Show the plot
plt.show()
