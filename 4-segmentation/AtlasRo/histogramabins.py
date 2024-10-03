import matplotlib.pyplot as plt
import csv
from collections import defaultdict

# Function to read the file and process the data
def read_file_and_process(filename):
    data = {}
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            name = row[0]
            count = int(row[2])
            if name in data:
                data[name] += count
            else:
                data[name] = count
    for i in data:
        if data[i] < 20:
            print(i)
    return data

# Function to group data into intervals and sum the counts
def group_data(data, interval_size):
    interval_data = defaultdict(int)
    for name, count in data.items():
        interval_start = (count // interval_size) * interval_size
        interval_end = interval_start + interval_size - 1
        interval_label = f"[{interval_start}-{interval_end}]"
        interval_data[interval_label] += count
    return interval_data

# Function to create the histogram
def create_histogram(data, interval_size):
    grouped_data = group_data(data,interval_size)
    intervals = sorted(grouped_data.keys(), key=lambda x: int(x.split('-')[0][1:]))
    counts = [grouped_data[interval] for interval in intervals]

    plt.figure(figsize=(10, 6))
    bars = plt.bar(intervals, counts, color='skyblue')
    plt.xlabel('Intervalos de tamaños de fascículos')
    plt.ylabel('Número de fibras')
    plt.title('Histograma agrupado en intervalos sin filtrar')
    plt.xticks(rotation=90)

    # Add text annotations above the bars
    for bar, count in zip(bars, counts):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(count),
                 ha='center', va='bottom')

    plt.tight_layout()
    plt.show()

# Main function to execute the process
def main(filename,interval_size):
    data = read_file_and_process(filename)
    create_histogram(data,interval_size)

# Example usage
filename = 'atlas_info_final.txt'  # Replace with your file name
interval_size = 1 # Replace with desired interval size
main(filename, interval_size)