import matplotlib.pyplot as plt
import csv

# Function to read the file and process the data
def read_file_and_process(filename):
    data = {}
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            name = row[0].replace("MNI_atlas","")
            count = int(row[2])
            if name in data:
                data[name] += count
            else:
                data[name] = count
    return data

# Function to create the histogram
def create_histogram(data):
    # Sort the data by count (values)
    sorted_data = sorted(data.items(), key=lambda item: item[1])

    names = [item[0] for item in sorted_data]
    counts = [item[1] for item in sorted_data]

    plt.figure(figsize=(10, 6))
    bars = plt.bar(names, counts, color='skyblue')
    plt.xlabel('Fascículos')
    plt.ylabel('Numero de fibras')
    plt.title('Histograma de fibras')
    plt.xticks(rotation=90)
    # Add text annotations above the bars
    for bar, count in zip(bars, counts):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(count),
                 ha='center', va='bottom')
    plt.tight_layout()
    plt.show()

# Main function to execute the process
def main(filename):
    data = read_file_and_process(filename)
    create_histogram(data)

# Example usage
filename = 'atlas_info.txt'  # Replace with your file name
main(filename)