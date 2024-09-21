import matplotlib.pyplot as plt
import csv

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
    return data

# Function to create the histogram
def create_histogram(data):
    names = list(data.keys())
    counts = list(data.values())

    plt.figure(figsize=(10, 6))
    plt.bar(names, counts, color='skyblue')
    plt.xlabel('Elements')
    plt.ylabel('Number of Copies')
    plt.title('Histogram of Elements')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()

# Main function to execute the process
def main(filename):
    data = read_file_and_process(filename)
    create_histogram(data)

# Example usage
filename = 'atlas_info_final.txt'  # Replace with your file name
main(filename)