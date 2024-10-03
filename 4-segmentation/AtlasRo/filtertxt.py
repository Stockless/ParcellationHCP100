# Open file1 and read its lines into a set for fast lookup
with open('atlasfiltro.txt', 'r') as file1:
    file1_lines = set(line.strip() for line in file1)

# Open file2 and prepare file3 for writing the matching lines
with open('atlasinformation.txt', 'r') as file2, open('atlas_info.txt', 'w') as file3:
    for line in file2:
        # Split the line by spaces and get the first word
        first_word = line.split()[0]
        
        # If the first word is in file1's set of lines, write the line to file3
        if first_word in file1_lines:
            file3.write(line)