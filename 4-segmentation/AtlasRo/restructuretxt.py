import re

def process_file_lines(file_path, output_path):
    # Regular expression for matching the first word in structure 1 and structure 2
    structure_1_pattern = r"^(.*?)_(.*?)_(.*?)_(\d+)$"  # Matches text1_text2_text3_number
    structure_2_pattern = r"^(.*?)_(.*?)_(.*?)$"        # Matches text1_text2_text3

    with open(file_path, 'r') as file:
        lines = file.readlines()

    with open(output_path, 'w') as output_file:
        for line in lines:
            # Split the line by tabs
            parts = line.strip().split('\t')
            first_word = parts[0]  # The first word to be transformed

            # Match against structure 1
            match_1 = re.match(structure_1_pattern, first_word)
            if match_1:
                text1, text2, text3, number = match_1.groups()
                # Transform first word: text1_text2number_text3number
                transformed_first_word = f"{text1}_{text2}{number}_{text3}{number}"
                parts[0] = transformed_first_word  # Replace the first word with the transformed one
                print(f"Transformed: {first_word} -> {transformed_first_word}")
            else:
                # Optionally handle structure 2 or leave it unchanged
                match_2 = re.match(structure_2_pattern, first_word)
                if match_2:
                    print(f"Structure 2 found (not changed): {first_word}")
                else:
                    print(f"Unknown structure: {first_word}")

            # Rebuild the line with the transformed first word
            new_line = '\t'.join(parts)
            output_file.write(new_line + "\n")

input_file = 'atlas_info.txt'   # Replace with your file path
output_file = 'atlas_info_filtrado.txt' # Replace with output file path
process_file_lines(input_file, output_file)
