import os
import re

def rename_files_with_new_structure(directory):
    # Regular expression for matching filenames
    structure_1_pattern = r"^(.*?)_(.*?)_(.*?)_(\d+)(\..+)$"  # Matches text1_text2_text3_number.extension
    structure_2_pattern = r"^(.*?)_(.*?)_(.*?)(\..+)$"        # Matches text1_text2_text3.extension

    for filename in os.listdir(directory):
        old_file = os.path.join(directory, filename)
        
        # Match against structure 1
        match_1 = re.match(structure_1_pattern, filename)
        if match_1:
            text1, text2, text3, number, extension = match_1.groups()
            # Create new name: text1_text2number_text3number.extension
            new_filename = f"{text1}_{text2}{number}_{text3}{number}{extension}"
            new_file = os.path.join(directory, new_filename)
            os.rename(old_file, new_file)
            print(f"Renamed: {old_file} -> {new_file}")
        else:
            # Optionally handle structure 2 or skip it
            match_2 = re.match(structure_2_pattern, filename)
            if match_2:
                print(f"Structure 2 found (not renamed): {old_file}")
            else:
                print(f"Unknown structure: {old_file}")

directory = os.path.join(os.getcwd(),"bundles")  # Replace with your directory path
rename_files_with_new_structure(directory)
