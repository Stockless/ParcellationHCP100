import os

def rename_files_in_directory(directory, find_text, replace_text):
    # List all files in the directory
    for filename in os.listdir(directory):
        # Create the new file name by replacing the text
        new_filename = filename.replace(find_text, replace_text)
        
        # Get the full path of the old and new file names
        old_file = os.path.join(directory, filename)
        new_file = os.path.join(directory, new_filename)
        
        # Rename the file if the name has changed
        if old_file != new_file:
            os.rename(old_file, new_file)
            print(f'Renamed: {filename} -> {new_filename}')

directory = os.path.join(os.getcwd(),"bundles")  # Replace with your directory path
find_text = '__'  # Text you want to replace
replace_text = '_'  # Text to replace with

rename_files_in_directory(directory, find_text, replace_text)