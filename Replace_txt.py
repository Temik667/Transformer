import os

# Specify the directory where the files are located
directory = 'Netlists_new/AX'

# List all files in the directory
file_list = os.listdir(directory)

# Filter the files to find those ending with ".py.txt"
py_files = [file for file in file_list if file.endswith(".net")]

# Rename the files by removing the ".py" extension
for file in py_files:
    new_name = file.replace(".net", ".asc")
    os.rename(os.path.join(directory, file), os.path.join(directory, new_name))
    print(f"Renamed: {file} to {new_name}")
