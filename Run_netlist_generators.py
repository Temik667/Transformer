import os
import subprocess

def run_scripts_in_folder(folder_path):
    # Iterate over the files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith('.py'):
            script_path = os.path.join(folder_path, filename)  # Construct full path to script
            try:
                # Run the script using the `python` command
                subprocess.run(['python', script_path], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error executing script: {e}")

# Example usage
directory_path = 'generators_new/AX'

# Use the os.listdir() function to get a list of items in the directory
items = os.listdir(directory_path)

# Filter the items to only include directories
folders = [item for item in items if os.path.isdir(os.path.join(directory_path, item))]
for folder in folders:
    name = "generators_new/AX/{}".format(folder)
    run_scripts_in_folder(name)
