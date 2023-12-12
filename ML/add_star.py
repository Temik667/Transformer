import os

def add_star_to_first_line(file_path):
    with open(file_path, 'r+') as file:
        content = file.readlines()
        if content:
            content[0] = '*' + content[0]
            file.seek(0)
            file.writelines(content)

def add_star_to_first_line_in_directory(directory):
    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            add_star_to_first_line(file_path)


directory_path = 'netlists' 
add_star_to_first_line_in_directory(directory_path)