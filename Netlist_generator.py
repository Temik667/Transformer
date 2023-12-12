import os

output_directory = "generators_new/RAD/Disk_5"  # Replace with the actual output directory path

for index in range(0, 30):
    step = round(0.5 + index*0.5, 2)
    name = "RAD_5_disk_{}.py".format(step)
    file_path = os.path.join(output_directory, name)

    # Read the template file
    with open("Sample generators/RAD/RAD_5_disk.py", "r") as template_file:  # Replace with the path to your template file
        template_code = template_file.read()

    # Replace the line "ro = D_wire" with "ro = {index} * D_wire"
    modified_code = template_code.replace("ro = 0", "ro = {}".format(step))
    modified_code = modified_code.replace("""filename = "netlist_12awg_10disk_RADDISP_3Diameter_AccInductance_complex_5disks.txt" """, 
                                          """filename = "{}.txt" """.format(name))

    # Create a new file with the modified code
    with open(file_path, "w") as file:
        file.write(modified_code)

    print(f"Created file: {file_path}")