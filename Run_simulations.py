from PyLTSpice import SimRunner
from PyLTSpice import SpiceEditor
from PyLTSpice import LTspice, RawRead
import os
import numpy as np
# import calculate_index as calc

def main():
    directory = './Netlists_new'
    
    folders = [folder for folder in os.listdir(directory) if os.path.isdir(os.path.join(directory, folder))]
    print(folders)
    for dirs in folders:
        path = './Netlists_new/{}'.format(dirs)
        new_dir = "FR_" + dirs
        os.makedirs(new_dir, exist_ok=True)
        # running simulation for each netlist file
                    
        for file in os.listdir(path):
            if os.path.isfile(os.path.join(path, file)) and file.endswith(".asc"):
                runner = SimRunner(output_folder='./temp_batch', simulator=LTspice)
                netlist_file = os.path.join(path, file)
                netlist = SpiceEditor(netlist_file)
                raw, log = runner.run_now(netlist, run_filename=file)
                raw_read = RawRead(raw)
                
                # takes the magnitude and freq series
                vout = np.abs(raw_read.get_wave('V(n61p1)'))
                vout = np.reshape(vout, (len(vout), 1))
                freq = np.abs(raw_read.get_trace('frequency'))
                out = np.concatenate((vout, freq[:, np.newaxis]), axis=1)
                
                #saves the file 
                filename = "FR_{}.txt".format(file)
                file_path = os.path.join(new_dir, filename)
                
                with open(file_path, "w") as nt:
                    np.savetxt(nt, out)
            break
                    


if __name__ == '__main__':
    main()
    # calc.calc()
