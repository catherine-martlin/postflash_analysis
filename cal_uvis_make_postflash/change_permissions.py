import os
import glob

def change_permissions(path_to_files):
    """Change permissions of the output files and directories.
        Code structure borrowed from cal_uvis_make_darks/FileIO.py.
        
        Parameters
        ----------
            path_to_files : str
                Path to the directory that contains the files you want
                to update.
    """
    
    os.chdir(path_to_files)
    all_directories = glob.glob('*')

    for directory in all_directories:
        try:
            os.chmod(directory, 0o775)
        except:
            print('Can not update. Onto Next.')
        for root, subdirs, files in os.walk(directory):
            for subdir in subdirs:
                try:
                    os.chmod(os.path.join(root, subdir), 0o775)
                except:
                    print('Can not update. Onto Next.')
            for name in files:
                try:
                    os.chmod(os.path.join(root, name), 0o775)
                except:
                    print('Can not update. Onto Next.')
