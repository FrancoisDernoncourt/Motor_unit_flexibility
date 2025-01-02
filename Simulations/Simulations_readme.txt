The simulation script was run using Python 3.13.1

The .ipynb file requires an input .mat file to run if the variable excitatory_input_source is set to 'load_experimental_data'
The files used in the paper can be downloaded at https://e.pcloud.link/publink/show?code=kZVJnbZa8qAwdfN8883UGEs1IgJ0zMB5zRk#folder=12243207650&tpl=publicfoldergrid
Once downloaded, the variables 'path_for_input_files' and 'filename_for_input_file' should be modified to accomodate the downloaded file paths and file names.

The run_batch_simulations.ipynb script allows to run several iterations of a simulation with different input files (or with other changing parameters if the code is modified a bit).
For the batch simulation to work, the changing parameters in Simulation_script.ipynb should be modified to str(os.getenv('VAR_NAME')).