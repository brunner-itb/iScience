<h1>
Documentation
</h1>
Because github pages are publicly accessible, 
I moved the documentation to the "doc" branch. 
To download the documentation either 

1. Download the zip file from the "doc" branch

2. clone the repo with

`git clone -b doc http:// https://personal_access_token@github.com/lukaskiwitz/thesis`

To clone the repo from the cli you will need to set up a personal acces token in your github settings. Then you can
download any update with

`git pull`

3. Alternatively, you can build the documentation yourself, provided that you have a properly configured conda
   environment. Move to the docsrc folder and run

`make html`


When done open `index.html` in your web browser to access the title page.
<h1>
    Requirements
</h1>

1. All python dependencies (provided as an anaconda yml-file)
2. correctly set LD_LIBRARY_PATH environment variable
3. installed thesisproject-package (provided as wheel file for pip)

<h1>
Installation Instructions (with anaconda)
</h1>

1. Install anaconda:
`https://www.anaconda.com/products/individual`


2. Run

`
conda env create -f fenics.yml -n "enviroment_name"
`

to create a new anaconda environment with the necessary packages.
This could take a moment.

3.Run

`conda activate "enviroment_name"`

to activate the newly created environment.

4.Install the package using pip with

`
pip install wheel-file
`

if you got the wheel file or 

`
pip install -e "repository_dir"
`

if you cloned this repository.

5.Run

`
conda env config vars set LD_LIBRARY_PATH=~/anaconda3/envs/enviroment_name/lib
`


and if you want to use paraview rendering 

`conda env config vars set PARAVIEW_PATH=~path_to_paraview_dir
`


to set up necessary environment variables.

6.Reactivate the conda environment to apply changes

`conda activate enviroment_name`


<h1>
    Usage
</h1>

For clarification, a fully documented `example_parameters.py` file has be created at `model/scripts/paper_models/example_parameters.py`. 
Refer to it if questions regarding the parameters and how to change them arise. \
If a parameter is not chosen its value will revert to its standard value chosen in the 
parameter template. These templates can be found at `model/scenarios/box_grid.py`, together with their unit conversion.

To replicate Figure 1 of the correspoding publication, follow these steps:
1. Navigate to `model/scripts/paper_models/Brunner_etal_Figure1/run`.
2. Depending on which panel you want to recreate, choose a folder, e.g. `B-E`.  Here you will find all necessary files to reproduce the data needed for the plotting.
3. Open the parameters.py file and edit the path variable at the top to point to the desired top level directoy for your simulation results.
4. Execute `python run_all.py && python list_post_process.py && python create_df.py`. Note: This is a parameter scan which will do multiple scans in parallel. How many are done simultaneously can be adjusted in the `run_all.py` script, line 33.
5. These scripts will first run the scan, perform the post_processing and create the data frames for each sub-run.
6. To combine the dataframes from the sub-runs, first edit the path variable in the `combine_dfs.py` file to point to your results, then execute it. This will then combine them into one `cell_df.hf` and `global_df.h5` each.
6. After this has finished, navigate to `model/scripts/paper_models/Brunner_etal_Figure1/plot/B`, adjust the paths according to the path you chose in the `parameters.py` file and execute both files here to create the plots of panel B.

The first time any model is run it will need to
- create a mesh
- mark the boundary parts
- run the FFC compiler

this may take some time. For large meshes this can be considerable (1h or so), but 
subsequent runs should be much faster. 

<h1>
    Paraview Integration
</h1>
Only works for simulation data generated with the current version

1. Download an omesa build of paraview (necessary to run without X-server)
2. Set the PARAVIEW_PATH env. variable to your paraview directory. For instance in you conda env. or in the post_process.py
4. pass the keyword argument "render_paraview=True" to run_post_process() 
5. various settings can be passed. "example_min" contains a commented example
   

The cell type markers are floats and will be set according to "SimContainer.marker_lookup".
SimContainer now has the "markers" attribute, which lists entity 
properties (parameters or python attributes) for which to export a MeshFunction. 
Assigment of non numerical values is controlled by marker_lookup.






 

 

