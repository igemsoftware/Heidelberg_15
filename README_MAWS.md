# MAWS manual
## What is MAWS?
Have **you** ever wondered how you could detect that **one specific** protein or small molecule? Then you know that finding methods to detect novel compounds takes a lot of thought, time and money. Creating antibodies specific to your protein of interest is a costly endeavour, and selecting aptamers to your small molecules takes a lot of planning, a good nucleic acid library, weeks of work, and the possibility of failure due to the aptamer simply not being contained in the DNA pool. So would it not be great, if **you** could bypass all of those restrictions? If you could just take a **crystal structure** of your target, load it onto **your laptop**, push a **few keys**, and have results in **under a day**? 

MAWS is all of that!

The name of our software MAWS, is in itself a specification of exactly what you can do with it - Make Aptamers Without SELEX. And this guide will show you exactly what you have to do in a step-by-step fashion.

## Setting Things Up
In this section we will show you, what you need to do, to create an environment which allows MAWS to function.

### Installation
1. Check the specifications of your computer. MAWS was originally developed on a laptop with four cores and 12GB of RAM. Anything in that region should work pretty well, lower specifications could result in longer runtimes.
2. Download and install your favourite Linux distribution. We made MAWS on Arch Linux, but other distributions should run the same way.
3. Download and install your favourite Python distribution. Standard Python should be fine, the following steps will be shown using the "conda" package manager, which you can install with `pip install conda`.
4. Install "binstar" or the "anaconda-launcher" by running `conda install binstar` or `conda install anaconda-launcher`.
5. Add new channels where from you will download and install the needed packages with the command `conda config --add channels https://conda.binstar.org/omnia`.
6. Create an environment in which you will run MAWS. Just replace "iGEM" with the name you want for your environment: `conda create -n iGEM omnia openmm ambermini numpy mpmath`.
7. Activate your environment using `source activate iGEM`, replacing "iGEM" with the name of your environment.
8. Congratulations! You are now ready to run MAWS!

### Preprocessing
In order for MAWS to be able to understand your input files, it needs quite a bit of preprocessing. Here, we show you, how you can easily do that.

For proteins:
1. Check if your **.pdb** file contains one or more of the following:
  * Any solvent molecules. This includes:
    * Waters
    * Glycerol
    * DMSO
    * and many more ...
  * Hydrogens (you will probably not have them, but check nonetheless).
  * Small molecule ligands.
  * Prosthetic groups.
  * Non-standard amino acids.
2. Remove them. In the case of cofactors and prosthetic groups, these have to be parametrised. In the case of modified amino acids, check if they are contained in Amber parameters. In the case of solvent, it disturbs the calculations, and absolutely has to be removed. In the case of hydrogens, these are handled internally, and external software may assign them wrongly, so remove them all.
3. Parametrise all small molecules, ligands and prosthetic groups. To that end, remove them and store them in a separate **.pdb** file. Then run the following command: `antechamber -i [INPUT_NAME].pdb -fi pdb -o [OUTPUT_NAME].mol2 -fo mol2 -c bcc -s 2`, replacing `[INPUT_NAME]` and `[OUTPUT_NAME]` with the respective names of input and output files. These can then be reinserted into the structure by the following steps:
  1. Run the following command: `parmchk -i [INPUT_NAME].mol2 -f mol2 -o [OUTPUT_NAME].frcmod`. This will generate a **.frcmod** file, which contains the parameters of your small molecule or prosthetic group.
  2. Merge all **.frcmod** files into one.
  3. This can then be run using the MAWS.py script using the `--hybrid` option to point to the merged **.frcmod** file.

For small molecules, they have to be parametrised by running the command `antechamber -i [INPUT_NAME].pdb -fi pdb -o [OUTPUT_NAME].mol2 -fo mol2 -c bcc -s 2`.

## Running MAWS
After having completed the above steps, open the command line in your simulation directory and type:
* `python MAWS.py [options] [Path to file] [File name]`
* Options are set to reasonable defaults and you can get a list and description of all options by typing: `python MAWS.py -h`
* That's it! By following this guide, **you** can start computing aptamers, too. 


