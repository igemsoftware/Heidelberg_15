# JAWS manual
## What's JAWS?
Are you using ribozymes or DNAzymes? Or are you perhaps thinking about how you are going to switch that one gene? Then you have maybe even thought about selecting an aptazyme, which switches in response to a ligand. This used to take quite a lot of time, but with our software JAWS, you have your riboswitch or aptazyme in minutes, all on your very own laptop. You just need a DNAzyme or ribozyme and an aptamer (which you can get from our other software MAWS), and JAWS creates the aptazyme for you.

## Setup
In this section we will show you, how to set up everything needed to run JAWS.
1. Install your favourite distribution of Linux. We developed JAWS on Arch, but everything else should work fine.
2. Download and install ViennaRNA with its Python extension.

## Running JAWS
To run JAWS, open the command line and change directory to the one containing the JAWS.py script. Type the following command:
```
python JAWS.py [-h] [-p P_THRESHOLD] [-a ACTIVE_TOLERANCE]
[-i INACTIVE_TOLERANCE] [-e ENERGY_THRESHOLD]
[-m MAX_ENERGY_DIFFERENCE] [-s STEM_LENGTH] [-f SHIFT]
[-r {dna_matthews1999.par,dna_matthews2004.par,rna_turner1999.par,rna_turner2004.par,rna_andronescu2007.par}]
[-v VIENNA_PATH]
5PRIME aptamer 3PRIME
```

JAWS accepts the following parameters:
* `-h, --help` Print help and exit
* `-p P_THRESHOLD, --p-threshold P_THRESHOLD` Minimal probability of base pairing at which a base-pair is considered to be present (defaults to 10<sup>7</sup>)
* `-a ACTIVE_TOLERANCE, --active-tolerance ACTIVE_TOLERANCE` Number of nucleotides not conforming to the active structure to be tolerated (defaults to 10)
* `-i INACTIVE_TOLERANCE, --inactive-tolerance INACTIVE_TOLERANCE` Number of nucleotides not conforming to the inactive structure to be tolerated (defaults to 2)
* `-e ENERGY_THRESHOLD, --energy-threshold ENERGY_THRESHOLD` Maximal energy of structures allowed (defaults to -12)
* `-m MAX_ENERGY_DIFFERENCE, --max-energy-difference MAX_ENERGY_DIFFERENCE` Maximal energy difference between active and inactive state (defaults to 9)
* `-s STEM_LENGTH, --stem-length STEM_LENGTH` Length of the stem connecting the aptamer to the DNAzyme or ribozyme (defaults to 10)
* `-f SHIFT, --shift SHIFT` Number of nucleotides to be displaced upon binding of the ligand (defaults to 7)
* `-r PARAMS, --params PARAMS` ViennaRNA parameter set. Can be one of `dna_matthews1999.par,dna_matthews2004.par,rna_turner1999.par,rna_turner2004.par,rna_andronescu2007.par` to use one of the parameter sets included with ViennaRNA or a path to a custom parameter set (defaults to `dna_matthews2004.par`)
* `-v VIENNA_PATH, --vienna-path VIENNA_PATH` Directory containing ViennaRNA parameter files. Required only if using a parameter set included with ViennaRNA (defaults to `/usr/share/ViennaRNA`)
* `5PRIME` (positional parameter) 5' sequence of the DNAzyme/ribozyme
* `aptamer` (positional parameter) aptamer sequence
* `3PRIME` (positional parameter) 3' sequence of the DNAzyme/ribozyme

This provides the following output: Sequences fitting the desired switching behaviour, together with energies of active and inactive state, and the number of successful base pairs in the active and inactive state. Furthermore, it tells the user which sequences are good, and which are slightly subpar, to allow for an even easier user experience.
