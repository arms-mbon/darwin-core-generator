# darwin-core-generator

## Linux

### check dependencies

This tool requires python 3.12. You should be able to check with one of:

```sh
$ python --version
$ python3 --version
$ python3.12 --version
```

Tip: Check https://cloudbytes.dev/snippets/upgrade-python-to-latest-version-on-ubuntu-linux for how to (safely) get newer versions of python on ubuntu.  Mind that the ubuntu release you run comes with and laregly requires / expects a very specific version of python. While you can add newer versions alongside, one should keep the default / core version around, and only use other ones in virtual environment contexts.

### create a virtual env

While optional, we strongly recommend running python projects in virtual environments. This allows to avoid version conflicts between dependencies of possibly unrelated python projects. 

Below is how to quickly set that up with `virtualenv`

Note: before executing any of the below statements, do make sure the current active working directory of your shell is pointing to this project folder.

```sh
$ virtualenv --version            # check if you have the tool installed
$ pip install virtualenv          # install the tool for venv management
$ virtualenv .venv -p python3.12  # create the local venv in .venv and link it to the python of choice
``` 

After this one can activate / and deactivate the venv in the shell through:

```sh
$ source .venv/bin/activate    # activate the venv
(venv) $                       # the prompt of the shell changes to indicate a venv is active

...

(venv) $ deactivate            # deactivate the env, return to normal shell
$                              # the pro,mpt changes back to normal
```

Below further shell commands stated on a line starting with `(venv) $` are expected to be executed in an active virtual environment.  If you chose not to use these, you can simply ignore.

```sh
(venv) $ python --version      # e.g check the linked / associated python interpreter
Python 3.12.4                  # should be <4.0,>=3.12
```

### install dependencies

```sh
(venv) $ pip install -r requirements.txt
```

### run the main

```sh
(venv) $ python main.py
```


## Windows

For users on Windows we recommend simply exploiting these exact same commands as above and as they are supported in the "terminal" mode of the popular MS Virtual Code tool.

However, we provide slightly fuller instructions for those for whom the above set of commands are too advanced.

To run this code on a windows machine:
 * You need to clone this repository using e.g. Github Desktop (GHD). It is important that you have a complete copy of this entire repo and its folders and subfolders, else the code will not run
 * Make sure you have python 3.12 installed. Because of the wonderful way computers work, if you update to that version (download and install it directly from a website) and expect python to have been updated in your working environment, you can forget it. Operating systems are designed by twisted minds, so you should get an expert to do it for you. We had to delete old versions, install the new version - remembering to tell it to update my PATH - and then reboot. Always check by typing $python --version from whatever terminal you will run the code from 
 * Then I used GHD to clone this repo and I asked to open it (from the GHD) using the default editor, being VS Code for me.
 * In VS code I then selected the GitBash terminal, not the default powershell 
 * I want to run this in a virtual environment, so that it is isolated from everything else I do. There I need "venv", which fortunately I had, but if you don't then you can "pip install" it from the terminal in VSC or any other terminal
 * Now set up the virtual env, calling it "venv" (or you could call it "jane") with the command ``` python -m venv venv``` in your GitBash termin in VSC (the first venv=the module, second=the name). You will see a directory called "venv" has been created in this repo
 * Now to run that virtual environment, type ``` source venv/Scripts/activate```. Note: as long as you don't delete this entire directory/repo from your computer, this virtual envrionment will stay, so when you later exit it, go home, come back, you can start again in this venv. 
 * The first time running the code, I need to type ``` pip install requirements.txt``` to get the necessary libraries in place
 * Then I run the code with ``` python main.py```
 * When pushing back to GH (if that is what you want to do), remember to set the gitignore for the "venv" folder

Note: This code is heavily tuned to our use case, being taking a certain set of PEMA outputs, combined with event logsheet metadata 
and turning those into DwC. There are many ifs and buts in here that are specific to these data

## Note about the ids in the various files
Extended_final_table (col 1)
 * For 18S and ITS, it is Otu(a number)
 * For COI it is ASV_##:(long number)

Fasta files (these are read in via main.py)
 * For 18S and ITS it is the same Otu(number)
 * for COI it is (a long number)_(a shorter number), i.e. the second part of that from the Extended_final_table with a _(number) appended

Tax_assigments
 * For 18S and ITS there is no such file
 * For COI it is the same as in the fasta file

## Other notes
 * The emof schema can be optimised to the gene type.
 * This code is designed to be run over several input files whoes name differ by gene type (18S and ITS) and date (e.g. July2020), as specified in the main.py. In there are also the URLs of the fasta files for each input file, because these need to be downloaded from the MDA
 * The output is added to data/outputs
 * Checks for duplicates and missing information are performed - for duplicates, a warning is logged and the first only is taken
 * A check on ENA to find the sample accession numbers for the run accession numbers that we have in the omics input is done
 * ITS, 18S, and COI differ in that : the input files come from different repo folders; COI has more in the emof and takes that from a PEMA file that only exists for this one
 * A check is run on ENA to get the sample accession number for the given run accession number

