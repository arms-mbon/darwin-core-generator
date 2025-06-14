# darwin-core-generator

This code is tuned for data produced by the PEMA bioinformatics pipeline for 18S, ITS, and COI, where those data are accessed via GitHub and are formulated in a specific way. 

The code main.py is the file you will use to run the code (see further instructions below) and that calls on the other python codes in the folder arms_toolbox. Please note that currently we are using the more up-to-date code in the folder arms_toolbox_datarelease001 (this was first created for the ARMS-MBON data release #1, and is also being used for data release #2) instead: the two will be merged at some point in the future.  

## Local usage and testing

Note: we only give the linux-style shell commands here. For users on Windows we recommend simply exploiting these exact same commands as they are supported in the "terminal" mode of the popular MS Virtual Code tool.

Still YMMV and practical execution on your platform might uncover untested / unwanted effects. So pls let us know if important changes and clarifications should be made.


### check dependencies

This tool requires python 3.12

You should be able to check with one of:

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

To run this code on a windows machine I did the following
* You need to clone this repository using e.g. Github Desktop. It is important that you have a complete copy of this entire repo and its folders and subfolders, else the code will not run
* Make sure you have python 3.12 installed. Because of the wonderful way computers work, if you update to that version (download and install it directly from a website) and expect python to have been updated in your working environment, you can forget it.
   Operating systems are designed by twisted minds, so you should get an expert to do it for you. We had to delete old versions, install the new version - remembering to tell it to update my PATH - and then reboot. 
   Always check by typing $python --version from whatever terminal you will run the code from
* Then I used GHD to clone this repo and asked to open it (from the GHD) using y default editor, being VS Code. 
* In VS code I then selected the GitBash terminal, not the default powershell 
* I want to run this in a virtual environment, so that it is isolated from everything else I do. There I need "venv", which fortunately I had, but if you don't then you can pipe install it from the terminal in VSC or any other terminal. So go to the darwin-core-generator folder in the gitbash terminal
* Now set up the virtual env, calling it "venv" (or you could call it "jane" but then you will forget later what "jane" is). 
   Type $ python -m venv venv 
   in your GitBash termin in VSC; first venv=the module, second=the name
   You will see a directory called "venv" has been created in this repo
* Now to run that virtual environment, type 
   $source venv/Scripts/activate
   Note: as long as you don't delete this entire directory/repo from your computer, this virtual envrionment will stay, so when you later exit it, go home, come back, you can start again in your VCS in this venv. 
* I needed to 
   $pip install -r requirements.txt
   to get the necessary libraries in place
* Now I run the code with 
   $ python main.py
* Exit 
* push back to GH when ready. Remember to set a .gitignore to ignore the "venv" directory 

Note that the code works on files that are in the GH repo data_workspace, but its outputs are send to this repo to the data/outputs folder.





