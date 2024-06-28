# darwin-core-generator

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