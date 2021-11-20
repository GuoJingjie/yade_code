# 2019 © Remi Monthiller <remi.monthiller@gmail.com>
# Bedload Simulation

This repository was created for the study of the influence of shape in turbulent bedload transport.
This document aims to give the informations necessary to use this framework.

## Getting started

### Prerequisites

Be sure to use a linux operating system. At least Ubuntu and Debian distributions should work. 
You will need a slightly modified version of YADE DEM to use the code and run the examples.
To get it clone or download the project at [https://gitlab.com/remi.monthiller/trunk](https://gitlab.com/remi.monthiller/trunk).
Move to the branch ***hydro_force_engine*** and compile the code.

In order to compile the code and install all prerequisites, you may want to see the instructions 
given on the Yade dem website before [https://yade$-dem.org/doc/](https://yade-dem.org/doc/).
Be sure to install all prerequisites mentionned on the website or the code won't compile.
Take a look to the ***Installation*** section and the ***Compilation*** subsection. The process of compilation is the same.

When all the prerequisites are installed, you could run the following commands for example, to launch the compilation :

```
git https://gitlab.com/remi.monthiller/trunk.git
cd trunk
git checkout hydro_force_engine
cd ..
mkdir build install
cd build 
cmake -DCMAKE_INSTALL_PREFIX=../install ../trunk
make
make install
```

The executable should be in the folder ***../install/bin***.
For convenience, you may want to rename the executable to ***yade*** and add the install folder to your path :

```
cd ../install/bin
mv yade-nnn.git-nnn yade 
PATH=$PATH:$PWD
```

You may want to add this PATH modification in your ***.bashrc***.

### Get the framework

Clone or download the project and unzip it. 
Then move to its location with a terminal to get started. For example :

```
git clone https://github.com/LeDernier/bedload_simulation.git
cd bedload_simulation
```

### Repository Description

The following executables may only be functionals the yade executable has been renamed to yade and its folder added to your ***PATH*** in your ***.bashrc***.
If not, you may want to edit the following executables so that they work in your case.

#### Executables

* ***add_extension***
May be useful to run computations on a cluster.
Adds an extension to folders.
Use :
```
./add_extension extension folders
```
Example :
```
$ ls
folder1/ folder2/ folder3/
$ ./add_extension dir folder\*
$ ls
folder1.dir/ folder2.dir/ folder3.dir/
```

* ***clean_batch***
By default, a simulation starts from the last timestep saved. 
This executable is useful to get rid of all previously saved data in order to restart completely a simulation.
Deletes all the simulation data of the cases given in argument. Keeps the parameters unchanged and the exctracted data.
Use :
```
./clean_batch case_folders
```
Example :
```
$ ls case-\*/data
case-1/data:
0.0000.yade 1.0000.yade 2.0000.yade
case-2/data:
0.0000.yade 1.0000.yade 2.0000.yade
$ ./clean_batch case-\*
$ ls case-\*/data
case-1/data:
case-2/data:
```

* ***create_batch***
May be usefull to create a batch of simulation without manually change parameters of the simulation.
Copy a case folder for each parameter value given in argument and changes the parameter files accordingly.
Use :
```
./create_batch case parameter values
```
Example :
```
$ ls
case/
$ ./create_batch case A 1.0 2.0 3.0
$ ls
case/ case_A_1.0/ case_A_2.0/ case_A_3.0/
$ diff case_A_1.0/params.py case_A_2.0/params.py
< 	A = 1.0
---
> 	A = 2.0
```

* ***extract_batch***
Before being post processed, after a simulation, data should be extracted using this executable.
The parameters used for the extraction may be edited in the file ***params_extract.py***.
Extracts the data of the simulations given in argument and stores them in pkl or json format.
Use :
```
./extract_batch cases
```
Example :
```
$ ls
case-1/ case-2/
$ ls case-1
clean data/ framework.py params.py run
$ ./extract_batch case-\*
$ ls case-1/
clean data/ data.pkl framework.py params.py run
```

* ***post_proc_batch***
Executable usefull to post process data.
The parameters used for the post processing may be edited in the file ***params_post_proc.py***.
Reads the data in pkl format in the cases given arguments and use them to plot graphics.
Use :
```
./post_proc_batch cases
```
Example :
```
$ ls
case-1/ case-2/ post_proc/
$ ls post_proc/

$ ./post_proc_batch case-\*
$ ls post_proc/
case_phi.pdf case_vx.pdf case_qsx.pdf
```

* ***remove_extension***
May be useful after the use of ***add_extension***.
Removes an extension to folders.
Use :
```
./remove_extension folders
```
Example :
```
$ ls
folder1.dir/ folder2.dir/ folder3.dir/
$ ./remove_extension folder\*
$ ls
folder1/ folder2/ folder3/
```

* ***run_batch***
Runs a batch of simulations in the background.
Use :
```
./run_batch cases
```
Example :
```
$ ls
case-1/ case-2/
$ ls case-1
clean data/ framework.py params.py run
$ ./run_batch case-\*
$ ls case-1/
clean data/ nohup.out framework.py params.py run
```
The file ***nohup.out*** can give information about the simulation.
You may want to use the ***top*** program to check the ressources used.


* ***simple_post_proc.py***
A simple, commented post processing executable.
Use :
```
./simple_post_proc.py
```

#### Folders

* ***common***
Contains all the python files common to all simulation cases.
That is where to look if you want to extend the framework. 
Contains the following :
    * ***caseBatch.py***
    File executed when running ***run_batch***.
    * ***case.py***
    File executed when running ***run***.
    * ***defaultParams.py***
    Default simulation parameters.
    * ***extractData.py***
    File executed when running ***extract_batch***.
    * ***frameworkCreationUtils.py***
    File containing the utils needed for the framework creation.
    * ***measures.py***
    File containing the utils to "measure" physical quantities. Used in the post processing.
    * ***postProc.py***
    File executed when running ***post_proc_batch***
    * ***simulationDefinition.py***
    Used by ***caseBatch.py*** and ***case.py***, contains the utils that create the Yade engines from simulations parameters.
    * ***simulationPyRunners.py***
    Used by ***caseBatch.py*** and ***case.py***, contains the utils that create the python runners from simulations parameters.

* ***case-\****
Example of simulation cases. Contains the following :
    * ***clean***
    Executable to delete simulation data.
    Use : `./clean`
    * ***data/***
    Folder that contains the simulation data as *.yade files.
    * ***framework.py***
    Contains the parameters linked to the framework creation of the simulation.
    * ***params.py***
    Contains the parameters of the simulation.
    * ***run***
    Executable to run the simulation.
    Use : `./run`
* ***data-\****
Folders containing data that may be useful for post processing.

#### Others

* ***pp_\*.py***
Example of post processing parameter files.
* ***e_\*.py***
Example of extract parameter files.

### Running a first test

Move to an example case and run the case. If all goes well a control panel of the simulation should pop up.
For example case-wet:
```
cd case-wet
./run
```

### Running validations

The code is validated in comparison with experimental data. To check the validity of the code, you can run the following :
```
./run_batch case-exp*
```
Wait until the simulations end and then execute the following :
```
cp pp_exp6 params_post_proc.py
./post_proc_batch case-exp6
cp pp_exp14 params_post_proc.py
./post_proc_batch case-exp14
cp pp_exp20 params_post_proc.py
./post_proc_batch case-exp20
```

## Using the framework

Copy an example folder and edit the parameters, then run your simulation.
For example :

```
cp -r case-wet case-test
cd case-test
gedit params.py
./run
```

Wait until at least two seconds of simulations elpsed, and run the data extraction and the post processing :
```
cd ..
./extract_batch case-test
./post_proc_batch case-test
```

## Further improvements

- [ ] Cylinders case (especially the fluid coupling)
- [ ] Orientations post processing is not working well
- [ ] Simplification of the post processing
