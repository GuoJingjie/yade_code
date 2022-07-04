.. _tutorial-installation:

Installing Yade (for Windows and Mac users)
===================================

In preparation for the THM short-course, we ask the participants to have a Linux Debian distribution installed on their laptop prior to arrival. 

If you already have a debian distribution on your laptop, please follow the `installation instructions on our website <https://yade-dem.org/doc/installation.html#packages>`_. 

If you do not have a debian distribution on your laptop, you have three ways to get one:

Easiest way - Use our premade Virtual Machine (Windows)
-------------------------------
We created a full debian machine preloaded with Yade + Paraview + Kate. You can install this easily with the following steps:

#. Download and install `VirtualBox <https://www.virtualbox.org/wiki/Downloads>`_ for your OS.
#. Download this `yade_machine.ova file <https://u.pcloud.link/publink/show?code=XZWvTHVZPlgEiVh0iELW8ifMCGU8J0qivdHX>`_ (this step may take a few minutes, so please be patient. The file is 6 gb.)
#. Open VirtualBox and click "Tools>Import" 
#. Locate the 'yade_machine.ova' file that you downloaded, and click "Next"
#. Edit the system properties to suit your needs. Set the CPU count to half of your laptop CPUs and the RAM to half of your total laptop RAM. 
#. Click "Import" on bottom right.
#. Start the machine and it should bring you into the Ubuntu desktop where you can open a new terminal (Ctrl-alt-T) and type 

.. code-block:: bash

    yadedaily --version

To test that yade is already installed and ready to go. 

Login details (feel free to change these as soon as you are into your new VM):
user: yadeuser
password: yadeuser

Less easy way - Create your own Virtual Machine (MacOS)
----------------
This is if the direction above do not work for you. The end result is the same. 

#. Download and install `VirtualBox <https://www.virtualbox.org/wiki/Downloads>`_  for your OS
#. Download the `Ubuntu 20.04 <https://releases.ubuntu.com/20.04.4/?_ga=2.37188649.861267526.1652781821-22456328.1652781821>`_ image 
#. Open VirtualBox and select "Machine>New..."
#. Select "Type" as "Linux" and "Version" as "Ubuntu (64-bit)" (If you do not see Ubuntu 64-bit, contact me directly for assistance).
#. Select at least 4gb of ram (preferably 8gb), Select "Create a virtual hard disk now"
#. Name the machine and then click "Create" 
#. Choose 20-30gb of storage, leave the remaining options as default.
#. Click "create" 
#. Click "Start" and then find the Ubuntu 20.04 .iso that you downloaded from Step 2. 
#. Follow the installation instructions (this will take some time depending on your HDD speed)
#. Once Ubuntu is fully installed and you are inside the machine, go ahead and install yade by opening a terminal (Ctrl-alt-T to open a new terminal). 

.. code-block:: bash

    sudo bash -c 'echo "deb http://www.yade-dem.org/packages/ focal main" >> /etc/apt/sources.list.d/yadedaily.list'
    wget -O - http://www.yade-dem.org/packages/yadedev_pub.gpg | sudo tee /etc/apt/trusted.gpg.d/yadedaily.asc
    sudo apt-get update
    sudo apt-get install yadedaily

#. Download and install Paraview.

You are all set!

Hard way - Dual Boot (MacOS or Windows)
--------------------
`Find instructions here <https://opensource.com/article/18/5/dual-boot-linux>`_