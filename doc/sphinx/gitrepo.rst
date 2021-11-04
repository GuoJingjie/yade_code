.. _yade-gitrepo-label:

##############
Yade on GitLab
##############

************************************************
Fast checkout (read-only)
************************************************
 
Getting the source code without registering on GitLab can be done via a single command. It will not allow interactions with the remote repository, which you access the read-only way::

 git clone --recurse-submodules https://gitlab.com/yade-dev/trunk.git

************************
Branches on GitLab
************************

Most useful commands are listed in the sections below. For more details, see these git guides:

1. `ProGit <https://github.com/progit/progit2>`_ `online Book <https://git-scm.com/book/en/v2>`_,

2. `Guide on setting up git <https://help.github.com/articles/set-up-git>`_,

3. `Git "choose your own adventure" <https://sethrobertson.github.io/GitFixUm/fixup.html>`_,

4. `Guide on fixing the conflicts <https://medium.com/@porteneuve/fix-conflicts-only-once-with-git-rerere-7d116b2cec67>`_.

Setup
=====

1. Register on `gitlab.com <https://gitlab.com>`_

2. Add your `SSH key <https://gitlab.com/-/profile/keys>`_ to GitLab

3. Set your username and email through terminal

   ::
 
      git config --global user.name "Firstname Lastname"
      git config --global user.email "your_email@youremail.com"
 
   You can check these settings with ``git config --list``.


4. To fork the repository (optional), click the “Fork” button on the `gitlab page <https://gitlab.com/yade-dev/trunk>`_, and also fork the `YADE data files <https://gitlab.com/yade-dev/yade-data>`_.

   .. note:: By default gitlab will try and compile the forked repository, and it will fail if you don't have runners attached to your account. To avoid receiving failure notifications go to repository settings (bottom of left panel->general->permissions) to turn of pipelines. 

5. Set Up Your Local Repo through terminal. The argument ``--recurse-submodules`` is to make sure that ``./data`` directory is filled with the recent data from `yade-data <https://gitlab.com/yade-dev/yade-data>`_ (the path is `relative <https://gitlab.com/yade-dev/trunk/-/blob/master/.gitmodules#L3>`_ to your gitlab profile):

   ::

      git clone --recurse-submodules git@gitlab.com:username/trunk.git
   
   This creates a new folder, named trunk, that contains the whole code (make sure username is replaced by your GitLab name). If you already have a cloned yade repository with ``./data`` directory in it, then you can populate your existing repository using command::

      git submodule update --init --recursive

6. Configure remotes

   ::

      cd to/newly/created/folder
      git remote add upstream git@gitlab.com:yade-dev/trunk.git
      git remote update
    
   Now, your "trunk" folder is linked with two remote repositories both hosted on gitlab.com, the original trunk from yade-dev (called "upstream" after the last command) and the fork which resides in your personal account (called "origin" and always configured by default). Through appropriate commands explained below, you will be able to update your code to include changes commited by others, or to commit yourself changes that others can get.
   
   Holding a fork under personnal account is in fact not strictly necessary. It is recommended, however, and in what follows it is assumed that the above steps have been followed.

Older versions
========================

In case you want to work with, or compile, an older version of Yade which is not tagged, you can create your own (local) branch of the corresponding daily build. Look `here <https://answers.launchpad.net/yade/+question/235867>`_ for details.

Committing and updating
========================

Inspecting changes
------------------

After changing the source code in the local repository you may start by inspecting them with a few commands. For the "diff" command, it is convenient to copy from the output of "status" instead of typing the path to modified files. ::

   git status
   git diff path/to/modified/file.cpp

Pushing changes
---------------

Depending on the remote repository you want to push to, follow one of the methods below.

1. Push to yade-dev

   Merging changes into yade-dev's master branch cannot be done directly with a push, only by merge request (see below). It is possible however to push changes to a new branch of yade-dev repository for members of that group. It is `currently <https://gitlab.com/gitlab-org/gitlab-ce/issues/23902>`_ the only way to have merge requests tested by the gitlab CI pipeline before being effectively merged. To push to a new yade-dev/branch::

      git branch localBranch
      git checkout localBranch
      git add path/to/new/file.cpp  #Version a newly created file
      git commit path/to/new_or_modified/file.cpp -m 'Commit message'  #stage (register) change in the local repository
      git pull --rebase upstream master #get updated version of sources from yade-dev repo and apply your commits on the top of them
      git push upstream localBranch:newlyCreatedBranch #Push all commits to a new remote branch.

   The first two lines are optional, if ignored the commits will go the to the default branch, called "master".
   In the last command ``localBranch`` is the local branch name on which you were working (possibly ``master``) and ``newlyCreatedBranch`` will be the name of that branch on the remote. Please choose a descriptive name as much as you can (e.g. "fixBug457895").

.. note:: If you run into any problems with command ``git pull --rebase upstream master``, you :ref:`always can revert<yade-rebasing>` or even better `fix the conflicts <https://medium.com/@porteneuve/fix-conflicts-only-once-with-git-rerere-7d116b2cec67>`_.


2. Push to personnal repository

   After previous steps proceed to commit through terminal, "localBranch" should be replaced by a relevant name::

      git branch localBranch
      git checkout localBranch
      git add path/to/new/file.cpp  #Version a newly created file
      git commit path/to/new_or_modified/file.cpp -m 'Commit message'  #stage (register) change in the local repository
      git push  #Push all commits to the remote branch
  
   The changes will be pushed to your personal fork.
   

Updating
--------

You may want to get changes done by others to keep your local and remote repositories synced with the upstream::

 git pull --rebase upstream master #Pull new updates from the upstream to your branch. Eq. of "bzr update", updating the local branch from the upstream yade-dev/trunk/master
 git push  #Merge changes from upstream into your gitlab repo (origin)

If you have local uncommited changes this will return an error. A workaround to update while preserving them is to "stash"::

 git stash #backup and hide changes
 git pull --rebase upstream master
 git push
 git stash pop #restore backed up changes


.. _yade-auto-rebase:

Auto rebase
-----------

We promote "rebasing" to avoid confusing logs after each commit/pull/push cycle. It can be convenient to setup automatic rebase, so it does not have to be added everytime in the above commands::

    git config --global branch.autosetuprebase always

Now your file ``~/.gitconfig`` should include::

	  [branch]
	    autosetuprebase = always

Check also ``.git/config`` file in your local trunk folder (rebase = true)::

	  [remote "origin"]
	    url = git@gitlab.com:yade-dev/trunk.git
	    fetch = +refs/heads/*:refs/remotes/origin/*
	  [branch "master"]
	    remote = origin
	    merge = refs/heads/master
	    rebase = true


.. _yade-rebasing:

Pulling a rebased branch
------------------------

If someone else rebased on the gitlab server the branch on which you are working on locally, the command ``git pull`` may complain that the branches have diverged, and refuse to perform operation, in that case this command::

 git pull --rebase upstream branchName

Will match your local branch history with the one present on the gitlab server.

If you are afraid of messing up your local branch you can always make a copy of this branch with command::

 git branch backupCopyName

If you forgot to make that backup-copy and want to go back, then make a copy anyway and go back with this command::

 git reset --merge ORIG_HEAD

The ``ORIG_HEAD`` backs up the position of ``HEAD`` before a potentially dangerous operation (merge, rebase, etc.).

A tutorial on `fixing the conflicts <https://medium.com/@porteneuve/fix-conflicts-only-once-with-git-rerere-7d116b2cec67>`_ is a recommended read.

.. note:: If you are lost about how to fix your git problems try `a git choose your own adventure <https://sethrobertson.github.io/GitFixUm/fixup.html>`_.

**************
Merge requests
**************

Members of yade-dev
===================

If you have tested your changes and you are ready to merge them into yade-dev's master branch, you'll have to make a "merge request" (MR) from the gitlab.com interface (see the "+" button at the top of the repository webpage). Set source branch and target branch, from yade-dev/trunk/newlyCreatedBranch to yade-dev/trunk/master. The MR will trigger a `pipeline <https://gitlab.com/yade-dev/trunk/pipelines>`_ which includes compiling, running regression tests, and generating the documentation (the `newly built <https://yade-dev.gitlab.io/trunk>`_ documentation is accessible via settings->pages or by clicking on the "Browse" button in the "Job artifacts" (in the right pane) in the ``doc_18_04`` build from the pipeline; then navigating to path ``Artifacts/install/share/doc``).
If the full pipeline succeeds the merge request can be merged into the master branch.

.. note::
   In case of MR to yade-dev's master from another branch of yade-dev, the pipeline will use group runners attached to yade-dev (the group runners are kindly provided by `3SR <https://www.3sr-grenoble.fr/?lang=en>`_, `UMS Gricad <https://gricad.univ-grenoble-alpes.fr/>`_ and `Gdańsk University of Technology <https://pg.edu.pl/>`_).

New developers
==============

Welcome! At start it is very convenient to work on a local fork of YADE in your own gitlab profile. When you are confident that your changes are ready to be merged into official YADE release, please open a Merge Request (MR) in the following way:

	1. Make sure that your work is in a separate branch, not in the ``master`` branch. You can "copy" your branch into another branch with command ``git checkout -b myNewFeature``. Please make sure that the amount of changes as compared to the master branch is not large. In case of larger code improvements it is better to split it into several smaller merge requests. This way it will be faster for us to check it and merge.

	2. Push your branch to the repository on your gitlab profile with command such as::

		git push --set-upstream origin myNewFeature

	3. You should see something like::

		remote:
		remote: To create a merge request for myNewFeature, visit:
		remote:   https://gitlab.com/myProfileName/trunk/-/merge_requests/new?merge_request%5Bsource_branch%5D=myNewFeature
		remote:

	4. When you visit the link mentioned above, you will have to select "Change branches" and make sure that correct target branch is selected. Usually that will be ``yade-dev/trunk:master``, because this is the official YADE repository.

	5. Fill in the title and description then click "Create merge request" at the bottom of the page.

	6. After we review the merge request we can click on it to run in our Continuous Integration (CI) pipeline. This pipeline can't start automatically for security reasons. It will be merged after the pipeline checks pass.


Alternatively, create a patch from your commit via::

 git format-patch origin  #create patch file in current folder)

and send to the developers mailing list (yade-dev@lists.launchpad.net) as attachment. In either way, after reviewing your changes they will be added to the main trunk.

When the pull request has been reviewed and accepted, your changes are integrated in the main trunk. Everyone will get them via ``git fetch``.

**********************
Guidelines for pushing
**********************

These are general guidelines for pushing to ``yade-dev/trunk``.

1. Set autorebase globaly on the computer (only once see above), or at least on current local branch. Non-rebased pull requests will not be accepted on the upstream. This is to keep history linear, and avoid the merge commits.  

2. Inspect the diff to make sure you will not commit junk code (typically some "cout<<" left here and there), using in terminal:

   ::

    git diff file1
  
   Or using your preferred difftool, such as kdiff3:
  
   ::
  
    git difftool -t kdiff3 file1

   Or, alternatively, any GUI for git: gitg, git-cola... 

3. Commit selectively:

 ::

  git commit file1 file2 file3 -m "message" # is good
  git commit -a -m "message"                # is bad. It is the best way to commit things that should not be commited

4. Be sure to work with an up-to-date version launching:

 ::

  git pull --rebase upstream master

5. Make sure it compiles and that regression tests pass: try ``yade --test`` and ``yade --check``.


**Thanks a lot for your cooperation to Yade!**
