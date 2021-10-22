# ICAR and Git

> Note: This document has been conveniently adapted from the [SUMMA and Git workflow](https://github.com/NCAR/summa/blob/master/docs/howto/summa_and_git_howto.md)

## Introduction

Git is a version control tool that we use to...

 * track changes over time in the ICAR source code,
 * facilitate parallel development between multiple model developers and research projects,
 * encourage peer code review and bug tracking, and
 * to distribute official and bleeding edge releases of the model source code.

## ICAR for model users (non developers)
In general, if you plan to apply the model rather than work directly on the source code, you may just download the model from the [ICAR GitHub page](https://github.com/NCAR/icar). Additional archived tags and releases (including all previous ICAR versions) are available [here](https://github.com/NCAR/icar/tags). You can simply compile the code and start working. If you want to work with the code and contribute to ICAR development, read on.

## ICAR for developers
If you plan on contributing to model development or would like a systematic way to incorporate updates to the ICAR source code, we encourage you to use Git. The following sections are designed to get you started using Git and working with the ICAR source code repository.

## Quickstart Cheatsheet
More details are provided below, but if you just want a reminder of the key commands to set up your repository (replace `https://` with `ssh://git@` if appropriate) here is the basic setup:

```shell
git clone https://github.com/<username>/icar.git
cd icar
git remote add upstream https://github.com/NCAR/icar.git
```

Then from the develop or main branch, periodically run:
```shell
git fetch upstream
git merge upstream/<branch_name>
```

or:
```shell
git pull upstream <branch_name>
```

### Git resources
If you are not familiar with Git yet, we encourage you to spend a few minutes getting acquainted with the system before you starting working with the ICAR source code and Git. It's not difficult to use and a few minutes of learning about Git will go along way in helping you manage your code development.

There are a number of good Git learning resources that will provide a basic introduction to the version control system.
* http://git-scm.com/about
* https://help.github.com/

### Getting the code
The basics steps to get the ICAR source code repository are as follows. This is basically a ICAR specific rendition of https://help.github.com/articles/fork-a-repo

#### Step 0: A few definitions

Using Git involves a few different copies of the archive (there is more than one way to use Git -- this is how we propose you use it):

1. **Truth repo:** this is the main copy of the archive, maintained on the GitHub server at https://github.com/NCAR/icar. You will not be able to edit this or push code to it directly. We don't either, we work on our own copy of this repo and only integrate the changes into the truth repo when they have been tested and reviewed.

2. **Your fork of the truth repo:** this is your version of the archive, maintained on the GitHub server by you. This is generally **not** where you edit code; it is more of a transfer point between your clone (see below) and the truth repo. Git keeps track of the differences between your fork and the truth repo. To move code changes from your fork to the truth repo, you must make a "pull request" to the ICAR truth repo. If it's decided that your code has been sufficiently tested and reviewed and that it is a useful addition to ICAR, then the administrators of the truth repo will do the actual pull (don't worry if some of the terminology does not yet make sense at this point).

3. **Your clone(s) of your fork of the truth repo:** a clone is a copy of your fork, stored on your local machine, where you can edit and test the code. You can have one or more clones. Each clone is a snapshot of a particular version of the code in your fork (e.g., you select which branches, and which versions of the code on those branches, to copy to your clone). You can change which version of the code each of your clones points to at any time. Your clone is where you edit files and commit changes. You can then push code changes from your clone back up to your fork.

Note that you do not have to use GitHub to use Git. You can clone the Truth Repo directly and then store your clone on a Git server other than GitHub. This means that you can bypass the fork step if you know what you are doing. However, GitHub does offer some nice features (and it is free for open source projects), so in the following we'll assume you'll be using GitHub.

One more note: In most of the following we are assuming that you will do much of your work at a Unix-like prompt (or the Mac OS X Terminal). If that is not the case, the general workflow is the same, but details may differ.

#### Step 1: Fork the repo

If you do not have a github account, create one. It's free for public repos (like the ICAR one). Once you have an account, go to the [NCAR/icar](https://github.com/NCAR/icar) page and click "Fork" in the upper right hand corner.

#### Step 2: Clone your fork

You've successfully forked the ICAR repository, but so far it only exists on GitHub. From this point on, anything you do to this fork will not affect the Truth Repo, so you can even delete the fork or modify it without worrying that you will make affect the main ICAR repo. To be able to work on the project, you will need to clone it to your local machine.

Run the following code on your local machine:

```shell
git clone https://github.com/<username>/icar.git
```

or (depending on your setup)

```shell
git clone ssh://git@github.com/<username>/icar.git
```

where `<username>` is your GitHub username.
This clones your fork of the repository into the current directory. The clone will live in a directory called "icar".

#### Step 3: Configure remotes

When a repository is cloned, it has a default remote called `origin` that points to **your** fork on GitHub, **not** the original repository it was forked from. To keep track of the original repository, you need to add another remote. You can name this anything you want, but the name `upstream` is descriptive and an informal convention.

Change to the newly cloned "icar" directory:

```shell
cd icar
```

If you type

```shell
git status
```

you should get some informational message, such as

```shell
On branch main
Your branch is up-to-date with 'origin/main'.

nothing to commit, working directory clean
```

If instead, you get the message

```shell
fatal: Not a git repository (or any of the parent directories): .git
```

then you are not in a Git repository (note that the top level directory in a Git repo will have a `.git` subdirectory). Make sure you are in the right place. If you are then proceed.

Assign the truth repo to a remote tracking branch called `upstream`, which will allow you easily pull in changes that are made in the truth repo and keep your clone in-sync with that version. Once you get more familiar with Git, you will be able to control which updates to include:

```shell
git remote add upstream https://github.com/NCAR/icar.git
```

or (depending on your setup):

```shell
git remote add upstream ssh://git@github.com/NCAR/icar.git
```

#### Step 4. Sync your clone with the truth repo

##### 4.a. Fetch information from the truth repo

Before starting to edit the code, pull in any new changes to the truth repo that have been made by other people since you first created the clone. You may want to do this periodically.

```shell
git fetch upstream
```

If you have already made changes to the code, this command by itself will not overwrite your files. For updates from the truth repo to show up in your files, you must do a **merge**.

##### 4.b. Merge changes

Determine which branches you will need to work with. At the very least, this will include the main branch. If you are working on a hotfix or a feature branch that already exists, you will need this branch as well; the ICAR administrator has likely given you the name of the appropriate branch to use. Alternatively, you may want to create a new branch (e.g., if you are the first person to work on a new feature or bug fix). For more information about the branches in the ICAR archive, see the [Git workflow for ICAR](https://github.com/NCAR/icar/blob/main/docs/howto/icar_git_workflow.md).

For each branch, merge any changes from the truth repo into your local version:

```shell
git checkout <branchname>

git merge upstream/<branchname>
```

where `<branchname>` is the name of the branch you want to update. Note that merging can be intimidating at first, because you are likely to get some messages about merge conflicts. Don't despair, there are many good internet resources that will explain this in detail.

### Working with the code

#### Making changes

##### 1. Select a branch

Change your active branch to the desired branch:

```shell
git checkout <branchname>
```

where `<branchname>` is the name of the branch. If you are creating a new branch, you can create and change to the branch in one step:

```shell
git checkout -b <new-branchname>
```

where `<new-branchname>` is the name of your new branch.

##### 2. Make changes

You can edit the code using any editor or development environment you prefer. You can also create new files, and move, rename, or delete existing files. You will not be able to push these changes to your fork until you **commit** them. If you are moving or renaming file, it is best to use `git mv` (so Git can track the move). If you are permanently deleting a file, use `git rm`.

It is a good idea to **compile and test** your changes on your local machine before you commit them. This avoids extra commits to fix typos, etc., and makes sure that you are not left with half-broken code when you come back to the project (or when you are collaborating with someone else).

At any point during the process of changing the code, you can pull in any changes that other people have made via the fetch/merge procedure described above. At any point you can also discard the changes and revert to the previous version.

#### Committing changes

Before committing your changes, remove any extraneous files that have been created during compiling and testing (e.g., *.o files). An easy way to do this is to type `make clean` in icar's `code` subdirectory.

##### 1. Register your changes for commit (staging)

To register the changes to (or creation of) a specific file, there are two steps that need to be completed in Git. First, changes need to be moved to a staging area (staged) and then they all files in the changing area are committed. You cannot commit a change without staging first and staging alone does not commit your changes. To move the committed changes to a repo that is different from the clone you are working on (for example, your fork), you will need to do one more step and push the changes to the remote repo:

```shell
git add <filename(s)>
```

where `<filename(s)>` are the names of the files with changes you want to commit. You can use wildcards, but be careful, because this often results in including files that you did not mean to include.

To register moving or renaming any files:

```shell
git mv <old> <new>
```

where `<old>` and `<new>` are the old and new filenames.

To register the deletion of a file:

```shell
git rm <filename>
```

After you have added (staged) all the files that you want to include in a commit, check that what you have is all you want:

```shell
git status
```

This will tell you which files have changed, but have not been staged and also tells you how to unstage files. Note that you do not need to include all the changes you have made in one single commit. It is often most useful if you keep your commits small and very specific, because it makes it easier to undo things if you later decide that a previous code change was not nearly as useful as you once thought.

##### 2. Commit the changes

Now that you have added (or staged) the changes you want to include, you have to actually commit them to your local repo:

```shell
git commit
```

This will bring up a commit log in your default editor. The list of files whose changes will be committed (i.e. were staged via `git add`) is shown in the header at the top of the file. If you disagree with this list, exit the editor and add and unstage files as necessary to correct the list. Then try "git commit" again. If satisfied with the list of changed files, add a description of the set of changes (including a brief description of the problem that motivated the changes). Save and exit.

Note: writing good commit messages is essential for making sense of your code changes later. Or perhaps more importantly, for someone else to make sense of your code changes later. Generally, in Git it is recommended to describe the nature of the change in a brief one-liner, then leave a blank line, following by a more detailed description of the changes as necessary. The reason is that some of the tools that let you browse past changes (e.g. `git log`), use this one-liner to summarize code changes.

##### Pushing commits to your fork

After committing your changes, you should push them to your fork (which has the alias `origin`) stored on GitHub. If you don't remember what all your remotes are, simply type (`git remote -v`):

```shell
git push origin <branchname>
```

where `<branchname>` is the name of the branch where you made the commits.


##### Making a pull request

To make your changes visible to other users/developers, your changes must be incorporated into the truth repo. To do this, you must create a pull request on the GitHub server.

**NOTE:** We ask that you perform at least some basic tests on your code before you issue a pull request. Make sure the code compiles and runs for at least the test cases you have been working with. If it is a bug fix, make sure that it actually fixes the bug. If possible, try to make sure that it doesn't create a new bug. We are working on generating some standard tests that everyone can download and run for this purpose; until then, please test the code using your own input files.

The ICAR administrator and other developers will review your pull request and decide if/how they want to incorporate your changes into the code. They are likely to suggest some changes (code style, content, etc.).

##### Useful git configuration
If you frequently make (and delete) remote branches and work from multiple computers, it might be useful to configure git to automatically prune remote branches when they have been deleted by running:

```shell
git config --local fetch.prune true
```

This prevents you from deleting a remote branch from computer A, then accidentally creating it again from computer B.

### Git workflow

For us to leverage Git to its full potential, we have implemented a Git-oriented workflow. This requires developers to adhere to a few rules regarding branch names and merge requests. A full description of the workflow we use can be found [here](https://github.com/NCAR/icar/blob/main/docs/howto/icar_git_workflow.md).  A useful tool to help use this workflow is called [git-flow](https://github.com/petervanderdoes/gitflow-avh).
