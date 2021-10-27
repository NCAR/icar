# Git workflow for ICAR

> Note: This document has been conveniently adapted from the [Cookbook for working with Git and VIC](https://github.com/UW-Hydro/VIC/wiki/Git-Workflow)

The basic workflow described here follows a workflow originally outlined by [Vincent Driessen](http://nvie.com/posts/a-successful-git-branching-model/). The workflow is built around the Git version control system. A basic description of the branching strategy and release management used for ICAR is presented here. We use a central truth repository ([https://github.com/NCAR/icar](https://github.com/NCAR/icar)) that contains our main branches and the official release version of ICAR. All active development takes place on forks and clones of this repo.

As you will note when reading below, this workflow can be a very branchy workflow.  Particularly if you work from multiple computers, this can lead to trying to remember a lot of git commands.  A useful tool to help with all these commands is called [git-flow](https://github.com/petervanderdoes/gitflow-avh). This allows you to do simpler things like:

```shell
git flow feature start <some_new_feature>
<work work work>
git flow feature finish <some_new_feature>
```

## Main Branches

There are two main branches: **main** and **develop**. The first one is the official release, the second one is the one undergoing active development.

 1. **main** -- The main branch represents the official release of the code. This branch is updated by new releases from the develop/release branches and by hotfixes. If you are not directly interested in model development, but want to use ICAR for your own modeling project, then this is typically the branch that you want to use.

 2. **develop** -– The develop branch represents the bleeding edge of the code. We recommend that all new development begins by branching from the develop branch.

 Both of the main branches are published on the Github page and are controlled by the ICAR admin group. The repository is organized into several branches for different purposes.

## Feature Branches

In general, any new feature branch should be based on the develop branch. Feature branches are used by developers to make self-contained changes to the ICAR source code. Keeping the changes self-contained makes it much easier to merge new features into the main source code. For example, if you want to include new features consisting of a different radiation scheme and a new advection scheme, you'd typically want to do that on two separate feature branches (or perhaps more). This will keep pull requests (to have your changes included in the main code) tractable. We merge completed feature branches into the develop branch for inclusion in the next release. A developer may have many feature branches, which may or may not be intended for inclusion in the main source code. Because context-switching is relatively easy with Git, the development of features using this strategy provides a clean and simple method for switching between features and the main branches during the development process.

## Support Branches

Sometimes ICAR development is driven by projects that require very specific modifications to the code that would not be appropriate for inclusion in a major release. The use of support branches allows for the continued development of the trunk while "supporting" project-specific versions of the ICAR code. Instead of completely removing these changes, which may be useful to others, we put the project-specific version of ICAR in a support branch and continue developing. Support branches are essentially branches that are not expected to be merged back into the development branch.

## Admin Branches

Although anyone could create these branches, they are designed for the preparation of a release of the main branch or a hotfix that cannot wait to the next major release. These branches should therefore only be used by members of the admin group.

 1. **release** -– The release branch supports the preparation of a new release. It includes any changes needed for the next public release or minor bug fixes.

 2. **hotfix** -- The hotfix branch facilitates mid-release bug fixes of the main branch. The key point of the hotfix branch is that it does not incorporate any new features from the develop branch, rather it is a branch off the main that addresses a specific issue or set of issues. When the hotfix is applied, the development branch is updated to reflect the hotfix changes.

## Naming Conventions
 * Main branch – main
 * Develop branch – develop
 * Feature branch – feature/{feature_name}
 * Hotfix branch – hotfix/{hotfix_name}
 * Release branch – release/{release_name}
 * Support branch – support/ICAR.{base_release_number}.{feature_branch_name}
 * Release name – ICAR.{major.minor.patch}
 * Support release name - ICAR.{base_release_number}.{feature_branch_name}.{##}

## User Permissions
Using Github to host the central or truth repository of our models allows us to easily control contributor permissions. Currently we split permission levels into 3 levels, Owners, Model Admins, and Developers.

 1. Owners have full access to all repositories and have admin rights to the organization.

 2. Model Admins have full access to specific repositories. They may push, pull, or make administrative changes to those repositories associated with their model. However, they should generally not push to the truth repo directly. Instead, they should fork, clone, edit locally, update their fork and then issue a pull request. This pull request should preferably be reviewed by someone else before it is merged.

 3. Developers have read-only access (pull, clone, fork) to any of the publically listed repositories under the NCAR name. If a developer would like a feature branch merged into the main repository, a pull request must be submitted and a Model Admin may merge it in.

## Workflow examples

### New feature

You have developed a novel way to parameterize convection and would like to add this as a process alternative to ICAR. We'll assume that you already have an account on GitHub and that you have the requisite software (Fortran compiler) and libraries (NetCDF, FFTW) already installed.

The process would be as follows:

 * Navigate to the main [ICAR repo](https://github.com/NCAR/icar)

 * Fork the repo by clicking on the 'Fork' button in the upper right corner

 * Navigate to your fork

 * Clone the fork to your local machine

 * Add the main ICAR repo as the upstream remote, so you can easily merge changes that are made in the main ICAR repo into your own local repo


```shell
git remote add upstream ssh://git@github.com:NCAR/icar.git
```

 * Checkout the `develop` branch
```shell
git checkout develop
```

 * Create and checkout the `feature/crazy_convection` branch (or whatever the appropriate name would be). If you create this branch while you are on the `develop` branch, the new branch will be based on `develop` (you can also specify this explicitly to git).
```shell
git checkout -b feature/crazy_convection
```
 * Push this new branch to your remote on GitHub
```shell
git push
```

 * Now make as many changes as you need to, commit them to your local repo and push them to your remote on GitHub. This is just like any other work you would do using Git. Once everything is working and everything is sufficiently tested, you will be ready to share your code with others.


 * Before you do that, merge any changes that have been made in the develop branch in the main ICAR repo into the `feature/crazy_convection` branch of your local repo. Assuming you are already on the `feature/crazy_convection` branch:

```shell
git fetch upstream
git merge upstream/develop
```

 * Resolve any merge conflicts

 * Push your latest version to your remote on GitHub

 * Issue a pull request. You do that on GitHub. Make sure that you make the pull request with respect to the correct branches. On your end this should be the `feature/crazy_convection` branch and on the other end the `develop` branch.

 * You changes will be reviewed and merged or more likely there will be some back-and-forth with suggested changes and clarifications.
