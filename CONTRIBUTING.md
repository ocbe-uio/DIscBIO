# Contributing to DIscBIO

What follows is a workflow for contributing to DIscBIO. By following the suggestions below, you make it easier for us to review and incorporate your contribution.

These guidelines assume you are officially listed as a contributor to the package. If you are not, there is an extra step in the beginning: you should create a fork of the upstream repository before creating a local clone.

#TODO: add instructions for forking

## Setup a local git repository

Assuming you already have git installed in your local machine, open a terminal and issue

```
git clone https://github.com/ocbe-uio/DIscBIO.git
```

on a terminal window. This will create a local copy of the DIscBIO repository in your computer.

## Create a separate, local branch for your contribution

Please avoid making changes directly to the `dev` channel. It is often considered a good idea to make changes to a separate branch and only merging them with a long-lived branch such as `dev` when you are reasonably sure your changes pass all the currently-implemented unit tests.

Let's start by creating and checking out a feature branch called `newFeature`:

```
git checkout -b newFeature
```

Now you are ready to change files. After you are done making changes, commit them to your feature branch

```
git add <changed files>
git commit -m "<What was changed>"
```

Try to group your commits by subject. If your changes involve new documentation and new code, it is best to make one commit containing the new code and another with the documentation changes.

## Test your changes

Before merging back to the development channel, make sure your changes do not break the package by running the unit tests. Inside an interactive R session, issue the following command:

```
devtools::check()
```

If there any errors or warnings, please make a new commit with the changes and retest. Once the package checks without errors or warnings, you are ready for the next step.

## Merge your feature branch back to `dev`

The final step is to merge your changes into your local copy of the development channel and then push it back upstream (to the remote repository)

```
git checkout dev
git merge newFeature
git push origin
```

If you are working on a fork, you should instead push to your fork and open a pull request from it to the upstream repository.