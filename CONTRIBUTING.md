# Contributing to DIscBIO

What follows is a workflow for contributing to DIscBIO. By following the suggestions below, you make it easier for us to review and incorporate your contribution.

If you are officially listed as a contributor to the package, you can contribute directly to the `dev` channel. If you are not, you should create a fork of the upstream repository before creating a local clone. Instructions for forking a project can be found [here](https://guides.github.com/activities/forking/).

## Setup a local git repository

Assuming you already have git installed in your local machine, open a terminal and issue a command to clone the remote (i.e. on GitHub) repository. If you are an official contributor to the package, the code below should work. Otherwise, you should write your fork address instead of `ocbe-uio`.

```
git clone https://github.com/ocbe-uio/DIscBIO.git
```

This will create a local copy of the DIscBIO repository in your computer.

## Create a separate, local branch for your contribution

Please avoid making work-in-progress changes directly to the `dev` channel. It is generally considered a good idea to make changes to a separate branch and only merging them with a permanent branch such as `dev` when you are reasonably sure your changes pass all the currently-implemented unit tests.

Let's start by creating and checking out a feature branch called `newFeature`:

```
git checkout -b newFeature
```

This branch exists only in your machine, not in the remote repository.

After checking out the new branch, you are ready to change files. After you are done editing files, commit the changes to your feature branch:

```
git add <changed files>
git commit -m "<Meaningful message about what has changed>"
```

Try to group your commits by subject. If your changes involve new documentation and new code, it is best to make one commit containing the new code and another with the documentation changes.

## Test your changes

Before merging your changes with the development channel, make sure your changes do not break the package by running the unit tests. Inside an interactive R session, issue the following command:

```
devtools::test()
```

If there are any test failures, please address them in a new commit. Once your changes pass these tests, perform a general check of the package to see if the code and documentation are also proper:

```
devtools::check()
```

If there any errors or warnings, please make a new commit with fixes and retest. Once the package checks without errors or warnings, you are ready for the next step.

## Merge your feature branch back to `dev`

The final step is to merge your changes into your local copy of the development channel and then push it back upstream (to the remote repository)

```
git checkout dev
git merge newFeature
git push origin
```

If you are working on a fork, you should instead push to your fork and open a pull request from it to the upstream repository.

## Style guide

Writting functional code is essential, but writing code which is functional and pleasing to the eye goes a long way into making your contribution understandable and thus more useful for the scientific community. It can be quite a controversial topic, since aesthetics are very much a subjective matter. As Hadley Wickham puts it, "good coding style is like correct punctuation: you can manage without it, butitsuremakesthingseasiertoread".

You are free to write code using the style you prefer, but we recommend reading and following [The tidyverse style guide](https://style.tidyverse.org/index.html), particularly [chapters 2](https://style.tidyverse.org/syntax.html) and [3](https://style.tidyverse.org/functions.html).
