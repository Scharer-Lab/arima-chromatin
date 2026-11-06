# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact


## Shared R script dependencies

The affected scripts load shared R code from local Scharer-Lab repository
checkouts rather than a machine-specific absolute path. Clone the following
repositories as sibling directories using the names shown:

- [pca](https://github.com/Scharer-Lab/pca)

Set `SCHARER_LAB_REPOS` to the absolute path of their parent directory before
running or sourcing the scripts. For example:

```sh
export SCHARER_LAB_REPOS="$HOME/scharer-lab"
```

The target files must be present in those checkouts. Access to private
repositories requires the appropriate GitHub permissions. This does not install
R packages, reference genomes, external tools, or project data.
