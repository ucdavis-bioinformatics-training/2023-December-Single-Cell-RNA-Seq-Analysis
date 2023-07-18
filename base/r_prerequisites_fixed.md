# R and RStudio Prerequisites

The analyses in this course are done in R, using the RStudio environment. Some familiarity with R will help you make the most of your time. While we will provide R code that can be run even if you have never used R before, we find that participants who have spent some time exploring R before the course get the maximum value out of our workshops. There is not be enough time to provide full instruction in introductory R, and still get through all of the material we need to cover.

Beginners without any previous knowledge will be able to complete this course, and achieve a more thorough understanding of the techniques and analyses covered, but may not be able to conduct an experiment on their own without further work in R.

## Outline

1. What are R and RStudio?
2. Exploring the RStudio interface
3. R basics
4. Installing packages
5. Introduction to the tidyverse
6. Overview of single cell analysis

## What are R and RStudio?

[R](http://r-project.org/) is a language and environment for statistical computing and graphics. R is a popular choice for bioinformatics data analysis because it:

* is free and open source
* compiles and runs on Linux, Windows and MacOS
* has a large and active user community
* is easily extensible
* produces publication-quality plots, including mathematical symbols and formulae where needed
* provides thorough documentation

[RStudio](http://rstudio.org/) is a full featured integrated development environment (IDE) for R. It's features include:

* syntax highlighting
* code completion
* smart indentation
* workspace browser and data viewer
* embedded plots
* package management in the form of "projects"

RStudio and its team have contributed to many R packages, including the tidyverse family of packages, which we will be using in this workshop.

## Exploring the RStudio interface

**1\.** Getting started

Let's start RStudio

<img src="figures/RStudio_open.png" alt="RStudio_open" width="800px"/>

**2\.** Open a new RScript File

File -> New File -> R Script

<img src="figures/RStudio_newfile.png" alt="RStudio_newfile" width="800px"/>

Then save the new empty file as Intro2R.R

File -> Save as -> Intro2R.R

## R basics

**3\.** Basics of your environment

The R prompt is the '>' , when R is expecting more (command is not complete) you see a '+'

<img src="figures/Prompt.png" alt="Prompt" width="500px"/>

**4\.** Writing and running R commands

In the source editor (top left by default)
type

    getwd()

Then on the line Control + Enter (Linux/Windows), Command + Enter (Mac) to execute the line.


**5\.** The assignment operator ( <- ) vs equals ( = )

The assignment operator is used assign data to a variable

    x <- 1:10
    x

<div class="r_output">[1]  1  2  3  4  5  6  7  8  9 10
</div>

In this case, the equal sign works as well

    x = 1:10
    x

<div class="r_output">[1]  1  2  3  4  5  6  7  8  9 10
</div>

But you can, but should **NEVER EVER DO**

    1:10 -> x
    x

<div class="r_output">[1]  1  2  3  4  5  6  7  8  9 10
</div>

The two act the same in most cases. The difference in assignment operators is clearer when you use them to set an argument value in a function call. For example:

    median(x = 1:10)
    x

<div class="r_output">Error: object 'x' not found
</div>

In this case, x is declared within the scope of the function, so it does not exist in the user workspace.

    median(x <- 1:10)
    x

<div class="r_output">[1]  1  2  3  4  5  6  7  8  9 10
</div>

In this case, x is declared in the user workspace, so you can use it after the function call has been completed. There is a general preference among the R community for using <- for assignment (other than in function signatures)

## Installing packages

## Introduction to the tidyverse

## Overview of single cell analysis
