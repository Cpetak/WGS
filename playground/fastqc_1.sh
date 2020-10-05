#!/bin/bash

mkdir fastqc

# fastqc is a program we are calling that is already available, followed by input and output location
fastqc example.fastq -o fastqc/ # because I'm going to be in myresults already


