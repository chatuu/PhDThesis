#!/bin/bash

printf "\n Creating New Direcory with the name: %s" $1
mkdir $1

cp *.h *.C run.sh startGridJob.sh ./$1

