#!/bin/sh

target="$1"
msg="$2"

pdflatex $target.tex

git add *
git commit -m "$msg"
git push -u origin master