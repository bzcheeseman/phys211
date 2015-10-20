#!/bin/sh

target="$1"
current_exp="$2"
msg="$3"

find "$current_exp" "$target" | pdflatex

git add *
git commit -m "$msg"
git push -u origin master