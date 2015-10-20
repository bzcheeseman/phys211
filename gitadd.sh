#!/bin/sh

current_exp="$1"
target="$2"

if [ !"$1" | !"$2"] ; then
	echo "Usage: ./gitadd.sh <current_exp> <target> <commit message>"
fi

cd "$current_exp/$target"
pdflatex "$target.tex"

cd "$(id -un)/desktop/phys211"

if [ !"$3" ] ; then
	git commit -am "$(date +%Y%m%d-%H%M%S)"
else
	git commit -am "$3" 
fi

git push -u origin master