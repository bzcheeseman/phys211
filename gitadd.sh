#!/bin/sh

target="$1"
msg="$3"
if [ "$3" ] ; then
	git add "$3"/*
	git commit -m "$msg"
	git push -u origin master
fi
if [ !"$3" ] ; then
	git commit -m "$msg"
	git push -u origin master
fi