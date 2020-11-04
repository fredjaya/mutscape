#!/bin/bash

awk -F '\t' '/^\# PSC\t/;/^PSC/' $1 | \
	sed s/\#\ /''/ | sed -E s/\\[[0-9]+\\]//g \
	> $2
