#!/bin/bash

for file in work/*/*/
  do echo $file && \ 
  less $file/.command.log_usage 
done
