#!/bin/bash
# script to clean files from a nautilus simulation
# version 1.1

function clean {
    rm *.out
    rm *.tmp
    rm ab/*
    rmdir ab
    rm struct/*
    rm *.percentage
    rm *.reaction
    rm *.pdf
    rmdir struct
    # To delete the stderr and stdout of a bash scheduler of the server. 
    # The last "." is very important, in order to avoid suppression of the submission script itself.
    rm *.sh.* 
}





echo "deleting files..."
clean
echo "done"


