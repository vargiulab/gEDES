#!/bin/bash
ls -lths *opt1-rescored*.pdb | sort -t "_" -nk 6 | awk '{if ($1 != "0") printf "mol addfile %s\n", $NF}' | sed '0,/addfile/s//new/' > load-opt1-rescored.tcl

