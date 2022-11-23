#!/bin/bash

gunzip $1.gz
seqtk sample -s100 $1 100000 > test_R1.fq
gzip $1

gunzip $2.gz
seqtk sample -s100 $2 100000 > test_R2.fq
gzip $2
