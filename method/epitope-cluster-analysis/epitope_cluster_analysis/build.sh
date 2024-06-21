#!/bin/bash

rm -rf .git

rm  ./build.sh

cd ..
tar -czvf IEDB_Cluster-2.0.tar.gz epitope_cluster_analysis
