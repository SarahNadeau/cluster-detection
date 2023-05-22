#!/bin/bash

# This script is to convert a nextstrain nexus file (e.g. downloaded from Auspice) into something treeio::read.nhx can handle.

grep "^  tree one = " nextstrain__timetree.nexus | sed "s/  tree one = //g" > tmp.nhx
sed "s/&/&&NHX:/g" tmp.nhx > tmp2.nhx

# Note: this does not handle unique/custom fields
sed "s/,region/:region/g" tmp2.nhx | \
sed "s/,country=/:country=/g" | \
sed "s/,division=/:division=/g" | \
sed "s/,location=/:location=/g" | \
sed "s/,num_date_CI={\([^}]*\)}//g" | \
sed "s/,div=/:div=/g" > tmp3.nhx

# Edge case: get rid of parentheses that mess up parsing
sed "s/Georgia (Asia)/Georgia Asia/g" tmp3.nhx > nextstrain__timetree.nhx

rm tmp.nhx 
rm tmp2.nhx
rm tmp3.nhx
