#!/bin/bash

function timestamp { PASSED=$(echo $(date +%s) - $2 | bc); echo "* Process '"$1"' done ["$(date -d@$PASSED -u +%H:%M:%S)"]"; }

for var in '1-1' '1-2' '1-5' '1-6' '1-9'
do

start_time=$(date +%s);

sleep 1

timestamp sample-"$var" $start_time;
	
done

start_time=$(date +%s);

sleep 1

timestamp sample-1-3-4 $start_time;

start_time=$(date +%s);

sleep 1

timestamp sample-1-7-8 $start_time;

start_time=$(date +%s);

sleep 1

timestamp SEAL $start_time;
