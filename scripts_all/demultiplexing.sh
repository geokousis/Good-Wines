#!/bin/bash

ls -1 params-pool*.txt > tmp_pools.txt
while read -r pool_file; do
    echo "$pool_file"
    ipyrad -f -p "$pool_file" -s 12 -c 25
done < tmp_pools.txt
rm tmp_pools.txt
