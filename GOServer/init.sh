#!/bin/bash

if [[ ! -d ./temp/ ]]
then
    echo "folder already exists"
else
    python3 ./download.py
fi
