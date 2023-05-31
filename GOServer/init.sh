#!/bin/bash

if [[ ! -d ./temp/ ]]
then
    pip3 install gdown
    gdown 1NFQqeo4af6It3ieXOhe1EV5hxh7T3Alp
    unzip MongoDumpAll.zip
fi