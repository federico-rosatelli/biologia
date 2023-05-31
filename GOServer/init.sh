#!/bin/bash

if [[ ! -d ./temp/ ]]
then
    pipx install gdown
    pipx ensurepath
    gdown 1NFQqeo4af6It3ieXOhe1EV5hxh7T3Alp
    unzip MongoDumpAll.zip
fi