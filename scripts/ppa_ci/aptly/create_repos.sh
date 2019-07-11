#!/bin/bash

set -e
for i in bionic buster stretch
do
    aptly repo create -distribution=$i -component=main yadedaily-$i
done