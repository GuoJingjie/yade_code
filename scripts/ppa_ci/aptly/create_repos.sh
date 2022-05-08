#!/bin/bash

set -e
for i in bionic buster stretch xenial bullseye focal jammy bookworm
do
    aptly repo create -distribution=$i -component=main yadedaily-$i
done
