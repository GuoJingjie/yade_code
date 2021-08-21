#!/bin/bash

set -e
for i in stretch buster bullseye bookworm bionic xenial focal
do
    aptly repo create -distribution=$i -component=main yadedaily-$i
done
