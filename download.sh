#!/bin/bash

# downloads 1.7Gb dataset files
wget https://hmgubox.helmholtz-muenchen.de/f/c998398ccff04dd39237/?dl=1 -O 2018_Angelidis.tar
mkdir data; cd data
tar -xvf ../2018_Angelidis.tar
