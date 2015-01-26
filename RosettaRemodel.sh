#!/bin/sh

remodel.macosclangrelease -database rosetta_2014.22.56873_bundle/main/database -s 1bxw.pdb -remodel:blueprint blueprint_OR4.txt -chain A -linmem_ig 10 -ex1 -ex2 -ignore_zero_occupancy false