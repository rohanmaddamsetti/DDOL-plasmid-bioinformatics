#!/usr/bin/env python

""" parse-MOBtyper-results.py by Rohan Maddamsetti.

This script takes ../results/unique_mob_results.txt, which are  the MOB-typer
 results that Hye-in Son generated,  and produces a csv file with the most
 relevant infomation for analysis in R.

Usage: python parse-MOBtyper-results.py > ../results/mobility-results.csv

"""
MOBtyper_results_file = "../results/unique_mob_results.txt"

with open(MOBtyper_results_file) as MOBresults:
    print("Plasmid,PredictedMobility")
    for i, line in enumerate(MOBresults):
        if i == 0: continue ## skip the header
        line = line.strip()
        fields = line.split("\t")
        sample_id = fields[0]
        SeqID = sample_id.split("_")[-1]
        plasmid_mobility = fields[13]
        print(f"{SeqID},{plasmid_mobility}")
            
