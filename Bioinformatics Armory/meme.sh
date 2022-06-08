#!/bin/bash

meme rosalind_meme.txt -text -protein | grep "regular expression" -A 2 | tail -1
