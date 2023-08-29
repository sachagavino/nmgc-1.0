#!/bin/bash

mkdir 001AU
mkdir 002AU
mkdir 003AU
mkdir 004AU
mkdir 005AU
mkdir 010AU
mkdir 020AU
mkdir 030AU
mkdir 040AU
mkdir 060AU
mkdir 100AU
mkdir 200AU
mkdir 300AU
mkdir 400AU

cp initial_files/* 001AU
cp initial_files/* 002AU
cp initial_files/* 003AU
cp initial_files/* 004AU
cp initial_files/* 005AU
cp initial_files/* 010AU
cp initial_files/* 020AU
cp initial_files/* 030AU
cp initial_files/* 040AU
cp initial_files/* 060AU
cp initial_files/* 100AU
cp initial_files/* 200AU
cp initial_files/* 300AU
cp initial_files/* 400AU

cp physical_structure/001AU.txt 001AU/1D_static.dat
cp physical_structure/002AU.txt 002AU/1D_static.dat
cp physical_structure/003AU.txt 003AU/1D_static.dat
cp physical_structure/004AU.txt 004AU/1D_static.dat
cp physical_structure/005AU.txt 005AU/1D_static.dat
cp physical_structure/010AU.txt 010AU/1D_static.dat
cp physical_structure/020AU.txt 020AU/1D_static.dat
cp physical_structure/030AU.txt 030AU/1D_static.dat
cp physical_structure/040AU.txt 040AU/1D_static.dat
cp physical_structure/060AU.txt 060AU/1D_static.dat
cp physical_structure/100AU.txt 100AU/1D_static.dat
cp physical_structure/200AU.txt 200AU/1D_static.dat
cp physical_structure/300AU.txt 300AU/1D_static.dat
cp physical_structure/400AU.txt 400AU/1D_static.dat

cp physical_structure/parameters_001AU.txt 001AU/parameters.in
cp physical_structure/parameters_002AU.txt 002AU/parameters.in
cp physical_structure/parameters_003AU.txt 003AU/parameters.in
cp physical_structure/parameters_004AU.txt 004AU/parameters.in
cp physical_structure/parameters_005AU.txt 005AU/parameters.in
cp physical_structure/parameters_010AU.txt 010AU/parameters.in
cp physical_structure/parameters_020AU.txt 020AU/parameters.in
cp physical_structure/parameters_030AU.txt 030AU/parameters.in
cp physical_structure/parameters_040AU.txt 040AU/parameters.in
cp physical_structure/parameters_060AU.txt 060AU/parameters.in
cp physical_structure/parameters_100AU.txt 100AU/parameters.in
cp physical_structure/parameters_200AU.txt 200AU/parameters.in
cp physical_structure/parameters_300AU.txt 300AU/parameters.in
cp physical_structure/parameters_400AU.txt 400AU/parameters.in
