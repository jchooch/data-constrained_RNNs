%%% see_activity_fit_example.m

clear all;
close all;
clc;

load('ensemble/output_1.mat');

see_activity_fit(R, Adata, [1, 101, 201, 301])