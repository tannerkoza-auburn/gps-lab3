%{
    DGPS with dynamic base length. Trimble on roof of Woltosz and Novatel
    on roof of car. Will attempt Lambda-method RTK and compare to
    traditional single-difference DGPS solution (code-only).
%}

clear; clc; close all;

%% --- Import Data
load('RCVRT.mat'); % trimble (static on woltosz roof)
load('RCVR1.mat'); % novatel on car