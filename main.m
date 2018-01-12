% MAIN - calculate articulatory strategies and save output files. 
% 
% Tanner Sorensen
% Signal Analysis and Interpretation Laboratory
% Feb. 14, 2017

addpath util
configStruct = config;
if ~exist(configStruct.graphicsPath,'dir'), mkdir(configStruct.graphicsPath); end
articulatory_strategies(configStruct)