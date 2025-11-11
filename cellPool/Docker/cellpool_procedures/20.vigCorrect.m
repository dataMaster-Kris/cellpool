%This script is adapted from an example usage script available with BaSiC.
%https://github.com/marrlab/BaSiC/blob/master/examples/example_brainWSI.m

%The run command includes definition of variable 'files'

clc; close all;
addpath(genpath('BaSiC-master'));

files = string(transpose(split(files, ' ')));

uniqueChannels = unique(string(arrayfun(@(x) extractAfter(x, "-"), files)));
uniqueChannelSuffixes = unique(string(arrayfun(@(x) extractBefore(x, "."), uniqueChannels)));
thisField = string(arrayfun(@(x) extractBefore(x, "-C"), files(1)))
thisField = extractAfter(thisField, "_")
 
for i = 1:length(uniqueChannels)

	this_files = find(endsWith(files, uniqueChannels(i)));
	for j = 1:length(this_files)
		IF(:,:,j) = imread(files(this_files(j)));
	end
	
	%% estimate flatfield and darkfield
	% For fluorescence images, darkfield estimation is often necessary (set
	% 'darkfield' to be true)
	[flatfield,darkfield] = BaSiC(IF,'darkfield','true');

	%% image correction
	IF_corr = zeros(size(IF));
	for j = 1:length(this_files)
    		IF_corr(:,:,j) = (double(IF(:,:,j))-darkfield)./flatfield;
	end
	
	%% save images for next steps
	for j = 1:length(this_files)
		outName = extractBefore(files(this_files(j)), ".")
    		imwrite(uint16(IF_corr(:,:,j)), outName + '.vigCorr.tiff');
	end

	%Save flatfield and darkfield
	thisCh = char(uniqueChannelSuffixes(i))
	imwrite(flatfield, fullfile(thisField + '-' + thisCh + '.flat-field.tiff'));
	imwrite(darkfield, fullfile(thisField + '-' + thisCh + '.dark-field.tiff'));

end

