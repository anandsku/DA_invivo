function [threshold]=calculate_optimal_threshold(density_target,density_distractors,dists_vec,distractor_factor,varargin)
%% Syntax
%
% [threshold]=calculate_optimal_threshold(density_target,density_distractors,dists_vec,distractor_factor,varargin)
%
%% Inputs  
%
% density_target -  the probability density distribution for target distances 
% 
% density_distractors -the probability density distribution for distractor distances 
%
% dists_vec - the domain on which the above two distributions is
% calculated. this is the discretized distance axis. 
% 
% distractor_factor - the factor by which the probability density function
% of the distractor distances should be multiplied while determining the
% threshold.**1
%
%
% supp_inputs.target_factor - the factor by which the probability density function
% of the target distances is multipled. similar to distractor_factor.**2. default=1.   
% 
% 
% supp_inputs.densities_ratio_comparison_standard -  the number which is
% subtracted from the ratio of density values. when this value is 1 (the
% default), it corresponds to picking the crossing of the distributions as
% the threshold. default=1.  
% 
% 
% 
%% Computation/Processing     
% 
% It simply computes the point (along the distance axis) at which the density_target
% and the density_distractors cross. That is an optimal value for the threshold.  
%
% 
%
%% Outputs  
% 
% threshold - optimal value of threshold for discriminting between target
% chunks and distractor chunks. 
%
%
%% Assumptions
%
%
%
% % % Triple percentage sign indicates that the code is part of the code
% template and may be activated if necessary in later versions. 
%% Version and Author Identity Notes  
% 
% Last modified by Anand S Kulkarni on 
% 
% previous version:
% next version: 
%% Related procedures and functions 
% 
%
%
%
%% Detailed notes
%
%**1 - this effectively determines where the threshold gets placed by
% shifting the crossing point of the two distributions (in the valley). lower 
% values shift the crossing point to the right and higher values shift it
% to the left. 
%
%**2 - the target factor is by default kept at 1. Moving around of threhsold is
% achieved by controlling the distractor_factor itself. 
%
%% Processing inputs and beginning stuff

% putting in a stop for easier debugging
dbstop if error

% processing mandatory inputs
narg_min=4;

if nargin<narg_min
     error(['The number of inputs should at least be ' narg_min])
end

% processing supplementary inputs

% Assigning default values to supplementary inputs
supp_inputs.densities_ratio_comparison_standard=1;
supp_inputs.target_factor=1;
supp_inputs.write_to_disk_q=0; % should the function write a file to disk containing its output  
supp_inputs.disk_write_dir='';

supp_inputs=parse_pv_pairs(supp_inputs,varargin);

%

%% Body of the function

% multiplying densities by their scaling factor
density_target=supp_inputs.target_factor*density_target;
density_distractors=distractor_factor*density_distractors;

% finding out the peaks of both the densitie distributions
inflections_target_inds=find(exp(diff(log(diff(density_target))))<0)+1;
inflections_distractors_inds=find(exp(diff(log(diff(density_distractors))))<0)+1;
in_peak_ind=inflections_target_inds(end);
out_peak_ind=inflections_distractors_inds(1);

if in_peak_ind>out_peak_ind 
    tm=out_peak_ind;
    out_peak_ind=in_peak_ind;
    in_peak_ind=tm;
end

% ratio of the density functions between their respective peak. we subtract one from it 
densities_ratio=(density_target(in_peak_ind:out_peak_ind)./density_distractors(in_peak_ind:out_peak_ind))...
                -supp_inputs.densities_ratio_comparison_standard;
% if the ratios at either end are incompatible with the chosen densities_ratio_comparison_standard            
% if supp_inputs.densities_ratio_comparison_standard>1
%    if densities_ratio(1)<0 
%        threshold=-1;
%        return
%    end
% else
%     if densities_ratio(end)> 0
%         threshold=-1;
%        return
%     end
% end

% calculating the index of the point where the distributions cross
crossing_ind=find(exp(diff(log(densities_ratio)))<0);
cross_pre_ind=crossing_ind+in_peak_ind-1;
cross_post_ind=cross_pre_ind+1;

% calculating threshold by linear interpolation
try
    threshold=polyxpoly([dists_vec(cross_pre_ind),dists_vec(cross_post_ind)],[densities_ratio(crossing_ind),densities_ratio(crossing_ind+1)],...
                         [dists_vec(cross_pre_ind),dists_vec(cross_post_ind)],[0,0]);
catch errr
    if strcmpi(errr.message, 'Undefined function ''polyxpoly'' for input arguments of type ''double''.')
        threshold=intersections([dists_vec(cross_pre_ind),dists_vec(cross_post_ind)],[densities_ratio(crossing_ind),densities_ratio(crossing_ind+1)],...
                         [dists_vec(cross_pre_ind),dists_vec(cross_post_ind)],[0,0]);
    else
        rethrow(errr)
    end    
end
if isempty(threshold)
   error('An optimal threshold could not be found. The shape of the distributions may be abnormal') 
end
