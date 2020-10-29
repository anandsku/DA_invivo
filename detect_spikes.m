function [spk]=detect_spikes2(tr,fs,def_params,varargin)
%% Syntax
%
% [spk]=detect_spikes2(tr,fs,def_params,varargin)
%
%% Inputs  
% tr - trace in volts
% fs - acquisition frequency in Hz
% def_params - sufflix of the file containing spike detection parameters,
% if not specified, the code uses default detection params. i.e. the
% parameter file used is 'spike_det_params_default.mat.' 
%
%% Computation/Processing     
% 
%
%
% 
%
%% Outputs  
% spk is structure containing the following fields
% - pks_t - time stamps of spike peaks
% - pks - voltage value of filtr at pks_t
% - thrs_t - time stamps of spike thresholds
% - thrs - voltage value of filtr at thrs_t
% - spkofs_t - time stamps of spike offset
% - spkofs - voltage value of filtr at spkofs_t
% - spk_hts - spike height. ie voltage difference between spike threshold and peak

% - hidv_t - time stamps when minimum dvdt was crossed
% - hidv - voltage value of filtr at hidv_t
% - wof_t - time stamps of ends of spike detection windows
% - wof - voltage value of filtr at wof_t

% - filtr - filtered trace
% - dvdt - dvdt of filtr
% - d2vdt - d2vdt of filtr
% - ft_obj- digital filter object
% - det_params - detection params
% - pkval - minimum peak value for spike detection used

% - thrs_bug - logical vector indicating spikes for which threshold could not be determined
% - spkofs_bug - logical vector indicating spikes for which spike offset could not be determined

%% Assumptions
%
%
%
%
% % % Triple percentage sign indicates that the code is part of the code
% template and may be activated if necessary in later versions. 
%% Version and Author Identity Notes  
% 
% Last modified by The Big Foot on 1/1/1400
% 
% previous version:
% next version: 
%% Related procedures and functions 
% 
%
%
%
%% Detailed notes
% This code uses algorithm obtained from Strahi from Roeper lab
%
%
%
%% Processing inputs and beginning stuff

% putting in a stop for easier debugging
% dbstop if error

% processing mandatory inputs
narg_min=2;


if nargin<3
    def_params='default';
end

ff=mfilename('fullpath');
[pt,~,~]=fileparts(ff);
paramfile=fullfile(pt,['spike_det_params_' def_params '.mat']);
if exist(paramfile,'file')
   load(paramfile,'supp_inputs') 
else
    error('Specified spike detection parameters file cannot be found')
end

if nargin<narg_min
    error(['The number of inputs should at least be ' narg_min])      
else
    % processing supplementary inputs
    if ~isempty(varargin)
        if iscell(varargin{1})
            if strcmpi(varargin{1}{1},'unpackpvpairs')
                 varargin=varargin{1}(2:end);
            end
        end
    end
    supp_inputs=parse_pv_pairs(supp_inputs,varargin);
end
%% Body of the function
if ~isvector(tr) || ~isnumeric(tr)
    error('The trace should be a numeric vector. It does not seem to be one or both of them.')
end

lentr=length(tr);
dt=1/fs;
windlen=ceil(supp_inputs.max_aft_thr*fs);
if supp_inputs.dofilt
    [filtr,ft_obj]=lowpass(tr,supp_inputs.filtf,fs);  
else
    filtr=tr;
    ft_obj=[];
end
if supp_inputs.filter_acf==1
    acf=supp_inputs.acf;
    [filtr,~]=bandstop(filtr,[acf-5,acf+5],fs);
end
% due to weird effect of filtering at the edges and the resultant false
% positive spikes, we have decided to replace the first and last 1 ms with
% actual trace
filtr(1:round(fs*.001))=tr(1:round(fs*.001));
filtr(end:-1:end-round(fs*.001))=tr(end:-1:end-round(fs*.001));


dvdt=diff(filtr)/dt;
dvdt=[dvdt;mean(dvdt)];% concatenating mean to equalize length with filtr
d2vdt=diff(diff(filtr));
d2vdt=[d2vdt;mean(d2vdt);mean(d2vdt)];

tmp_spikes=find(dvdt>supp_inputs.hidvdt);

spk=struct;
spk.filtr=filtr;
spk.ft_obj=ft_obj;
spk.det_params=supp_inputs;
spk.dvdt=dvdt;
spk.d2vdt=d2vdt;


if isnan(supp_inputs.pkval)
    spk.pkval=max(filtr);
else
   spk.pkval=supp_inputs.pkval; 
end


spk.pkval=spk.pkval-supp_inputs.max_below_hipeak;

spk.pks_t=[];
spk.pks=[];
spk.hidv_t=[];
spk.hidv=[];
spk.spk_hts=[];
spk.thrs_t=[];
spk.thrs=[];
spk.spkofs_t=[];
spk.spkofs=[];
spk.thrs_bug=[];
spk.spkofs_bug=[];
spk.wof_t=[];
spk.wof=[];





if isempty(tmp_spikes)
    warning('No spikes found')    
    return
end

currind=1;
ctr=0;
while true    
   loc=tmp_spikes(currind);
   endp=loc+windlen;
   if endp>lentr-1
       endp=lentr-1;
   end
   vvec=filtr(loc:endp);
   dvvec=dvdt(loc:endp);
   if max(vvec)>=spk.pkval && min(dvvec)<supp_inputs.lodvdt && abs(max(vvec)-filtr(loc))>supp_inputs.min_spk_ht && max(vvec)<supp_inputs.hipeak && max(vvec)>supp_inputs.lopeak
       
       ctr=ctr+1;
%        {ctr,dvvec(1),min(dvvec)}
       [maxval,maxind]=max(vvec);
       maxind=(maxind+loc-1)/fs;
       spk.pks_t=[spk.pks_t;maxind]; % spike location
       spk.pks=[spk.pks;maxval]; % spike peak value
       spk.hidv_t=[spk.hidv_t;loc/fs];
       spk.hidv=[spk.hidv;filtr(loc)];
       currind=find(tmp_spikes>endp,1);
   else
       currind=currind+1;       
   end  
   if isempty(currind) || currind>length(tmp_spikes)
       break
   end
end

spk.wof_t=spk.hidv_t+supp_inputs.max_aft_thr;
spk.wof_t(spk.wof_t>lentr/fs)=lentr/fs;
spk.wof=filtr(round(spk.wof_t*fs));

% NOW FOR DETECTING SPIKE ONSET AND OFFSET
onwin=round(supp_inputs.onwin*fs);
ofwin=round(supp_inputs.ofwin*fs);

spksloc=round(spk.pks_t*fs);
nospks=length(spk.pks_t);
    
spk.thrs_t=NaN(nospks,1);
spk.thrs=NaN(nospks,1);
spk.thrs_bug=ones(nospks,1);


spk.spkofs_t=NaN(nospks,1);
spk.spkofs=NaN(nospks,1);
spk.spkofs_bug=ones(nospks,1);

spk.spk_hts=NaN(nospks,1);


for ii=1:nospks   
   
    spkind=spksloc(ii);
   
     if spkind<=onwin
        warning('The onset and offset of the spike could not be determined due to proximity to trace beginning.');
        continue         
     elseif spkind+ofwin>length(d2vdt)
         warning('The onset and offset of the spike could not be determined due to proximity to trace end.');
         continue
     end
     prevec2=(d2vdt(spkind-onwin:1:spkind-1));
     [maxpre2,~]=max(prevec2);
     
     prevec=(dvdt(spkind-onwin:1:spkind-1));
     [maxpre,~]=max(prevec);
     
     postvec2=abs(d2vdt(spkind+1:1:spkind+ofwin));
     [maxpost2,~]=max(postvec2);
     
     postvec=abs(dvdt(spkind+1:1:spkind+ofwin));
     [maxpost,~]=max(postvec);
     
     uplen=round(supp_inputs.upwin*fs); %for 2, uplen would have to be 1/10 of a ms
     
     strkon=1;
     currind=1;
     strklen=0;
     exok=-1;
     while true
         if prevec2(currind)>maxpre2*supp_inputs.onpklim && prevec(currind)>maxpre*supp_inputs.onpklim      
             currind=currind+1;
             strklen=strklen+1;
         else
             currind=currind+1;
             strkon=currind;
             strklen=0;
         end
         
         if currind>=length(prevec2)
            exok=0;
            break
         end
         
         if strklen>=uplen  
            exok=1;
            break             
         end
         
     end
     if exok==1
        spk.thrs_t(ii)=(spkind-onwin+strkon-1)/fs;
        spk.thrs(ii)=filtr(spkind-onwin+strkon-1);
        spk.spk_hts(ii)=filtr(spkind)-filtr(spkind-onwin+strkon-1);
        spk.thrs_bug(ii)=0;
     else
%         spk.thrs_t(ii)=spkind; with changes in initiailization values, this
%         is unnecessary
%         spk.thrs_bug(ii)=1;
%         spk.spk_hts(ii)=0;
        warning('The onset of the spike could not be determined.');
     end
     
     % for post 
%      mindif=0.001*fs;
     strkon=1;
     currind=1;
     strklen=0;
     exok=-1;
     while true
         if postvec(currind)<maxpost*supp_inputs.ofpklim && postvec2(currind)<maxpost2*supp_inputs.ofpklim      
             currind=currind+1;
             strklen=strklen+1;
         else
             currind=currind+1;
             strkon=currind;
             strklen=0;
         end
         
         if currind>=length(postvec2)
            exok=0;
            break
         end
         
         if strklen>=uplen %&& strkon>=mindif   
            exok=1;
            break             
         end
         
     end
     if exok==1
        spk.spkofs_t(ii)=(spkind+strkon+1)/fs;
        spk.spkofs(ii)=filtr(spkind+strkon+1);
        spk.spkofs_bug(ii)=0;       
     else        
        warning('The offset of the spike could not be determined.');
     end
   
end













