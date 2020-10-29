% main
% pre-defined variables
clear
close all
main_dir='/Users/anand/Dropbox (Paladinilab)/Paladinilab Team Folder/Anand/papers/invivo1/code';
fn1='trace.mat';
fs=20000;

% figure 
f1=figure;
ax1=subplot(2,1,1);
ax2=subplot(2,1,2);
hold(ax1,'on')
hold(ax2,'on')

% loading voltage trace
tracefile=fullfile(main_dir,fn1);
load(tracefile)% loads a variable called trace

% detecting spikes and calculating spike threshold
spk=detect_spikes(tr,fs,'default','filtf',1000);

timevec=(1/fs:1/fs:1/fs*length(spk.filtr));
plot(ax1,timevec,spk.filtr*1000)
plot(ax1,spk.pks_t,spk.pks*1000,'*r')
plot(ax1,spk.thrs_t,spk.thrs*1000,'or')

% calculating vmin and vthr
vthr=spk.thrs;
vmin=NaN(length(vthr)-1,1);
for ii=1:length(spk.pks_t)-1
    curspk=round(spk.pks_t(ii)*fs);
    nxspk=round(spk.pks_t(ii+1)*fs);
    vmin(ii)=min(spk.filtr(curspk:nxspk));    
end

result_cell=detect_modes(vmin,vthr,0,ax2);

% axes decoration
xlabel(ax1,'Time (s)')
ylabel(ax1,'Voltage (V)')
xlabel(ax2,'Voltage (V)')
ylabel(ax2,'Counts')





