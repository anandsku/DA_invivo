function result_cell=detect_modes(vhyp,vthr,vthrbi_flag,plotax)
  % vhyp - vector containing all Vmin values for a given trace
  % vthr - vector containing all spike threshold values for a cell
  % vthrbi_flag - flag variable indicating that the distribution is bimodal
  % in Vthr (Vthr BC>0.53)
  % plotax - axes object where results should be plotted
  % result_cell - indicates the values of mode1, mode2 and threshold value
  % used to demarcate rebound and plateau events
  
  binwidth=0.00025;

   
   m3=skewness(vhyp,0);
   m4=kurtosis(vhyp,0)-3;
   nn=length(vhyp);
   bc_hyp=(m3^2+1)/(m4+3*(nn-1)^2/((nn-2)*(nn-3)));
       
   m3=skewness(vthr,0);
   m4=kurtosis(vthr,0)-3;
   nn=length(vthr);
   bc_thr=(m3^2+1)/(m4+3*(nn-1)^2/((nn-2)*(nn-3)));
   
   yyaxis(plotax,'right')
   vlo=min(vthr);
   vhi=max(vthr);
   eds=(vlo-binwidth:binwidth:vhi+binwidth);
   [cts,~]=histcounts(vthr,eds); 
   l1=histogram(plotax,'BinEdges',eds*1000,'BinCounts',cts);
   
   
   yyaxis(plotax,'left')
   vlo=min(vhyp);
   vhi=max(vhyp);
   eds=(vlo-binwidth:binwidth:vhi+binwidth);
   [cts,~]=histcounts(vhyp,eds); 
   l2=histogram(plotax,'BinEdges',eds*1000,'BinCounts',cts) ;
   
   bc_hypstr=sprintf('%0.3f',bc_hyp);
   bc_thrstr=sprintf('%0.3f',bc_thr);
   
  
  
  txx=sprintf(['Vmin BC = ' bc_hypstr ' Vthr BC = ' bc_thrstr]);
  title(plotax,txx)
          
   

  
 
    %% fitting two Gaussians
  fitcell=cell(0);
  
  if vthrbi_flag==1
      fitcell{1}=vhyp;
      fitcell{2}=vthr;
      lrswitch=1;
  else
       fitcell{1}=vhyp; 
       lrswitch=0;
  end
  
result_cell={'Mode1','Mode2','Mid-threshold'};
  for lk=1:length(fitcell)
   fitvec=fitcell{lk};
  pdf_normmixture = @(fitvec,p,mu1,mu2,sigma1,sigma2) ...
                         p*normpdf(fitvec,mu1,sigma1) + (1-p)*normpdf(fitvec,mu2,sigma2);

%    pdf_normmixture = @(fitvec,mu1,mu2,sigma1,sigma2) ...
%                           normpdf(fitvec,mu1,sigma1) + normpdf(fitvec,mu2,sigma2);

   pdf_uni=@(fitvec,p,mu1,sigma1) ...
                         p*normpdf(fitvec,mu1,sigma1);
                    
    pStart=.5;
    muStart=quantile(fitvec,[.15 .85]);
    
      
    sigmaStart=sqrt(var(fitvec) - .1*diff(muStart).^2);
    start=[pStart muStart sigmaStart sigmaStart];
    
    lb = [0 -Inf -Inf 0 0];
    ub = [1 Inf Inf Inf Inf];
    options = statset('MaxIter',1000, 'MaxFunEvals',2000);

    paramEsts=mle(fitvec, 'pdf',pdf_normmixture, 'start',start, ...
                          'lower',lb, 'upper',ub,'options',options);
                      
   vlo=min(fitvec);
   vhi=max(fitvec);
   eds=(vlo-binwidth:binwidth:vhi+binwidth);   
   pdfgrid = pdf_normmixture(eds,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
   pdfgrid1=pdf_uni(eds,paramEsts(1),paramEsts(2),paramEsts(4));
   pdfgrid2=pdf_uni(eds,1-paramEsts(1),paramEsts(3),paramEsts(5));
   
   [~,maxpdfind]=max(pdfgrid);
   maxx=eds(maxpdfind);

if max(pdfgrid2)>max(pdfgrid1)
    [midp]=calculate_optimal_threshold(pdfgrid1',pdfgrid2',eds',1);
else
    [midp]=calculate_optimal_threshold(pdfgrid2',pdfgrid1',eds',1);
end
  
if midp==-1
   error('rr') 
end
   
   core_stats=[];
   fringe_stats=[];
   
       if lk==2 || lrswitch==1
%            [~,maxx2ind]=max(pdfgrid(midpind:end));
%            maxx2=eds(midpind+maxx2ind-1);    
            maxx2=paramEsts(3);
            core_stats=[paramEsts(2),paramEsts(4)];
            fringe_stats=[paramEsts(3),paramEsts(5)];
             maindist=fitvec(fitvec<midp);
              q1=prctile(maindist,25);
              q3=prctile(maindist,75);
              iq=q3-q1;
              iqfac2=(maxx2-q3)/iq;
              iqfac=(midp-q3)/iq;

       else
%               [~,maxx2ind]=max(pdfgrid(1:midpind));
%               maxx2=eds(maxx2ind);
              maxx2=paramEsts(2);
            core_stats=[paramEsts(3),paramEsts(5)];
            fringe_stats=[paramEsts(2),paramEsts(4)];
              maindist=fitvec(fitvec>midp);
              q1=prctile(maindist,25);
              q3=prctile(maindist,75);
              iq=q3-q1;
              iqfac2=(q1-maxx2)/iq;
              iqfac=(q1-midp)/iq;
       end
       midp=(midp+maxx2)/2;
%        midp=midp-(abs(maxx2)-abs(midp))/3; % p37_2    
       plot(plotax,eds*1000,pdfgrid*binwidth*length(fitvec),'k-');
       plot(plotax,eds*1000,pdfgrid1*binwidth*length(fitvec),'g-');
       plot(plotax,eds*1000,pdfgrid2*binwidth*length(fitvec),'r-');
       l3=plot(plotax,[midp*1000,midp*1000],ylim(plotax),'k--');
       l4=plot(plotax,[maxx2*1000,maxx2*1000],ylim(plotax),'m--');
       l5=plot(plotax,[maxx*1000,maxx*1000],ylim(plotax),'m--');
       meth='2gauss';
       result_cell=[result_cell;{maxx,maxx2,midp}];   
       legend(plotax,[l2,l1,l5,l4,l3],{'Vmin','Vthr','Mode 1','Mode 2','Threshold'})
              

    
  end 
 
end



