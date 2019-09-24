function [DF,SF,WF,SF2]=checkFilterFeasibility(N,K,DF,SF)

%CHECKFILTERFEASIBILITY checks whether prescribed filters are feasible.
%Otherwise it switches to different implementations
%   [DF,SF]=BUILDSHEARLETSYSTEM(N,K,DF,SF)
%   * N are the dimensions of the shearlet
%   * K specifies the levels of shearing for each scale
%   * DF is the directional filter
%   * SF is the scaling filter
%   ** DF is the assigned directional filter
%   ** SF is the assigned scaling filter
%   ** DF is the assigned wavelet filter
%   ** SF2 is the assigned scaling2 filter
%

WF=mirrorModulation(SF);%Wavelet filter
SF2=SF;
success=0;
for l=1:7
    [DFp,SFp]=assignFilters(DF,SF,l);
    
    L=length(SF);
    %First check
    lcheck1=L;%This is actually the wavelet filter but it is same size as scaling
    for k=1:(length(K)-1);lcheck1=L+2*(lcheck1-1);end   
    if any(lcheck1>N);continue;end
    
    %Second check
    M=size(DFp);
    lcheck2=L;
    for k=1:max(K);lcheck2=L+2*(lcheck2-1);end
    lcheck2=lcheck2+(M(1)-1)*2^(max(K)+1);
    if any(lcheck2>N) || any(M(2)>N);continue;end
    success=1;
    break;
end
DF=DFp;SF=SFp;
assert(success>0,'The specified shearlet system is not available for data of size%s . Try decreasing the number of scales an shearings',sprintf(' %d',N));
if l>1;fprintf('The specified shearlet system was not available for data of size%s. Filters were automatically set to configuration %d (see checkFilterFeasibility.m)\n',sprintf(' %d',N),l);end

function [DFp,SFp]=assignFilters(DF,SF,l)
    if l==1%Different configurations
        DFp=DF;
        SFp=SF;
    elseif l==2
        DFp=buildDirectionalFilter('dmaxflat4');
        SFp=buildQuadratureMirrorFilter('maxflat2');
    elseif l==3
        DFp=buildDirectionalFilter('cd');
        SFp=buildQuadratureMirrorFilter('maxflat2');
    elseif l==4
        DFp=buildDirectionalFilter('cd');
        SFp=buildQuadratureMirrorFilter('Coiflet0');
    elseif l==5
        DFp=buildDirectionalFilter('cd');
        SFp=buildQuadratureMirrorFilter('Daubechies1');
    elseif l==6
        DFp=buildDirectionalFilter('oqf_362');
        SFp=buildQuadratureMirrorFilter('Daubechies1');
    elseif l==7
        DFp=buildDirectionalFilter('oqf_362');
        SFp=buildQuadratureMirrorFilter('Daubechies0');
    end
end

end
