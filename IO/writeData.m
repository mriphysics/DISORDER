function rec=writeData(rec)

%WRITEDATA   Writes reconstructed data to disk and assigns geometry
%information
%   REC=WRITEDATA(REC)
%   * REC is a reconstruction structure with writing information
%   ** REC is a reconstruction structure with updated information
%

%GENERATE THE FILE NAME
fileName=generateNIIFileName(rec);
if rec.Fail || ~any(rec.Dyn.Typ2Wri)
    if nargout==0;rec=[];end%To save memory
    return;
end
rec.Dyn.Typ2Rec(ismember(rec.Dyn.Typ2Rec,[3 6]))=[];
typ2Rec=rec.Dyn.Typ2Rec;

%REMOVE THE OVERDECODING FROM THOSE DATASETS THAT MAY NEED TO
for n=typ2Rec';datTyp=rec.Plan.Types{n};    
    if ~ismember(n,[5 9 10 12 16]);rec.(datTyp)=removeOverencoding(rec.(datTyp),rec.Alg.OverDec);end
end   

%GENERATE THE GEOMETRICAL INFORMATION
[MS,MT]=generateNIIInformation(rec.Par.Mine.APhiRec);

%CORRECTIONS FOR OVERENCODING
vROI=zeros(4,1);vROI(4)=1;
overDec=rec.Alg.OverDec;
overDec(overDec>1)=1;overDec=abs(overDec);        
for n=1:3
    if overDec(n)>1
        N=size(rec.x,n);Nor=round(N/overDec(n));
        vROI(n)=vROI(n)-ceil((N-Nor)/2);
    end
end
vROI=MT*vROI;
MT(1:3,4)=vROI(1:3);

%GENERATE STRUCTURES FOR WRITING AND WRITE IMAGES
c=1;
for n=typ2Rec';datTyp=rec.Plan.Types{n};datName=rec.Plan.TypeNames{n};    
    if rec.Dyn.Typ2Wri(n)
        if n==8;rec.(datTyp)=single(abs(rec.(datTyp)));end%Mask      
        x{c}=rec.(datTyp);MSV{c}=MS;MTV{c}=MT;suff{c}=datName;        
        c=c+1;
    end
end
if c>1;writeNII(fileName,suff,x,MSV,MTV);end

%WRITE SURROGATE INFORMATION
if rec.Dyn.Typ2Wri(17)
    T=gather(rec.T);
    if isfield(rec,'DISORDER');MotionInfo=rec.DISORDER;else MotionInfo=[];end
    save(sprintf('%s_%s.mat',fileName,strcat(rec.Plan.TypeNames{17})),'T','MotionInfo');
end

if nargout==0;rec=[];end%To save memory

function [MS,MT]=generateNIIInformation(A)
    MT=A;
    MS=sqrt(sum(MT(1:3,1:3).^2,1));
    MTT=[-1 0 0 0;0 -1 0 0;0 0 1 0;0 0 0 1];
    MT=MTT*MT;
    MS(4)=1;
end

end