function writeRaw(rec)

%WRITERAW   Writes raw data before SENSE reconstruction
%   WRITERAW(REC)
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the enoding information (rec.Enc)
%

fileName=generateNIIFileName(rec);
if isfield(rec,'x') && ismember(rec.Alg.SaveRaw,1:2);rec=rmfield(rec,'x');end
if isfield(rec,'X') && ismember(rec.Alg.SaveRaw,1:2);rec=rmfield(rec,'X');end
if isfield(rec,'S') && rec.Dyn.Batch~=1;rec=rmfield(rec,'S');end
if isfield(rec,'M') && rec.Dyn.Batch~=1;rec=rmfield(rec,'M');end
if ismember(rec.Alg.SaveRaw,[2 4])
    rec.Names=[];
    rec.Par.Filename=[];
    rec.Par.Scan.Date=[];
    rec.Par.Labels.Date=[];
    rec.Par.Labels.Time=[];
end 
rec=gatherStruct(rec);
   
save(sprintf('%s%03d.mat',fileName,rec.Dyn.Batch),'rec','-v7.3');

function x=gatherStruct(x)
    if isstruct(x)
        y=fieldnames(x);
        for n=1:length(y);x.(y{n})=gatherStruct(x.(y{n}));end
    elseif iscell(x)
        x=gatherCell(x);
    elseif isnumeric(x)
        x=gather(x);
    end
end

function x=gatherCell(x)
    for n=1:length(x)
        if iscell(x{n});x{n}=gatherCell(x{n});elseif isnumeric(x{n});x{n}=gather(x{n});end
    end
end
 
end
