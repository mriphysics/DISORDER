function [theta,NS,recTypSS]=generateParametersExp(expN,quick)

%GENERATEPARAMETERSEXP   Generates the parameters for a given synthetic
%DISORDER experiment
%   [THETA,NS,RECTYPSS]=GENERATEPARAMETERSEXP(EXPN,QUICK)
%   * EXPN is the experiment descriptor
%   * {QUICK} is a flag for quick experiments, it defaults to 0
%   ** THETA is a vector with each column specifying the rotation parameter.
%   row 1 is inter-shot rotation, row 2 is inter-shot translation, row 3 
%   is intra-shot rotation, row 4 is intra-shot translation. Note only row
%   1 functionality was tested in the DISORDER manuscript, but the code is
%   designed to perform other simulations (although no error control is
%   implemented)
%   ** NS is a vector with first row specifying the number of shots (or 
%   of transform states), second row specifying the SENSE factors along
%   each direction, third row specifying the number of within shot
%   subdivisions, fourth/fifth row specifying the multiband subdivisions, 
%   sixth row specifying some random acceleration, seventh row specifying
%   the number of averages. Note only first and second row functionalities
%   where used in the DISORDER manuscript, but the code is designed to
%   perform other simulations (although no error control is implemented)
%   ** RECTYPSS is the reconstruction type. This starts by including 
%   different fields separated by "_", first the encoding type 
%   ({'Sequential','Checkered','Random','Random-checkered'}), then the 
%   number of spatial multiresolution levels ({'1','2'}), then the number 
%   of temporal multirresolution levels ({'1','2'}). Note temporal
%   multirresolution was not tested in the DISORDER manuscript, but the
%   code may accept it perhaps with minimum modifications
%

if nargin<2 || isempty(quick);quick=0;end

if strcmp(expN,'05')
    theta=[2 5 10];
    NS=[4 64];
    recTyp={'Sequential','Checkered','Random','Random-checkered'};%Type of reconstruction
    for n=1:length(recTyp);recTyp{n}=strcat(recTyp{n},'_1','_1');end
elseif strcmp(expN,'06')
    theta=[5 10 20];
    NS=[16 256];  
    recTyp={'Checkered_1','Random_1','Random-checkered_1','Checkered_2','Random_2','Random-checkered_2'};
    for n=1:length(recTyp);recTyp{n}=strcat(recTyp{n},'_1');end    
elseif strcmp(expN,'07')
    theta=[0 5 10 20];
    if ~quick
        NS=[64 16;
            1 2];
    else
        NS=[64 25 16;
             1 1.6 2];
    end
    recTyp={'Checkered_2','Random_2','Random-checkered_2'};
    for n=1:length(recTyp);recTyp{n}=strcat(recTyp{n},'_1');end
end
if quick;theta=theta(:,1:end-1);NS=NS(:,1:end-1);end

theta(end+1:4,:)=0;
NS(end+1:7,:)=1;
recTyp=recTyp(1:end);
recTypS=cell(length(recTyp));
for n=1:length(recTyp);recTypS{n}=strsplit(recTyp{n},'_');end
recTypSS=cell(length(recTypS),length(recTypS{1}));
c=1;
for s=1:length(recTypS{1})
    for n=1:length(recTypS)
        recTypSS{c}=recTypS{n}{s};
        c=c+1;
    end
end
