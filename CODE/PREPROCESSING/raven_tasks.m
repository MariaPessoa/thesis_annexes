% cd('C:\Users\liamo\Documents\BIOINF\PRETHESIS\RAVEN\installation');
% checkInstallation;
initCobraToolbox(false); % false, as we don't want to update
changeCobraSolver ('ibm_cplex', 'all'); 
setRavenSolver('cobra');

root_folder = "C:\Users\liamo\Documents\BIOINF\PRETHESIS\"
essentialTasks = parseTaskList('C:\Users\liamo\Documents\BIOINF\PRETHESIS\Human-GEM\data\metabolicTasks\metabolicTasks_Essential.xlsx');
fullTasks = parseTaskList('C:\Users\liamo\Documents\BIOINF\PRETHESIS\Human-GEM\data\metabolicTasks\metabolicTasks_Full.xlsx');

%H1-CSMS task evaluation
h1_csms = load(root_folder+'DATA\Human1_Publication_Data_Scripts\tINIT_GEMs\run_tINIT_outputs\GTEx\tINIT_GTEx_outputs.mat');
%cd('C:\Users\liamo\Documents\BIOINF\PRETHESIS\Human-GEM\code') 
import compareModelField.* %extracted from compareMultipleModels
[id,binary_matrix] = compareModelField(h1_csms.INIT_output.model,"rxns");

writematrix(binary_matrix,root_folder+'DATA\H1_CSMS_RC.csv');
writecell(h1_csms.INIT_output.id,root_folder+'DATA\H1_CSMS_RC_ids.csv');
writecell(id,root_folder+'DATA\H1_CSMS_RC_rxns.csv');

useModels = sort({'adipose tissue', 'brain','breast', 'blood', 'colon', 'kidney', 'liver', 'lung', ...
       'muscle', 'pancreas', 'stomach'});
keep = ismember(h1_csms.INIT_output.id, useModels);
h1_csms11 = h1_csms.INIT_output.model(keep);

for i = 1:numel(h1_csms11) %from compareMultipleModels again
    fprintf('\n Checking model # %.0f',i)
    taskReport{i} = checkTasks(h1_csms11{i},[],false,false,false,fullTasks);
end 

taskMatrix = zeros(length(taskReport{1}.ok),numel(taskReport));
for i = 1:numel(taskReport)
    taskMatrix(:,i) = taskReport{i}.ok;
end
% writematrix(taskMatrix,root_folder+'DATA\H1_CSMS_11esstasks_MAT.csv');
writematrix(taskMatrix,root_folder+'DATA\H1_CSMS_11fulltasks_MAT.csv');

%base model task evaluation
cd('C:\Users\liamo\Documents\BIOINF\PRETHESIS\Human-GEM\code') 
import addBoundaryMets.*
ihuman = importModel('C:\Users\liamo\Documents\BIOINF\PRETHESIS\Human-GEM\model\Human-GEM.xml',false) %sbml only
model = ravenCobraWrapper(ihuman);
ihuman = addBoundaryMets(ihuman); %13096
% writecell(ihuman.rxns,'C:\Users\liamo\Documents\BIOINF\PRETHESIS\DATA\orderReactions.csv')

%get reactions required for tasks
% [taskReport, essentialReactions, taskStructure]=checkTasks(ihuman, [], false, false, true, essentialTasks);

% size(essentialReactions) %reactions x tasks
% writematrix(essentialReactions,'C:\Users\liamo\Documents\BIOINF\PRETHESIS\DATA\essentialReactions_v2.csv')

%task evaluation 
[taskReport, essentialReactions, taskStructure]=checkTasks(ihuman, [], false, false, false, essentialTasks);
[taskReport2, essentialReactions2, taskStructure2]=checkTasks(ihuman, [], false, false, false, fullTasks);

writetable(struct2table(taskReport),root_folder+'DATA\Human-GEMxml_essMAT.csv');
%writetable(struct2table(taskReport2),root_folder+'DATA\Human-GEMxml_fullMAT.csv');

