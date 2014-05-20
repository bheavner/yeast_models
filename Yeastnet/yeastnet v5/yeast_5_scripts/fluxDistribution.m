function fluxDistribution(filename)

model = readCbModel(filename,inf);

% define reactions of interest
rxnNames = {
    % glycolysis
    'glucose transport'
    'hexokinase (D-glucose:ATP)'
    'glucose-6-phosphate isomerase'
    'phosphofructokinase'
    'fructose-bisphosphate aldolase'
    'glyceraldehyde-3-phosphate dehydrogenase'
    'phosphoglycerate kinase'
    'phosphoglycerate mutase'
    'enolase'
    'pyruvate kinase'
    'pyruvate decarboxylase'
    'alcohol dehydrogenase, reverse rxn (acetaldehyde -> ethanol)' % should be 'alcohol dehydrogenase (ethanol)'
    'ethanol transport'
    ' ' % TCA
    'pyruvate dehydrogenase'
    'citrate synthase'
    'citrate to cis-aconitate(3-)'
    'cis-aconitate(3-) to isocitrate'
    'isocitrate dehydrogenase (NAD+)'
    'oxoglutarate dehydrogenase (lipoamide)'
    'oxoglutarate dehydrogenase (dihydrolipoamide S-succinyltransferase)'
    'succinate-CoA ligase (ADP-forming)'
    'succinate dehydrogenase (ubiquinone-6)'
    'fumarase'
    'malate dehydrogenase'
    ' ' % PPP
    'glucose 6-phosphate dehydrogenase'
    '6-phosphogluconolactonase'
    'phosphogluconate dehydrogenase'
    'ribose-5-phosphate isomerase'
    'ribulose 5-phosphate 3-epimerase'
    'transketolase 1'
    'transaldolase'
    'transketolase 2'
    ' ' % other
    };

FBAsolution = optimizeCbModel(model,[],'one');
fprintf('aerobic growth\nrate:\t%.2f\n\n',FBAsolution.f);

for k = 1:length(rxnNames)
   ind = strcmp(rxnNames{k},model.rxnNames);
   if sum(ind)>1, disp(rxnNames{k}); end
   fprintf('%.2f\t%s\n',FBAsolution.x(ind),rxnNames{k});
end
       
ind = strcmp('oxygen exchange',model.rxnNames); model.ub(ind) = 0;

ind = ismember(model.rxnNames,{...
    'lipid pseudoreaction [no 14-demethyllanosterol, no ergosta-5,7,22,24(28)-tetraen-3beta-ol]'
    'ergosterol exchange'
    'lanosterol exchange'
    'zymosterol exchange'
    'phosphatidate exchange'
    });
model.ub(ind) = inf;

ind = strcmp('lipid pseudoreaction',model.rxnNames); model.ub(ind) = 0;

FBAsolution = optimizeCbModel(model,[],'one');
fprintf('anaerobic growth\nrate:\t%.2f\n\n',FBAsolution.f);

for k = 1:length(rxnNames)
   ind = find(strcmp(rxnNames{k},model.rxnNames));
   fprintf('%.2f\t%s\n',FBAsolution.x(ind),rxnNames{k}); %#ok<FNDSB>
end