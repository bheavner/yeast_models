% ks 26 may 2011

clc
% clearvars

strMatchExact = @(str,strArray)find(strcmp(str,strArray));

%% yeast 4

filename = 'yeast_4.05';

m = TranslateSBML([filename,'.xml']); % from http://yeast.sf.net/

% glucose uptake => 1
j = strMatchExact('r_1293',{m.reaction.id}); 
m.reaction(j).kineticLaw.parameter(2).value = 1;
% oxygen uptake => inf
j = strMatchExact('r_1435',{m.reaction.id}); 
m.reaction(j).kineticLaw.parameter(2).value = inf;

% ATP requirement => 0
j = strMatchExact('r_0249',{m.reaction.id}); 
m.reaction(j).kineticLaw.parameter(1).value = 0;
m.reaction(j).kineticLaw.parameter(2).value = inf;

% save and test
OutputSBML(m,[filename,'_modified.xml']);
% analyse_model([filename,'_modified']);

%% iIN800

copyfile('iIN800.xml','iIN800_valid.xml'); % from http://129.16.106.142/models.php?c=S.cerevisiae, *invalid SBML*
replaceTextInFile('<html xmlns="http://www.w3.org/1999/xhtml">','<body xmlns="http://www.w3.org/1999/xhtml">','iIN800_valid.xml');
replaceTextInFile('</html>','</body>','iIN800_valid.xml')

filename = 'iIN800_valid';

m = TranslateSBML([filename,'.xml']); 

% add a _b to all boundary metabolites

J = find([m.species.boundaryCondition]);
oldIDs = {m.species.id};
for k = J(:)'
    m.species(k).id = [m.species(k).id,'_b'];
end
newIDs = {m.species.id};
for k = 1:length(m.reaction)
    s = {'reactant','product','modifier'};
    for kk = 1:length(s)
        ss = s{kk};
        for kkk = 1:length(m.reaction(k).(ss))
            j = strMatchExact(m.reaction(k).(ss)(kkk).species,oldIDs);
            m.reaction(k).(ss)(kkk).species = newIDs{j};
        end
    end
end

% boundary metabolites to end
[dummy,j] = sort([m.species.boundaryCondition]); %#ok<ASGLU>
m.species = m.species(j);

% assign growth
j = strMatchExact('Growth',{m.reaction.name});
m.reaction(j).kineticLaw.parameter(3).value = 1;
m.reaction(j).product(:) = [];

% define growth medium

J = strmatch('Uptake of',{m.reaction.name});
for k = J(:)'
    m.reaction(k).kineticLaw.parameter(1).value = 0;
    str = strrep(m.reaction(k).name,'Uptake of ','');
    switch str
        case 'alpha-D-glucose'
            m.reaction(k).kineticLaw.parameter(2).value = 1;            
            
        case {'oxygen','NH3','H(+)','phosphate','potassium','sodium','sulfate'}
            % no iron nor water in model
            m.reaction(k).kineticLaw.parameter(2).value = inf;            
            
        otherwise
            m.reaction(k).kineticLaw.parameter(2).value = 0;
    end
end

J = strmatch('Excretion of',{m.reaction.name});
for k = J(:)'
    m.reaction(k).kineticLaw.parameter(1).value = 0;
    m.reaction(k).kineticLaw.parameter(2).value = inf;
end

% add gene associations
gene_map = iIN800_gene_map;

for k = 1:length(m.reaction)
    j = strMatchExact(m.reaction(k).id(3:end),gene_map(:,1));
    if ~isempty(j)
        gene = gene_map{j,2};
        gene = strrep(gene,':',' and ');
        gene = strrep(gene,'-','_');
    else
        gene = '';
    end
    m.reaction(k).notes = strrep(m.reaction(k).notes,'<notes>','');
    m.reaction(k).notes = strrep(m.reaction(k).notes,'</body>','');
    m.reaction(k).notes = strrep(m.reaction(k).notes,'</notes>','');
    m.reaction(k).notes = [m.reaction(k).notes,'<p>GENE_ASSOCIATION:',gene,'</p></body>'];
end

% save and test
OutputSBML(m,[filename,'_modified.xml']);
% analyse_model([filename,'_modified']);


%% iFF708

copyfile('iFF708.xml','iFF708_valid.xml'); % from http://129.16.106.142/models.php?c=S.cerevisiae, *invalid SBML*
replaceTextInFile('<html xmlns="http://www.w3.org/1999/xhtml">','<body xmlns="http://www.w3.org/1999/xhtml">','iFF708_valid.xml');
replaceTextInFile('</html>','</body>','iFF708_valid.xml')

filename = 'iFF708_valid';

m = TranslateSBML([filename,'.xml']); % from http://129.16.106.142/models.php?c=S.cerevisiae

% add a _b to all boundary metabolites

J = find([m.species.boundaryCondition]);
oldIDs = {m.species.id};
for k = J(:)'
    m.species(k).id = [m.species(k).id,'_b'];
end
newIDs = {m.species.id};
for k = 1:length(m.reaction)
    s = {'reactant','product','modifier'};
    for kk = 1:length(s)
        ss = s{kk};
        for kkk = 1:length(m.reaction(k).(ss))
            j = strMatchExact(m.reaction(k).(ss)(kkk).species,oldIDs);
            m.reaction(k).(ss)(kkk).species = newIDs{j};
        end
    end
end

% boundary metabolites to end
[dummy,j] = sort([m.species.boundaryCondition]); %#ok<ASGLU>
m.species = m.species(j);

% assign growth
j = strMatchExact('Growth',{m.reaction.name});
m.reaction(j).kineticLaw.parameter(3).value = 1;
m.reaction(j).product(:) = [];

% define growth medium

J = strmatch('Uptake of',{m.reaction.name});
for k = J(:)'
    m.reaction(k).kineticLaw.parameter(1).value = 0;
    str = strrep(m.reaction(k).name,'Uptake of ','');
    switch str
        case 'alpha-D-glucose'
            m.reaction(k).kineticLaw.parameter(2).value = 1;            
            
        case {'oxygen','NH3','H(+)','phosphate','potassium','sodium','sulfate'}
            % no iron nor water in model
            m.reaction(k).kineticLaw.parameter(2).value = inf;            
            
        otherwise
            m.reaction(k).kineticLaw.parameter(2).value = 0;
    end
end

J = strmatch('Excretion of',{m.reaction.name});
for k = J(:)'
    m.reaction(k).kineticLaw.parameter(1).value = 0;
    m.reaction(k).kineticLaw.parameter(2).value = inf;
end

% add gene associations
gene_map = iFF708_gene_map;

for k = 1:length(m.reaction)
    j = strMatchExact(m.reaction(k).id(3:end),gene_map(:,1));
    if ~isempty(j)
        gene = gene_map{j,2};
        gene = strrep(gene,':',' and ');
        gene = strrep(gene,'-','_');
    else
        gene = '';
    end
    m.reaction(k).notes = strrep(m.reaction(k).notes,'<notes>','');
    m.reaction(k).notes = strrep(m.reaction(k).notes,'</body>','');
    m.reaction(k).notes = strrep(m.reaction(k).notes,'</notes>','');
    if isempty(m.reaction(k).notes)
        m.reaction(k).notes = '<body xmlns="http://www.w3.org/1999/xhtml">';
    end
    m.reaction(k).notes = [m.reaction(k).notes,'<p>GENE_ASSOCIATION:',gene,'</p></body>'];
end

% save and test
OutputSBML(m,[filename,'_modified.xml']);
% analyse_model([filename,'_modified']);

%% iMM904

copyfile('iMM904.xml','iMM904_known.xml'); % from http://129.16.106.142/models.php?c=S.cerevisiae, *invalid SBML*
replaceTextInFile('or  UNKNOWN','','iMM904_known.xml');

filename = 'iMM904_known';

m = TranslateSBML([filename,'.xml']); % from http://gcrg.ucsd.edu/In_Silico_Organisms/Yeast

% glucose uptake => 1
j = strMatchExact('R_D_Glucose_exchange',{m.reaction.name}); 
m.reaction(j).kineticLaw.parameter(1).value = -1;

% ATP requirement => 0
j = strMatchExact('R_ATP_maintenance_requirement',{m.reaction.name}); 
m.reaction(j).kineticLaw.parameter(1).value = 0;
m.reaction(j).kineticLaw.parameter(2).value = inf;

% save and test
OutputSBML(m,[filename,'_modified.xml']);
% analyse_model([filename,'_modified']);