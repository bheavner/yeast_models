function modelToReconstruction(filename)

% a function for converting an SBML model to a reconstruction, using SBO
% terms. reliant on the SBML toolbox.

% kieran smallbone and ben heavner: 18 aug 11

if ~strncmp(fliplr(filename),fliplr('.xml'),length('.xml'))
    filename = [filename,'.xml'];
end

m = TranslateSBML(filename);

% only real reactions
J = ismember([m.reaction.sboTerm],[395,397]);
m.reaction(J) = [];

% no modelling parameters
for k = 1:length(m.reaction)
    m.reaction(k).kineticLaw(1:end) = [];
end

% % no modelling notes
% for k = 1:length(m.speciesType)
%     m.speciesType(k).notes = '';
% end
% for k = 1:length(m.reaction)
%     m.reaction(k).notes = '';
% end

% find unused species
used_species = {};
for k = 1:length(m.reaction)
    used_species = union(used_species,{m.reaction(k).reactant.species});
    used_species = union(used_species,{m.reaction(k).product.species});
    used_species = union(used_species,{m.reaction(k).modifier.species});
end
J = ~ismember({m.species.id},used_species);
m.species(J) = [];

% ... and species types
J = ~ismember({m.speciesType.id},{m.species.speciesType});
m.speciesType(J) = [];

% ... and compartments
J = ~ismember({m.compartment.id},{m.species.compartment});
m.compartment(J) = [];

% save
filename = strrep(filename,'.xml','_reconstruction.xml');
OutputSBML(m,filename);

% validate
[~,errors] = TranslateSBML(filename,1,0);
if isempty(errors)
    disp('Results: This document is valid SBML!');
else
    for k = 1:length(errors)
        fprintf('%s:\tline %g:\t(SBML Validation Rules #%g)\t%s',...
            strtrim(errors(k).severity),errors(k).line,errors(k).errorId,errors(k).message);
    end
end