function tractOut=wma_loadTck(filepath)
%tract=wma_loadTck(filepath)
%
% Loads a tck and converts it to vistasoft format.  Requires mrtrix3
%

% read_mrtrix_tracts from https://github.com/MRtrix3/mrtrix3/blob/master/matlab/read_mrtrix_tracks.mtracks = read_mrtrix_tracks (filepath);
mrtrxHeaderFields=fieldnames(tracks);
dataInd=contains(mrtrxHeaderFields,'data');
dataTypeInd=contains(mrtrxHeaderFields,'datatype');

tractOut=[];
tractOut.name='tract';
tractOut.colorRgb=[  20    90   200];
tractOut.thickness=-0.5000;
tractOut.visible=1;
tractOut.seeds=[];
tractOut.seedRadius=0;
tractOut.seedVoxelOffsets=[];
tractOut.params{1,1}='mrtrix_header';
tractOut.params{1,2}={mrtrxHeaderFields{~(dataInd&~dataTypeInd)}};
tractOut.fibers=cellfun(@transpose, tracks.data, 'UniformOutput', false)';
tractOut.query_id=-1;

end
