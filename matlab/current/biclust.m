function res = biclust(data, method, varargin)
% Wrapper function to various biclustering algorithms
%
  % Using a quite useful hack at http://stackoverflow.com/questions/2775263/how-to-deal-with-name-value-pairs-of-function-arguments-in-matlab, much better than that ugly InputParser
  %
  % Inputs:
  % data - Input matrix data
  % method - Method to be used for biclustering
  % arguments for method used
         
  if nargin<2
    error('Wrong number of inputs');
  end

% count arguments
  nArgs = length(varargin);
  if round(nArgs/2)~=nArgs/2
    error('biclust: needs propertyName/propertyValue pairs')
  end

% define default options

  switch method
    case 'cc'
      options = struct('alpha',1.5, 'delta',1.0, 'numbiclust',100, 'random', ...
                       false);
    case 'plaid'
      options = struct('cluster','b', 'fit_model',['m' 'a' 'b'], 'background',true, ...
                       'background_layer',NaN, 'background_df',1, 'row_release', ...
                       0.7, 'col_release',0.7, 'shuffle', 3, 'max_layers',20, ...
                       'iter_startup',5, 'iter_layer',10, 'verbose',true);
    case 'opsm'
      options = struct('partialmodelsize',10);

    case 'isa'
      options = struct('thr_row',1:0.5:3, 'thr_col',1:0.5:3, 'nseeds',100, ...
                       'direction', {{'updown', 'updown'}});
    case 'bimax'
      options = struct('minnorows',0, 'minnocols',0, 'numbiclust', 100);
    case 'kSpectral'
      options = struct('normalization','log', 'numbereigenvalues',3, 'minr',2, ...
                       'minc',2, 'withinvar',1);
    case 'floc'
      options = struct('numbiclust', 20,'prow', 0.5, 'pcolumn', 0.5, 'resth', ...
                       [], 'minrow', 8, 'mincol', 6, 'niter', 500, 'blocrow', ...
                       [], 'bloccolumn', []);
    case 'bayes'
      options = struct('numclusters',100, 'normchoice',3, 'alpha',10);
    case 'xmotif'

      options = struct('ns',10, 'nd',10, 'sd',5, 'alpha',0.05, 'number',100);
    case 'itl'
      options = struct('numclust',100, 'display', 0); %check
      
    case 'las'
      options = struct('numbc', 100, 'threads', 1,...
                       'iterationsperbc',1000, 'scorethreshold',1);
    case 'bsgp'
      options = struct('n',100, 'display', 0);
    case 'qubic'
      options = struct('quantile',0.06,'rank',1,'consistency',0.95,...
                        'output_bicluster',100,'filter_proportion',1);
    case 'fabia'
      options = struct('bicluster_no',5,'sparseness_factor',0.1,'iteration_no',500,...
                        'spl',0,'spz',0.5,'non_negative',0,'random',1,'center',2,'norm',1,...
                        'scale',0,'lap',1,'nL',0,'lL',0,'bL',0);
      case 'fuzzyBiclust'
         disp('fuzzy algorithm running');

    otherwise
      error('method: Wrong input');
  end
  
% read the acceptable names
switch method
    case 'fuzzyBiclust'
    disp('fuzzy algorithm running');
    otherwise
    
  optionNames = fieldnames(options);

      
  for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
    inpName = lower(pair{1}); % make case insensitive

    if any(strmatch(inpName,optionNames))
      %# overwrite options. 
      options.(inpName) = pair{2};
    else
      error('%s is not a recognized parameter name',inpName)
    end
  end
end
% call method accordingly
  switch method
    case 'cc'
      res = cc(data, options.alpha, options.delta, options.numbiclust, ...
               options.random);
    case 'plaid'
      res = plaid(data, options.cluster, options.fit_model, options.background, ...
                  options.background_layer, options.background_df, ...
                  options.row_release, options.col_release, options.shuffle, ...
                  options.max_layers, options.iter_startup, options.iter_layer, ...
                  options.verbose);
    case 'opsm'
      res = opsm(data, options.partialmodelsize);
    case 'isa'
      res = itersa(data, options.thr_row, options.thr_col, options.nseeds, ...
                   options.direction);
    case 'bimax'
      res = bimax(data, options.minnorows, options.minnocols, ...
                  options.numbiclust);
    case 'kSpectral'
      res = kSpectral(data, options.normalization, options.numbereigenvalues, ...
                      options.minr, options.minc, options.withinvar);
    case 'floc'
      res = floc(data, options.numbiclust, options.prow, options.pcolumn, ...
                 options.resth, options.minrow, options.mincol, options.niter, ...
                 options.blocrow, options.bloccolumn);
    case 'bayes'
      res = BBC(data, options.numclusters, options.normchoice, options.alpha);
    case 'xmotif'
      res = xmotif(data, options.ns, options.nd, options.sd, options.alpha, ...
                   options.number);
    case 'itl'
      res = ITL(data, options.numclust, options.display); % Check
    case 'las'
      res = LAS(data, options.numbc, options.threads,...
                options.iterationsperbc, options.scorethreshold);
    case 'bsgp'
      res = spectralCoClustering(data, options.n, options.display);
    case 'qubic'
      res = qubic(data,options.quantile,options.rank,options.consistency,...
                    options.output_bicluster,options.filter_proportion);
    case 'fabia'
      res = fabia(data,options.bicluster_no,options.sparseness_factor,options.iteration_no,...
                    options.spl,options.spz,options.non_negative,options.random,options.center,...
                    options.norm,options.scale,options.lap,options.nL,options.lL,options.bL);
      case 'fuzzyBiclust'
          res = fuzzyBiclust(data);
          
    otherwise
      error('method: Wrong input');
  end
end
