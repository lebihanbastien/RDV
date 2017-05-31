function RGB = rgb( index_or_name )
% RGB defined hard-coded colors in MATLAB RGB format.

colors= {
0    0    0    'black';
0.6  0    0    'dark red';
0.25 0.25 0.9  'cobalt blue';
0    0.5  0    'dark green';
1    0.6  0    'orange';
0.3  0.8  0    'avocado';
0    1    1    'cyan';
1    0    0    'red';
1    0.2  0.2  'light red';
0.5  0.5  1    'light blue';
0    1    0    'light green';
0.8  0.5  0    'brown';
0.5  0.5  0.5  'dark gray';
1    1    0.6  'cream';
0    0.6  0.3  'super green';
1    0.5  0.5  'peach';
1    1    0    'yellow';
0    0    0.8  'dark blue';
0.8  0.8  0.8  'gray';
0.5  0    0.9  'purple';
0.3  0.8  0    'avocado';
1    0.5  1    'magenta';
0    0.8  0.8  'aqua';
0.9  0.75 0    'gold';
1    1    1    'white';
};

% Extract names from rightmost column of colors cell array:
names= colors(:,4);

% Strip off rightmost column of colors cell array and convert remaining
% columns to a matrix:
RGBs= cell2num(colors(:,1:3));


% Section 2: Convert index_or_name to an RGB triple.

if nargin ~= 1
   error('This function must be called with exactly one argument.');
end

if isnumeric(index_or_name)
   index= index_or_name;


   if length(index)==3 & all(index>=0)

      % If contents of index_or_name are an RGB triple with elements in
      % [0,1], return them as output without modification:

      if all(index<=1)
         RGB= index;
         return
      end

      % If contents of index_or_name are an RGB triple with elements in
      % [0,255], scale by 1/255 and return as output:

      if all(index<=255)
         RGB= index/255;
         return
      end

   end

   if length(index) > 1
      error('When calling with a color index, specify a single number.');
   end
   if ismember(index,1:21)
      RGB= RGBs(index,:);
      return
   end
   error('A color index must be a whole number between 1 and 21.');
end

if isa(index_or_name,'char')
   index = find(strcmp(names, index_or_name));
   if ~isempty(index)
      RGB= RGBs(index,:);
   else
      fprintf(2, ['Warning: Unknown color name "%s".  ' ...
        'Substituting black.\n'], index_or_name);
      RGB= [0 0 0];
   end
   return
end

error('Input argument has unexpected data type.');

end

