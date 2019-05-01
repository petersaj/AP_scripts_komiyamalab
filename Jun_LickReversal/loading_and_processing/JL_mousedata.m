
function mice = JL_mousedata()

n = 1;
mice(n).name = 'JL041';
mice(n).gender = 'f';
mice(n).virus = 'GCaMP5g';
mice(n).location = 'PMM';
mice(n).switch = 'y';
mice(n).comment = 'lots of contamination';
mice(n).interneurons = 1:36;
mice(n).scope = 1;

n = 2;
mice(n).name = 'JL042';
mice(n).gender = 'f';
mice(n).virus = 'GCaMP5g';
mice(n).location = 'PMM';
mice(n).switch = 'y';
mice(n).comment = 'lots of contamination';
mice(n).interneurons = 1:32;
mice(n).scope = 1;

n = 3;
mice(n).name = 'JL043';
mice(n).gender = 'f';
mice(n).virus = 'GCaMP5g';
mice(n).location = 'PMM';
mice(n).switch = 'y';
mice(n).comment = 'ok';
mice(n).interneurons = 1:38;
mice(n).scope = 1;

n = 4;
mice(n).name = 'JL048';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP5g';
mice(n).location = 'PMM';
mice(n).switch = 'n';
mice(n).comment = 'lots of contamination';
mice(n).interneurons = 1:40;
mice(n).scope = 1;

%   n = 5;
%     mice(n).name = 'JL051';
%     mice(n).gender = 'm';
%     mice(n).virus = 'GCaMP6f';
%     mice(n).location = 'ALM';
%     mice(n).switch = 'only once';
%     mice(n).comment = 'cells got very sick early on';
%     mice(n).interneurons = 1:36;

n = 5;
mice(n).name = 'JL053';
mice(n).gender = 'f';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'ALM';
mice(n).switch = 'y';
mice(n).comment = 'only a small region of infection';
mice(n).interneurons = [1:28 65];
mice(n).scope = 1;

n = 6;
mice(n).name = 'JL055';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP5g';
mice(n).location = 'ALM';
mice(n).switch = 'n';
mice(n).comment = 'only a small region of infection';
mice(n).interneurons = [1:29];
mice(n).scope = 1;

n = 7;
mice(n).name = 'JL064';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP5g';
mice(n).location = 'ALM';
mice(n).switch = 'y';
mice(n).comment = '';
mice(n).interneurons = [1:37];
mice(n).scope = 1;

% from here on new mice
% scope 1 = MOM1
% scope 2 = Bscope

n = 8;
mice(n).name = 'JL068';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP5g';
mice(n).location = 'ALM';
mice(n).switch = 'n';
mice(n).comment = '';
mice(n).interneurons = [1:40];
mice(n).scope = 1;

n = 9;
mice(n).name = 'JL072';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'PMM';
mice(n).switch = 'y';
mice(n).comment = '';
mice(n).interneurons = [1:24];
mice(n).scope = 1;

n = 10;
mice(n).name = 'JL073';
mice(n).gender = 'f';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'PMM';
mice(n).switch = 'y';
mice(n).comment = '';
mice(n).interneurons = [1:16];
mice(n).scope = 1;

n = 11;
mice(n).name = 'JL074';
mice(n).gender = 'f';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'PMM';
mice(n).switch = 'y';
mice(n).comment = '';
mice(n).interneurons = [1:27];
mice(n).scope = 1;

n = 12;
mice(n).name = 'JL075';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'PMM';
mice(n).switch = 'n';
mice(n).comment = '';
mice(n).interneurons = [1:58];
mice(n).scope = 1;

n = 13;
mice(n).name = 'JL077';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'PMM';
mice(n).switch = 'n';
mice(n).comment = '';
mice(n).interneurons = [1:33];
mice(n).scope = 1;

n = 14;
mice(n).name = 'JL078';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'PMM';
mice(n).switch = 'y';
mice(n).comment = '';
mice(n).interneurons = [1:33];
mice(n).scope = 1;

n = 15;
mice(n).name = 'JL081';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'PMM';
mice(n).switch = 'y';
mice(n).comment = 'messed up the reversal, mouse did not learn properly after it. made roilabels for IN';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;

n = 16;
mice(n).name = 'JL082';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'PMM';
mice(n).switch = 'n';
mice(n).comment = 'made roilabels for IN';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;

n = 17;
mice(n).name = 'JL088';
mice(n).gender = 'f';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'ALM';
mice(n).switch = 'y';
mice(n).comment = 'made roilabels for IN';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;

% n = 18;
% mice(n).name = 'JL089';
% mice(n).gender = 'f';
% mice(n).virus = 'GCaMP6f';
% mice(n).location = 'ALM';
% mice(n).switch = 'n';
% mice(n).comment = 'stacks are messed up, cannot identify IN';
% mice(n).interneurons = find_interneurons(mice(n).name);

n = 18;
mice(n).name = 'JL093';
mice(n).gender = 'f';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'ALM';
mice(n).switch = 'y';
mice(n).comment = 'made roilabels for IN. Bad reversal animal, never managed to reverse more than once per day.';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;

n = 19;
mice(n).name = 'JL094';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'ALM';
mice(n).switch = 'y';
mice(n).comment = 'made roilabels for IN';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;

n = 20;
mice(n).name = 'JL100';
mice(n).gender = 'f';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'ALM';
mice(n).switch = 'n';
mice(n).comment = 'made roilabels for IN';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;

n = 21;
mice(n).name = 'JL101';
mice(n).gender = 'f';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'ALM';
mice(n).switch = 'n';
mice(n).comment = 'made roilabels for IN';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;

n = 22;
mice(n).name = 'JL092';
mice(n).gender = 'f';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'ALM';
mice(n).switch = 'n';
mice(n).comment = 'made roilabels for IN';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;

n = 23;
mice(n).name = 'JL098';
mice(n).gender = 'f';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'ALM';
mice(n).switch = 'n';
mice(n).comment = 'very weak labeling; made roilabels for IN';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;

n = 24;
mice(n).name = 'JL102';
mice(n).gender = '?';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'PMM';
mice(n).switch = 'n';
mice(n).comment = 'Entered by Andy: update later';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;


n = 25;
mice(n).name = 'JL113';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'PMM';
mice(n).switch = 'y';
mice(n).comment = '';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;


n = 26;
mice(n).name = 'JL116';
mice(n).gender = 'f';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'ALM';
mice(n).switch = 'y';
mice(n).comment = '';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;


n = 27;
mice(n).name = 'JL104';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'PMM';
mice(n).switch = 'n';
mice(n).comment = 'last two days are odor lick';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;


n = 28;
mice(n).name = 'JL107';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'PMM';
mice(n).switch = 'n';
mice(n).comment = '';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;


n = 29;
mice(n).name = 'JL119';
mice(n).gender = 'f';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'ALM';
mice(n).switch = 'n';
mice(n).comment = '';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;


n = 30;
mice(n).name = 'JL123';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'ALM';
mice(n).switch = 'y';
mice(n).comment = 'only one reversal, excluded last day';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;


n = 31;
mice(n).name = 'JL124';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'ALM';
mice(n).switch = 'y';
mice(n).comment = '';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;


n = 32;
mice(n).name = 'JL126';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'ALM';
mice(n).switch = 'y';
mice(n).comment = '';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;


n = 33;
mice(n).name = 'JL127';
mice(n).gender = 'm';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'ALM';
mice(n).switch = 'y';
mice(n).comment = '';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;


n = 34;
mice(n).name = 'JL129';
mice(n).gender = 'f';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'ALM';
mice(n).switch = 'n';
mice(n).comment = '';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;


n = 35;
mice(n).name = 'JL118';
mice(n).gender = 'f';
mice(n).virus = 'GCaMP6f';
mice(n).location = 'ALM';
mice(n).switch = 'n';
mice(n).comment = '';
mice(n).interneurons = find_interneurons(mice(n).name);
mice(n).scope = 2;

% this nested function reads the roilabel file to get the interneurons if
% they're not specified manually above
function interneurons = find_interneurons(mousename)
    mouse_dir = ['/usr/local/lab/People/Jun/Data/' mousename];
    roilabel_file = dir([mouse_dir filesep '*.roilabel']);
    
    if isempty(roilabel_file)
        warning(['No interneurons for ' mousename]);
        interneurons = [];
        return
    end
    
    load([mouse_dir filesep roilabel_file.name],'-mat','roi_labels');
    interneurons = find(cellfun(@(x) any(strcmp('IN',x)),roi_labels)');

end


end
