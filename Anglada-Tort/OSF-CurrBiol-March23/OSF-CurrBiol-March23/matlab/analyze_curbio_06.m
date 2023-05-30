function [data,raw_data,UN,NK,NI,unets]=analyze_curbio_06(EXPR)
EXPRS={'slider-2tones','btw-2tones','btw-3tones','wth-3tones','btw-5tones','btw-free','melcon-2tones','mem-3tones-10','mem-3tones-5','mem-3tones-match','btw-3tones-india','wth-3tones-india'};

if sum(strcmp(EXPRS,EXPR))~=1
    fprintf('NORI:wrong EXPR=%s could not find str 1\n',EXPR);
    EXPRS
end
assert(sum(strcmp(EXPRS,EXPR))==1);

fname_base='data/';


switch EXPR
    case 'btw-3tones'
        fname_seed='exp1_btw_3tones_590chains.csv';
        NI=2;
    case 'wth-3tones'
        fname_seed='exp11_wth_3tones_615chains.csv';
        NI=2;
    case 'btw-5tones'

        fname_seed='exp2_btw_5tones_159chains.csv'; NI=4;
    case 'btw-2tones'
        fname_seed='exp5_btw_2tones_398chains.csv';NI=1;

    case 'btw-free'
        fname_seed='exp3_btw_free.notes_216chains.csv'; NI=13;
    case 'btw-7tones'
        fname_seed='old/btw-7tones-b1.120chains.csv'; NI=6;
    case 'slider-2tones'
        fname_seed='exp4_btw_2tones_slider_369chains.csv'; NI=1;
    case 'melcon-2tones'
        fname_seed='exp6_melodic.pleasant_all.batches.csv';NI=1;

    case 'mem-3tones-10'
        fname_seed='exp8_btw_3tones_memory.10sec_240chains.csv';NI=2;
    case 'mem-3tones-5'
        fname_seed='exp7_btw_3tones_memory.5sec_240chains.csv';NI=2;


    case 'mem-3tones-match'
        fname_seed='exp9_btw_3tones_memory.control_240chains.csv';NI=2;

    case 'btw-3tones-india'
        fname_seed='exp10_btw_3tones_india_120chains.csv';NI=2;
    case 'wth-3tones-india'
        fname_seed='exp12_wth_3tones_india_223chains.csv';NI=2;
    otherwise
        fprintf('NORI:wrong EXPR=%s could not find str 2\n',EXPR);
        assert(1==0)
end
fname=sprintf('%s%s',fname_base,fname_seed);

fprintf('reading fname: %s...\n',fname)
t=readtable( fname);

for ll=1:size(t,1)
    if isnumeric(t.network_id(ll))
        t.network_id_str{ll}=sprintf('%d',t.network_id(ll));
    else
        t.network_id_str{ll}=sprintf('%s',t.network_id{ll});
    end
end

unets=unique(t.network_id_str);
UN=length(unets);

if strcmp(EXPR,'melcon-2tones')
    fprintf('reorganizing  data (working on melodic consnance experiment!)...\n')
    UN= size(t,1);
    NK=2;
    NI=1;
    unets=1:UN;
    data=nan(UN,NK,NI);

    data(:,1,1)=t.target_interval;
    data(:,2,1)=t.answer;

elseif strcmp(EXPR,'btw-free')

    NK=length(unique(t.degree));
    
           assert(min(t.degree)==0);
        assert(max(t.degree)==(NK-1));

    fprintf('reorganizing data...\n')
    data=nan(UN,NK,NI);
    for ll=1:size(t,1)
        raw_net=t.network_id_str(ll);
        raw_iter=t.degree(ll);
        net=find(strcmp(unets,raw_net),1);
        iter=raw_iter;
        vals=str2double(split(t.sung_intervals{ll},','));
        assert(~isnan(sum(vals)))
        for I=1:length(vals)
            data(net,iter+1,I)=vals(I);
        end

    end

        fprintf('done loading!\n')


    else
        NK=length(unique(t.degree));
        assert(min(t.degree)==0);
        assert(max(t.degree)==(NK-1));
        


        fprintf('reorganizing data...\n')
        data=nan(UN,NK,NI);
        for ll=1:size(t,1)
            raw_net=t.network_id_str(ll);
            raw_iter=t.degree(ll);
            net=find(strcmp(unets,raw_net),1);
            iter=raw_iter;

            for I=1:NI
              
                field=sprintf('sung_interval%d',I);
                if  ( (NI==1) && (~strcmp(EXPR,'btw-2tones')))
                    field=sprintf('interval');
                end
                if (strcmp(EXPR,'slider-2tones'))
                    field=sprintf('location');% even though thre are duplicates - they will be elimated
                end
                
                vals=getfield(t,field);
                
                data(net,iter+1,I)=vals(ll);

            end
            % check we are not missing intervals
            I=(NI+1);
            field=sprintf('target_interval%d',I);
            assert(~isfield(t,field));

        end
        fprintf('done loading!\n')

    end
    raw_data=t;

 