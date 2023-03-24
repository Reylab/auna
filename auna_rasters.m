function auna_rasters(handles,ZonesList,ConceptsList,merge_concepts, filename)

[filepath, filename, ext] = fileparts(filename);

LineFormat = struct;
LineFormat.Color = 'blue';
LineFormat.LineWidth = 0.35;

for ch_n = 1: length(handles.par.channels)
    
ch = handles.par.channels(ch_n);
classes = handles.par.classes{ch_n};

if iscell(classes)
    classes = cell2mat(classes); %fix when indexing handles.par.classes output cell
end

if merge_concepts
    concepts2use = 1;
    labels = {strjoin(ConceptsList,', ')};
    merge_events = [];
    for conci_merge = 1:length(ConceptsList)
        current_concept = ConceptsList{conci_merge};
        merge_events =  [merge_events handles.concepts.(current_concept)];
    end
    
else
    concepts2use = ConceptsList;
    labels = ConceptsList;
end


for cl_n = 1:length(classes)
    
    if strcmp(classes,'mu')
        clname = sprintf('ch %d-mu',ch);
        index = find(handles.ch_data{ch_n}.sp_index)*1e3/handles.sr;
    else
        clname = sprintf('ch %d-cl %d',ch,classes(cl_n));
        index = find(handles.ch_data{ch_n}.index(cl_n,:))*1e3/handles.sr;
    end

    nrows = length(concepts2use);
    ncols = length(ZonesList);
    fig = figure('Visible','off');
    for zi = 1:length(ZonesList)
        current_zone = ZonesList{zi};

        time_pre_ms = handles.zones.(current_zone).pre_time;
        time_pos_ms = handles.zones.(current_zone).post_time;
        limits = handles.zones.(current_zone).limits; 
        if any(isnan(limits))
            error(['Zone: ' current_zone ' has some undefined limit'])
        end

        for conci = 1:length(concepts2use)
            if merge_concepts
                events = merge_events((merge_events>limits(1)) & (merge_events<limits(2)));
            else
                current_concept = ConceptsList{conci};
                events = handles.concepts.(current_concept); %in ms
                events =  events((events>limits(1)) & (events<limits(2)));
            end
            subplot_ax = subplot(nrows,ncols,zi+(conci-1)*ncols);


            %rasters

            hold on;        
            spikes1 = cell(length(events),1);
            for ei=1:length(events)
                spikes1{ei}=index((index>(events(ei)-time_pre_ms)) ...
                    & (index<(events(ei)+time_pos_ms))) - events(ei);
            end
            if numel(spikes1)>0 
                lst = numel(spikes1);
                if ~all(cellfun(@isempty,spikes1))
                    plotSpikeRaster(spikes1,'PlotType','vertline','LineFormat',LineFormat);
                end
                set(subplot_ax,'YLim',[0.5 lst+0.5])
            end  

            set(subplot_ax,'XLim',[-time_pre_ms time_pos_ms]);  
            set(subplot_ax,'FontSize',7);
            axis off
            title([current_zone ': ' labels{conci}])

            limy=ylim; %
            ylimits = (limy(2));% 
            set(subplot_ax,'YLim',[0 ylimits]);
            line([0 0],[0 ylimits],'linestyle',':')
            line([1000 1000],[0 ylimits],'linestyle',':')
        end

    end

    figurename = [filepath filesep filename clname ext];
    saveas(fig,figurename);
    close(fig)
    end
end
end