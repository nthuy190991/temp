function create_shp_from_snake_result(version_name,verbose)

reso_grid=[0.15 0.15];

%% Load snake results
current_dir=pwd;
dir_result='3279_to_3404/';%strcat(version_name,'/');

cd(dir_result)
list=dir(strcat('*',version_name,'.mat'));

cd(current_dir)

for i=1:length(list)
        
    filename=list(i).name;
    load(strcat(dir_result,filename),'patch_corner','snake_results','snake_mask')
    
    str=strsplit(filename,'17_');
    tile_no=str{2}(1:4);
    
    strx=tile_no(1:2);
    stry=tile_no(3:4);
    
    xlims=[2e5+str2double(strx)*1e3, 2e5+str2double(strx)*1e3+1e3];
    ylims=[51e5+str2double(stry)*1e3, 51e5+str2double(stry)*1e3+1e3];
    
    if verbose
        figure, hold on
        for i=1:length(snake_results)
            s=snake_results{i,:};
            if ~isempty(s)
                x=(s(:,2)+patch_corner(i,1))*reso_grid(1)+xlims(1);
                y=ylims(2)-(s(:,1)+patch_corner(i,2))*reso_grid(2);
                plot(x,y,'b','LineWidth',2)
                % text(mean(x),mean(y),num2str(i))
            end
        end
    end
    
    
    %% Export results
    str=strsplit(filename,'_LiDAR2017_');
    tile_no=str{2}(1:4);
    
    [~,~]=create_shp_from_geo_coor(strcat(dir_result,'res_BE_LiDAR2017_',tile_no,'_',version_name),...
        xlims,ylims,reso_grid,1,verbose);
    
end
end