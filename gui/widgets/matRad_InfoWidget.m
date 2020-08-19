classdef matRad_InfoWidget < matRad_Widget
    
    properties
        
    end
    
    methods
        function this = matRad_InfoWidget(handleParent)
            if nargin < 1
                matRad_cfg = MatRad_Config.instance();
                handleParent = figure(...
                    'Units','normalized',...
                    'Position',[0.45 0.45 0.1 0.1],...
                    'Visible','on',...
                    'Color',matRad_cfg.gui.backgroundColor,... 
                    'IntegerHandle','off',...                    
                    'MenuBar','none',...
                    'Name','matRad Info',...
                    'NumberTitle','off',...
                    'HandleVisibility','callback',...
                    'Tag','figure1');
            end
            this = this@matRad_Widget(handleParent);
        end
    end
    
    methods (Access = protected)
        function this = createLayout(this)
            h94 = this.widgetHandle;
            matRad_cfg = MatRad_Config.instance();
            h95 = uicontrol(...
                'Parent',h94,...
                'Units','normalized',...
                'String','About',...
                'Position',[0.238095238095238 0.134831460674157 0.563492063492063 0.280898876404494],...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Callback',@(hObject,eventdata) btnAbout_Callback(this,hObject,eventdata),...
                'Tag','btnAbout',...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight','bold' );
            
            h96 = uicontrol(...
                'Parent',h94,...
                'Units','normalized',...
                ...%'String','v3.0.0',...
                'Style','text',...
                'Position',[0.227106227106227 0.752808988764045 0.523809523809524 0.191011235955056],...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Tag','text15',...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight','bold');
            
            h97 = uicontrol(...
                'Parent',h94,...
                'Units','normalized',...
                ...%'String','github.com/e0404/matRad',...
                'Style','text',...
                'Position',[0.0384615384615385 0.528089887640449 0.942307692307693 0.168539325842697],...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Tag','text31',...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight','bold' );
            
            this.createHandles();
            handles=this.handles;
            %Alter matRad Version string positioning
            vString = matRad_version();
            vPos = get(handles.text15,'Position');
            urlPos = get(handles.text31,'Position');
            btnPos = get(handles.btnAbout,'Position');
            
            %vPos([1 3]) = urlPos([1 3]);
            vPos([1 3]) = [0 1];
            vPos(4) = vPos(4)*1.25;
            btnPos(2) = 0.05;
            urlPos(2) = btnPos(2)+btnPos(4)+0.05;
            vPos(2) = urlPos(2) + urlPos(4) + 0.05;
            vPos(4) = 0.98 - vPos(2);
            
            set(handles.btnAbout,'Position',btnPos);
            set(handles.text31,'String','www.matRad.org','Position',urlPos,'Enable','inactive','ButtonDownFcn', @(~,~) web('www.matrad.org','-browser'));
            set(handles.text15,'String',vString,'Position',vPos);
            this.handles=handles;
        end
    end
    
    methods (Access = protected)
        
        function btnAbout_Callback(this, hObject, event)
            handles = this.handles;
            %msgbox({'https://github.com/e0404/matRad/' 'email: matrad@dkfz.de'},'About');
            
            matRad_cfg = MatRad_Config.instance();
            [~,matRadVer] = matRad_version;
            
            msg{1} = ['matRad ''' matRadVer.name '''']; %Name
            if matRad_cfg.eduMode
                msg{1} = [msg{1} ' Educational'];
            end
            msg{end+1} = sprintf('v%d.%d.%d',matRadVer.major,matRadVer.minor,matRadVer.patch); %Version Number
            if isdeployed
                msg{end+1} = 'Standalone Version';
            elseif ~isempty(matRadVer.branch) && ~isempty(matRadVer.commitID)
                msg{end+1} = sprintf('Git: Branch %s, commit %s',matRadVer.branch,matRadVer.commitID(1:8));
            end
            
            [env,envver]  = matRad_getEnvironment();
            msg{end+1} = sprintf('Environment: %s v%s %s',env,envver,version('-release'));
            
            msg{end+1} = 'Web: www.matrad.org';
            msg{end+1} = 'E-Mail: contact@matrad.org';
            
            msgbox(msg,'About matRad');
            
            
            this.handles = handles;
        end
    end
end
