function Grav3D
%%% Code by Oksum E. Pham. L.T.,(2020) 
%%% Main GUI builder function
clc;clear all;clear global;
crt_Gui
save('storeF.mat')
end
function crt_Gui 
% create the main window and menu compenents
f=emptW;set(f,'Name','Grav3D')
crt_menu(f)
end

function START_func(~,~,modeP) %%% START an inverion or forward procedure
try
switch modeP
   case 1
   [r0,z0,~,~,~,~,n,data,dx,dy]=call_inputs(modeP); % read inputs from GUI
   [ny0,nx0]=size(data);    
   hwtb = waitbar(.2,'processing');
   %%% run forward procedure and memorize outputs
   gcalc_FRWf=FORWARD_func(data,nx0,ny0,dx,dy,r0,z0,n); 
   save('storeF.mat','gcalc_FRWf','-append')
   load('storeF.mat','x_frw','y_frw')
   waitbar(.6,hwtb);pause(.1);waitbar(1,hwtb);delete(hwtb);
   drawnow;mapviewer(x_frw,y_frw,gcalc_FRWf,'mGal','','ax_21');
   set(findobj(gcf,'Tag','tgmsg2'),'string','Message [3]: Calculation performed')
   case -1
   [r0,z0,WH,SH,criterio,mxit,n,data,dx,dy]=call_inputs(modeP); % read inputs from GUI
   [ny0,nx0]=size(data);
   %%% run inverse procedure and memorize outputs
   [zcalc_INVf,rms_stor]=INVERSION_func(data,nx0,ny0,dx,dy,r0,z0,WH,SH,criterio,mxit);
   gcalc_INVf=FORWARD_func(zcalc_INVf,nx0,ny0,dx,dy,r0,z0,n);
   gdiff_INVf=data-gcalc_INVf;
   save('storeF.mat','zcalc_INVf','gcalc_INVf','gdiff_INVf','rms_stor','-append')
   load('storeF.mat','x_inv','y_inv')
   drawnow;mapviewer(x_inv,y_inv,-zcalc_INVf,'km','Calculated Depth','ax_12')
   h=findobj(gcf,'Tag','tbut3');set(h,'value',1);
   set(findobj(gcf,'Tag','tgmsg'),'string','Message [3]: Calculation performed')
panel_jump(1,2)
end
catch 
    hwtb = waitbar(0,'Process Failed ...');pause(2);delete(hwtb); 
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Subroutine functions....
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r0,z0,WH,SH,criterio,mxit,n,data,dx,dy]=call_inputs(modeP)
%%% read input parameters and settings from GUI and data from temp. file 
r0=0;z0=0;WH=0;SH=0;criterio=0;mxit=0;
n=10;
switch modeP
    case 1
    hp=findobj(gcf,'Tag','pnl_dt2');he=findall(hp,'style','edit');
    z0=str2double(get(he(1),'string'));z0=abs(z0);
    r0=str2double(get(he(2),'string'));
    load('storeF.mat','data_frw','dx_frw','dy_frw');
    data=data_frw;dx=dx_frw;dy=dy_frw;
    case -1
    hp=findobj(gcf,'Tag','pnl_dt1');he=findall(hp,'style','edit');
    mxit=str2double(get(he(1),'string'));mxit=abs(mxit);
    criterio=str2double(get(he(2),'string'));criterio=abs(criterio);
    SH=str2double(get(he(3),'string'));
    WH=str2double(get(he(4),'string'));
    WHSH=unique([WH SH]);WH=WHSH(1);SH=WHSH(2);
    z0=str2double(get(he(5),'string'));z0=abs(z0);
    r0=str2double(get(he(6),'string'));
    load('storeF.mat','data_inv','dx_inv','dy_inv');
    data=data_inv;dx=dx_inv;dy=dy_inv;    
end
end

function [zcalc,rmsset]=INVERSION_func(data,nx0,ny0,dx,dy,r0,z0,WH,SH,criterio,mxit)
%%% INVERSION PROCEDURE
[g0,nxe,nye]=paddData(nx0,ny0,data);% enlarge data
nxm=2*nxe; nym=2*nye;    
k=getfreqs(nxm,nym,dx,dy); % obtain wavenumbers   
LPF=LP_filt(k,WH,SH); % filter design
%%% run inversion sheme
[zcalc,rmsset]=INV_sheme(g0,r0,z0,k,LPF,criterio,nxe,nye,nx0,ny0,mxit);
end

function gcalc=FORWARD_func(data,nx0,ny0,dx,dy,r0,z0,n)
%%% FORWARD PROCEDURE
z=data-z0; % extract mean depth
[ze,nxe,nye]=paddData(nx0,ny0,z); % enlarge data
nxm=2*nxe; nym=2*nye;
k=getfreqs(nxm,nym,dx,dy); % obtain wavenumbers
%%% run forward gravity calculation
gcalc=FW_PRKR3D(ze,z0,r0,nxe,nye,nx0,ny0,n,k); % calculate gravity
end

function [matrix,nx,ny]=paddData(nx,ny,matrix)
% PADDING DATA
matrix(1,nx+floor(nx/2))=0;
matrix(ny+floor(ny/2),1)=0;
matrix=rot90(rot90(matrix));
matrix(1,nx+2*floor(nx/2))=0;
matrix(ny+2*floor(ny/2),1)=0;
matrix=rot90(rot90(matrix));
if (mod(nx,2)~=0) nx=nx-1; matrix(:,end)=[]; end
if (mod(ny,2)~=0) ny=ny-1; matrix(end,:)=[]; end
end

function [k]=getfreqs(nxm,nym,dx,dy) 
% WAVENUMBERS
dkx= 2.*pi./((nxm-1).*dx);
dky= 2.*pi./((nym-1).*dy);
nyqx= (nxm/2)+1;
nyqy= (nym/2)+1; 
[kx,ky]=meshgrid([(0:nxm/2) (nyqx+1:nxm)-(nxm+1)].*dkx,...
    [((0:nym/2)) (nyqy+1:nym)-(nym+1)].*dky);
k= sqrt(bsxfun(@plus,kx.^2,ky.^2));
k(1,1)=0.00000001;
end

function LPF=LP_filt(k,WH,SH)
%%% build filter
LPF=k.*0;
[nym,nxm]=size(k);     
k=k./(2*pi);  
for j=1:nym;
   for i=1:nxm;
      if k(j,i)<WH
      LPF(j,i)=1;  
elseif k(j,i)<SH
      LPF(j,i)=0.5.*(1+cos((((2*pi)*k(j,i))-(2*pi*WH))/(2*(SH-WH))));
else
LPF(j,i)=0;
      end
    end;
end;
end

function [z,r]=INV_sheme(g0,r0,z0,k,LPF,criterio,nx,ny,nx0,ny0,mxit)
%%% the iterative inversion sheme
hwtb = waitbar(0,'Iteration Starting, please wait..');
pause(1)
Fh=-fft2(g0)./(2.*pi.*6.67.*r0.*exp(-k.*z0));
Fh=Fh.*LPF;
Fh(1,1)=0;
h=real(ifft2(Fh)); 
h_old=h;
rms=1000;
iter=0;
m=2;
while rms>criterio
    Fh=Fh-((-k).^(m-1)).*fft2(h.^m)./factorial(m);
    Fh=Fh.*LPF;
    Fh(1,1)=0;
    h=real(ifft2(Fh));
    dh=h-h_old;
    dh2=dh.^2;
    rms=sqrt(sum(sum(dh2))./(numel(dh2)));
    iter=iter+1;
    h_old=h;
    r(iter)=rms;
    waitbar(iter/mxit,hwtb,sprintf('%12.0f',iter))
    if iter==mxit;break;end
    m=m+1;
end
z=h_old+z0;
z=z(ny/2+1:ny/2+ny0,nx/2+1:nx/2+nx0);
pause(.5);delete(hwtb);
end

function g=FW_PRKR3D(z,z0,r0,nx,ny,nx0,ny0,n,k)
%%% forward gravity calculation
hs=-2*pi*20/3.*r0.*exp(-abs(k).*z0);
tongF=0;
for m=1:n;
     tongF=tongF+((-abs(k)).^(m-1))./(factorial(m)).*fft2(z.^m);
end;
Fg=hs.*tongF;
g0=(ifft2(Fg));
g0(1,1)=0;
g=real(g0);
g=g(ny/2+1:ny/2+ny0,nx/2+1:nx/2+nx0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calc. functions
function [pwr,wn]= raps_data(T,dx)
% get radially averaged spectrum vs frequency 
[ny,nx]=size(T); 
nrow=2*floor(ny/2); ncol=2*floor(nx/2); 
maxdim=max([nrow ncol]);
np2=2^nextpow2(maxdim);  
rowdiff=round((np2-nrow)/2); 
coldiff=round((np2-ncol)/2);
T=T(1:nrow,1:ncol);
T=T-mean(T(:)); 
wf=tukeywin(nrow,.05)*tukeywin(ncol,.05)'; %truncation window 5% from edges
dw=T.*wf;
TT=zeros(np2); 
TT(rowdiff+1:rowdiff+nrow,coldiff+1:coldiff+ncol)=dw;  
spec =(abs(fftshift(fft2(TT))));
spec(spec<1.0e-13)=1.0e-13;
spec=log(spec);
% create cartesian coordinate system (unit) and
% shifting the origin to zero frequency.
[xo,yo]=meshgrid((1:np2)-np2/2,(1:np2)-np2/2); 
[~,L]=cart2pol(xo,yo);L=round(L);% calculate distances to the center
halfspec=(np2/2)-1; %considering the half of the spectrum
wn=(1:halfspec)/(np2*dx); % freq axis
for i=1:halfspec
pwr(i)=mean(spec(L==i));%find data of equal distances and get mean  
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot related functions
function mapviewer(x,y,matrix,unitt,tit,nam)
% plot of selected map
[ax,index]=axfinder(nam);
drawnow
set(gcf,'CurrentAxes',ax(index))
contourf(x,y,matrix,18);shading flat;rotate3d off;axis equal;axis tight
set(gca,'FontSize',10,'FontWeight','bold')
h=colorbar('eastoutside');
title(h,unitt,'FontWeight','bold');
set(h,'FontSize',10,'FontWeight','bold')
xlabel('X (km)');ylabel('Y (km)');title(tit)
xlim([min(min(x)) max(max(x))])
ylim([min(min(y)) max(max(y))])
box on
set(ax(index),'Tag',nam)
end

function  [ax,index]=axfinder(nam)
%axis selection
ax=findall(gcf,'type','axes');
tgs=get(ax,'Tag');
index=strcmpi(tgs,nam);index=find(index==1);
axes(ax(index));
end

function plotout_mapinv(src,event)
warning off
try
switch get(event.NewValue,'String')
    case 'Zcalc'
    load('storeF.mat','x_inv','y_inv','zcalc_INVf')    
    drawnow;mapviewer(x_inv,y_inv,-zcalc_INVf,'km','Calculated Depth','ax_12')    
    case 'RMS'
    load('storeF.mat','rms_stor')    
    rms_plot(rms_stor,'ax_12')  
    case 'Gcalc'
    load('storeF.mat','x_inv','y_inv','gcalc_INVf')    
    drawnow;mapviewer(x_inv,y_inv,gcalc_INVf,'mGal','Calculated Gravity','ax_12')
    case 'Gobs-Gcalc'
    load('storeF.mat','x_inv','y_inv','gdiff_INVf')    
    drawnow;mapviewer(x_inv,y_inv,gdiff_INVf,'mGal','Gravity Difference','ax_12')
    case 'Gobs'
    load('storeF.mat','x_inv','y_inv','data_inv')     
    drawnow;mapviewer(x_inv,y_inv,data_inv,'mGal','Observed Gravity','ax_12')
end
catch
end
end

function rms_plot(rms_it,nam) % plot the rms graphic 
[ax,index]=axfinder(nam);
plot(rms_it,'-ko','MarkerFaceColor','r','MarkerSize',5);
set(gca,'FontSize',10,'FontWeight','bold')
sent1=['RMS(1)=' num2str(rms_it(1))];
sent2=['RMS(end)=' num2str(rms_it(end))];
title ([sent1 '    ' sent2])
xlabel('Iteration number');ylabel('RMS (km)');
xlim([1 numel(rms_it)])
grid on
pbaspect([1 .9 1]);
set(ax(index),'Tag',nam)
end

function rapsplot(wn,pwr)
%%% plot of raps
[ax,index]=axfinder('ax_11');
drawnow
set(gcf,'CurrentAxes',ax(index))
plot(wn,pwr,'k','LineWidth',2);
axis on
pbaspect([1 1 1])
grid on
xlabel('k/2*pi');ylabel('Log(P)');
pbaspect([2 1 1]);
set(ax(index),'Tag','ax_11')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GUI related functions 
function f=emptW
f=figure('MenuBar','none','NumberTitle','off','DockControls','off',...
'Color','w','units','normalized','outerposition',[0 .03 1 .97]);
menu1 = uimenu('Parent',f,'Label','MODE'); 
uimenu('Parent',menu1,'Label','INVERSION','CallBack',{@selectmode,-1});
uimenu('Parent',menu1,'Label','FORWARD','CallBack',{@selectmode,1});
uimenu('Parent',f,'Label','INVERSION','ForegroundColor','b',...
    'enable','off','Tag','tgmode');
end
function crt_menu(f)
bg=uibuttongroup(f,'units','normalized','Position',[0 0.96 1 0.04],...
'BackgroundColor','w','BorderType','etchedout','Tag','BG1',...
'FontWeight','bold','SelectionChangeFcn',@panelup);
uicontrol(bg,'Style','togglebutton','String','DATA & SETTINGS','FontWeight','bold',...
'units','normalized','Position',[0 0 .2 1],'Tag','tbut1');
uicontrol(bg,'Style','togglebutton','String','OUTPUT',...
'FontWeight','bold','units','normalized', 'Position',[0.2 0 .2 1],'Tag','tbut2');
uicontrol(bg,'style','text','units','normalized','Position',[.4 .0 .6 1],...
'String','Message [1]: Import new gravity grid ...','BackGroundColor','w','ForegroundColor','r',...
'FontWeight','bold','Tag','tgmsg')
panel_output_menu(f)
panel_data_menu(f)
end

function panel_data_menu(f)
hp0=uipanel(f,'units','normalized','Position',[0 0 1 0.96],...
'BackgroundColor','w','Tag','panel_data');
hp12=uipanel(hp0,'units','normalized','Position',[0 0 1 1],...
'BackgroundColor','w','Tag','pnl_dt2','visible','off');
hp11=uipanel(hp0,'units','normalized','Position',[0 0 1 1],...
'BackgroundColor','w','Tag','pnl_dt1');
pnl_dt1_menu(hp11)%%INV panel Data
pnl_dt2_menu(hp12)%%FRW panel Data
end
function pnl_dt2_menu(pn) 
uicontrol(pn,'style','pushbutton','units','normalized','Position',[0 .93 .2 .06],...
'String','Import New Depth Grid','CallBack',{@call_data,1})
uicontrol(pn,'style','text','units','normalized','Position',[0 .9 .2 .03],...
'String','File Name','BackgroundColor','w','Tag','tgsorc2')
hp11dc=uipanel(pn,'units','normalized','Position',[.2 .93 .12 .07],...
'BackgroundColor','w','Title','Density Contrast (gr/cc)','FontWeight','bold',...
'TitlePosition','lefttop','ForegroundColor','b');
uicontrol(hp11dc,'style','edit','units','normalized','Position',[0.1 .2 .8 .6],...
'FontWeight','bold','BackgroundColor','w','string','0.4')
hp11md=uipanel(pn,'units','normalized','Position',[.32 .93 .12 .07],...
'BackgroundColor','w','Title','Mean Depth(km)','FontWeight','bold',...
'TitlePosition','lefttop','ForegroundColor','b');
uicontrol(hp11md,'style','edit','units','normalized','Position',[0.1 .2 .8 .6],...
'FontWeight','bold','BackgroundColor','w','string','20')
hpfg1=uipanel(pn,'units','normalized','Position',[0 .05 .5 .85],...
'BackgroundColor','w','Title','Depth Model','FontWeight','bold',...
'TitlePosition','centertop','ForegroundColor','b');
hpfg2=uipanel(pn,'units','normalized','Position',[.5 .05 .5 .85],...
'BackgroundColor','w','Title','Gravity Model','FontWeight','bold',...
'TitlePosition','centertop','ForegroundColor','b');
axes('Parent',hpfg1,'units','normalized','Position',...
[0.1 0.1 0.8 0.8],'Tag','ax_20');axis off
axes('Parent',hpfg2,'units','normalized','Position',...
[0.1 0.1 0.8 0.8],'Tag','ax_21');axis off
uicontrol(pn,'style','pushbutton','units','normalized','Position',[.2 .9 .24 .03],...
'FontWeight','bold','string','START FORWARD','CallBack',{@START_func,1})
%%%%%% export items
uicontrol(pn,'style','text','units','normalized','Position',[0.5 0.92 .5 .04],...
'BackgroundColor','w','ForegroundColor','r','Tag','tgmsg2','string',...
'Message [1]: Import new depth grid ...','FontSize',10)
c1 = uicontextmenu;c2 = uicontextmenu;
uicontrol(hpfg1,'UIContextMenu',c1,'string','Export (right click)',...
'units','normalized','ForeGroundColor','k',...
'Position',[0,.965,.2,.035]);
uimenu('Parent',c1,'Label','Image','Separator','on',...
'Callback',{@frw_outs,'data_frw','x_frw','y_frw','*.png','km','Depth Model'});
uicontrol(hpfg2,'UIContextMenu',c2,'string','Export (right click)',...
'units','normalized','ForeGroundColor','k',...
'Position',[0,.965,.2,.035]);
uimenu('Parent',c2,'Label','Image','Separator','on',...
'Callback',{@frw_outs,'gcalc_FRWf','x_frw','y_frw','*.png','mGal','Gravity Model'});
uimenu('Parent',c2,'Label','Data','Separator','on',...
'Callback',{@frw_outs,'gcalc_FRWf','x_frw','y_frw','*.grd'});
%%%%%%%%%%%%%%%%%%%
end

function pnl_dt1_menu(pn) 
uicontrol(pn,'style','pushbutton','units','normalized','Position',[0 .93 .2 .06],...
'String','Import New Gravity Grid','CallBack',{@call_data,-1})
uicontrol(pn,'style','text','units','normalized','Position',[0 .9 .2 .03],...
'String','File Name','BackgroundColor','w','Tag','tgsorc1')
hp11dc=uipanel(pn,'units','normalized','Position',[.2 .93 .12 .07],...
'BackgroundColor','w','Title','Density Contrast (gr/cc)','FontWeight','bold',...
'TitlePosition','lefttop','ForegroundColor','b');
uicontrol(hp11dc,'style','edit','units','normalized','Position',[0.1 .2 .8 .6],...
'FontWeight','bold','BackgroundColor','w','string','0.4')
hp11md=uipanel(pn,'units','normalized','Position',[.32 .93 .12 .07],...
'BackgroundColor','w','Title','Mean Depth(km)','FontWeight','bold',...
'TitlePosition','lefttop','ForegroundColor','b');
uicontrol(hp11md,'style','edit','units','normalized','Position',[0.1 .2 .8 .6],...
'FontWeight','bold','BackgroundColor','w','string','20')
hp11wh=uipanel(pn,'units','normalized','Position',[.44 .93 .12 .07],...
'BackgroundColor','w','Title','Filter P. WH','FontWeight','bold',...
'TitlePosition','lefttop','ForegroundColor','b');
uicontrol(hp11wh,'style','edit','units','normalized','Position',[0.1 .2 .8 .6],...
'FontWeight','bold','BackgroundColor','w','string','0.02','CallBack',{@setlinSHWH,'L1'})
hp11sh=uipanel(pn,'units','normalized','Position',[.56 .93 .12 .07],...
'BackgroundColor','w','Title','Filter P. SH','FontWeight','bold',...
'TitlePosition','lefttop','ForegroundColor','b');
uicontrol(hp11sh,'style','edit','units','normalized','Position',[0.1 .2 .8 .6],...
'FontWeight','bold','BackgroundColor','w','string','0.04','CallBack',{@setlinSHWH,'L2'})
hp11cr=uipanel(pn,'units','normalized','Position',[.68 .93 .12 .07],...
'BackgroundColor','w','Title','RMS Criterio','FontWeight','bold',...
'TitlePosition','lefttop','ForegroundColor','b');
uicontrol(hp11cr,'style','edit','units','normalized','Position',[0.1 .2 .8 .6],...
'FontWeight','bold','BackgroundColor','w','string','0.0001')
hp11mx=uipanel(pn,'units','normalized','Position',[.8 .93 .12 .07],...
'BackgroundColor','w','Title','Max. Iteration','FontWeight','bold',...
'TitlePosition','lefttop','ForegroundColor','b');
uicontrol(hp11mx,'style','edit','units','normalized','Position',[0.1 .2 .8 .6],...
'FontWeight','bold','BackgroundColor','w','string','100')
hpfg1=uipanel(pn,'units','normalized','Position',[0 .05 .3 .85],...
'BackgroundColor','w','Title','Observed Gravity','FontWeight','bold',...
'TitlePosition','centertop','ForegroundColor','b');
hpfg2=uipanel(pn,'units','normalized','Position',[.3 .05 .7 .80],...
'BackgroundColor','w','Title','Raps','FontWeight','bold',...
'TitlePosition','centertop','ForegroundColor','b');
axes('Parent',hpfg1,'units','normalized','Position',...
[0.1 0.1 0.8 0.8],'Tag','ax_10');axis off
axes('Parent',hpfg2,'units','normalized','Position',...
[0.1 0.1 0.8 0.8],'Tag','ax_11');axis off
uicontrol(pn,'style','pushbutton','units','normalized','Position',[.2 .9 .72 .03],...
'FontWeight','bold','string','START INVERSION','CallBack',{@START_func,-1})
end

function panel_output_menu(f)
%%%INV
hp0=uipanel(f,'units','normalized','Position',[0 0 1 0.96],...
'BackgroundColor','w','Tag','panel_output','visible','off');
bg1=uibuttongroup(hp0,'units','normalized','Position',[0 0.96 1 0.04],...
'BackgroundColor','y','BorderType','etchedout',...
'FontWeight','bold','SelectionChangeFcn',@plotout_mapinv);
uicontrol(bg1,'Style','togglebutton','String','Zcalc','FontWeight','bold',...
'units','normalized','Position',[0 0 .2 1],'Tag','tbut3','ForegroundColor','k');
uicontrol(bg1,'Style','togglebutton','String','RMS','ForegroundColor','k',...
'FontWeight','bold','units','normalized', 'Position',[0.2 0 .2 1],'Tag','tbut4');
uicontrol(bg1,'Style','togglebutton','String','Gcalc','ForegroundColor','k',...
'FontWeight','bold','units','normalized', 'Position',[0.4 0 .2 1],'Tag','tbut5');
uicontrol(bg1,'Style','togglebutton','String','Gobs-Gcalc','ForegroundColor','k',...
'FontWeight','bold','units','normalized', 'Position',[0.6 0 .2 1],'Tag','tbut6');
uicontrol(bg1,'Style','togglebutton','String','Gobs',...
'FontWeight','bold','units','normalized', 'Position',[0.8 0 .2 1],'Tag','tbut7');
axes('Parent',hp0,'units','normalized','Position',...
     [0.15 0.1 0.7 0.7],'Tag','ax_12');
axis off
%%%export items
hpexp=uipanel(hp0,'units','normalized','Position',[0 .89 1 0.06],...
'BackgroundColor','w','Tag','pnl_exp');
uicontrol(hpexp,'style','checkbox','units','normalized',...
'BackgroundColor','w','String','Zcalc','value',1,'Position',[0 0 .2 1])
uicontrol(hpexp,'style','checkbox','units','normalized',...
'BackgroundColor','w','String','RMS','value',1,'Position',[.2 0 .2 1])
uicontrol(hpexp,'style','checkbox','units','normalized',...
'BackgroundColor','w','String','Gcalc','value',1,'Position',[.4 0 .2 1])
uicontrol(hpexp,'style','checkbox','units','normalized',...
'BackgroundColor','w','String','Gobs-Gcalc','value',1,'Position',[.6 0 .2 1])
uicontrol(hpexp,'style','checkbox','units','normalized',...
'BackgroundColor','w','String','Gobs','value',1,'Position',[.8 0 .2 1])
c= uicontextmenu;
uicontrol(hp0,'UIContextMenu',c,'string','Export (right click)',...
'units','normalized','ForeGroundColor','b',...
'Position',[0 .86 1 0.03]);
uimenu('Parent',c,'Label','Save as Image(s)','Separator','on',...
'Callback',{@inv_outs,'*.png'});
uimenu('Parent',c,'Label','Save as numeric Data(s)','Separator','on',...
'Callback',{@inv_outs,'*.grd'});
end

function panelup(~,~)
h=get(findobj(gcf,'Tag','tbut1'),'value');
if h==1;
    drawnow;uistack(findobj(gcf,'Tag','panel_data'),'top');
    drawnow;set(findobj(gcf,'Tag','panel_output'),'visible','off')
    drawnow;set(findobj(gcf,'Tag','panel_data'),'visible','on')
end
if h==0;
    drawnow;uistack(findobj(gcf,'Tag','panel_output'),'top');
    drawnow;set(findobj(gcf,'Tag','panel_data'),'visible','off')
    drawnow;set(findobj(gcf,'Tag','panel_output'),'visible','on')
end
drawnow;refresh
end

function selectmode(~,~,modeP)
switch modeP
    case -1
    set(findobj(gcf,'Tag','tgmode'),'Label','INVERSION')
    set(findobj(gcf,'Tag','tbut2'),'visible','on')
    set(findobj(gcf,'Tag','tgmsg'),'visible','on')
    uistack(findobj(gcf,'Tag','pnl_dt1'),'top');
    drawnow;set(findobj(gcf,'Tag','pnl_dt1'),'visible','on')
    drawnow;set(findobj(gcf,'Tag','pnl_dt2'),'visible','off')
    case 1
    set(findobj(gcf,'Tag','tgmode'),'Label','FORWARD')
    set(findobj(gcf,'Tag','tgmsg'),'visible','off')
    set(findobj(gcf,'Tag','tbut2'),'visible','off')
    uistack(findobj(gcf,'Tag','pnl_dt2'),'top');
    drawnow;set(findobj(gcf,'Tag','pnl_dt2'),'visible','on')
    drawnow;set(findobj(gcf,'Tag','pnl_dt1'),'visible','off')
    panel_jump(2,1)
end
drawnow;refresh
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Export functions
function frw_outs(~,~,dt,dtx,dty,typ,unitt,tit) %% export forward output data to file
try
[filename, pathname] = uiputfile(typ,'Set filename');
sorcf=[pathname filename];
T=load('storeF.mat',dt,dtx,dty);
D =getfield(T, dt);
X=getfield(T,dtx);
Y=getfield(T,dty);
if strcmp(unitt,'km');D=-D;end
 switch typ
    case '*.png'
    img_export_map(X,Y,D,unitt,tit,sorcf)
    case '*.grd'
    xmin=min(X(:));xmax=max(X(:));ymin=min(Y(:));ymax=max(Y(:));
    grdout(D,xmin,xmax,ymin,ymax,sorcf)
 end
catch
end
end

function img_export_map(x,y,matrix,unitt,tit,sorcf)
tfg=figure('MenuBar','none','NumberTitle','off','Resize','off',...
'Color','w','units','normalized','outerposition',[0 0 1 1],...
'DockControls','off','visible','off');
contourf(x,y,matrix,18);shading flat;axis equal;axis tight
title(tit,'FontSize',16,'FontWeight','normal')
xlabel('X (km)','FontSize',16,'FontWeight','normal');
ylabel('Y (km)','FontSize',16,'FontWeight','normal');
xlim([min(min(x)) max(max(x))])
ylim([min(min(y)) max(max(y))])
set(gca,'FontSize',16,'FontWeight','normal')
h=colorbar('eastoutside');
title(h,unitt,'FontSize',16,'FontWeight','normal');
set(h,'FontSize',16,'FontWeight','normal');
print(tfg,sorcf, '-dpng', '-r300');
delete(tfg);
end

function inv_outs(~,~,typ) % export selected inversion outputs to file
try
fil=get(findobj(gcf,'Tag','tgsorc1'),'string');[~,fil]=fileparts(fil);
[filename, pathname] = uiputfile(typ,'Export',[fil '_outs']);
sorcf=[pathname filename];[~,sorc]=fileparts(sorcf);
h=get(findobj(gcf,'Tag','pnl_exp'),'children');%%checkboxes in pnl
r=flipud(cell2mat(get(h,'value')));%% values of checkboxes from left to right
dt=[{'zcalc_INVf'},{'rms_stor'},{'gcalc_INVf'},{'gdiff_INVf'},{'data_inv'}];
unitt=[{'km'},{'RMS(km)'},{'mGal'},{'mGal'},{'mGal'}];
tit=[{'Calculated Depth'},{'RMS'},{'Calculated Gravity'},...
    {'Gravity Difference'},{'Observed Gravity'}];
ext=[{'_zcalc'},{'_rms'},{'_gcalc'},{'_gdiff'},{'_gobs'}];
index=find(r==0);
dt(index)=[];unitt(index)=[];tit(index)=[];ext(index)=[];% eleminate none selected data
if numel(index)==5;wtb=msgbox('NO OUTPUT(s) SELECTED');pause(2);delete(wtb);end
load('storeF.mat','x_inv','y_inv');
xmin=min(x_inv(:));xmax=max(x_inv(:));
ymin=min(y_inv(:));ymax=max(y_inv(:));
hwtb = waitbar(0,'Exporting');
for i=1:numel(dt)
    waitbar(i/numel(dt),hwtb,['Exporting:' tit(i)])
     T=load('storeF.mat',cell2mat(dt(i)));
     D =getfield(T, cell2mat(dt(i)));
     if strcmp(unitt(i),'km');D=-D;end
 if strcmp(typ,'*.png')
 if strcmp(dt(i),'rms_stor');img_export_vctr(D,[pathname sorc cell2mat(ext(i)) '.png']);continue;end    
 img_export_map(x_inv,y_inv,D,cell2mat(unitt(i)),cell2mat(tit(i)),[pathname sorc cell2mat(ext(i)) '.png'])
 else
 if strcmp(dt(i),'rms_stor'); rmsvec=[(1:numel(D))' D'];
 save([pathname sorc cell2mat(ext(i)) '.dat'],'rmsvec','-ascii');continue;end    
 grdout(D,xmin,xmax,ymin,ymax,[pathname sorc cell2mat(ext(i)) '.grd'])
 end
 end
waitbar(i/numel(dt),hwtb,'Exporting Completed');pause(1.5);delete(hwtb)
catch
end
end

function img_export_vctr (rms_it,sorcf)
%%% export vector data
tfg=figure('MenuBar','none','NumberTitle','off','Resize','off',...
'Color','w','units','normalized','outerposition',[0 0 1 1],...
'DockControls','off','visible','off');
drawnow;plot(rms_it,'-ko','MarkerFaceColor','r','MarkerSize',5);
xlabel('Iteration number','FontSize',16);ylabel('RMS (km)','FontSize',16);
xlim([1 numel(rms_it)])
grid on
pbaspect([1 .8 1]);
sent1=['RMS(1)=' num2str(rms_it(1))];
sent2=['RMS(end)=' num2str(rms_it(end))];
title ([sent1 '    ' sent2])
drawnow;set(gca,'FontSize',16,'FontWeight','normal')
drawnow;print(tfg,sorcf, '-dpng', '-r300');
delete(tfg);
end
function grdout(matrix,xmin,xmax,ymin,ymax,namefile)
%%%%%%%%%%% Function for output of a GRID
%Get grid dimensions
aux=size(matrix);
nx=aux(2);ny=aux(1);
grdfile=fopen(namefile,'w');                % Open file
fprintf(grdfile,'%c','DSAA');               % Header code
fprintf(grdfile,'\n %i %i',[nx ny]);        % Grid size
fprintf(grdfile,'\n %f %f',[xmin xmax]);    % X limits
fprintf(grdfile,'\n %f %f',[ymin ymax]);    % Y limits
fprintf(grdfile,'\n %f %f',[min(min(matrix)) max(max(matrix))]); % Z limits
fprintf(grdfile,'\n');
for jj=1:ny                                 % Write matrix
    for ii=1:nx
        fprintf(grdfile,'%g %c',matrix(jj,ii),' ');
    end
    fprintf(grdfile,'\n');
end
fclose(grdfile);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Import Data Functions
function call_data(~,~,modeP) 
%interactively import gravity or depth data
clc
drawnow
[filename, pathname] = uigetfile('*.grd', 'Import Golden Software Binary/Text grid (*.grd)');
if ischar([pathname filename])
switch modeP
    case -1
[data_inv,x_inv,y_inv,nx_inv,ny_inv,~,~,~,~,dx_inv,dy_inv]=gridform([pathname filename]);
    errcode=error_checker(data_inv,dx_inv,dy_inv,nx_inv,ny_inv,-1);
    if errcode>0;return;end
    [pwr,wn]= raps_data(data_inv,dx_inv);
    save('storeF.mat','data_inv','x_inv','y_inv','dx_inv','dy_inv','-append');
    drawnow;rapsplot(wn,pwr);axis on
    drawnow;mapviewer(x_inv,y_inv,data_inv,'mGal','','ax_10');axis off
    set(findobj(gcf,'Tag','tgsorc1'),'string',filename)
    [~,~]=axfinder('ax_11');axis on;refresh
    STR='Message [2]: Loading performed ... Set input parameters and Start INVERSION';
    set(findobj(gcf,'Tag','tgmsg'),'string',STR)
    case 1
 [data_frw,x_frw,y_frw,nx_frw,ny_frw,~,~,~,~,dx_frw,dy_frw]=gridform([pathname filename]);
    errcode=error_checker(data_frw,dx_frw,dy_frw,nx_frw,ny_frw,1);
    if errcode>0;return;end
    save('storeF.mat','data_frw','x_frw','y_frw','dx_frw','dy_frw','-append');
    mapviewer(x_frw,y_frw,-(abs(data_frw)),'km','','ax_20')
    set(findobj(gcf,'Tag','tgsorc2'),'string',filename)
    STR='Message [2]: Loading performed ... Set input parameters and Start FORWARD';
    set(findobj(gcf,'Tag','tgmsg2'),'string',STR)
end
end
end

function [matrix,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=gridform(k);
%%grid data loader
fidc=fopen(k);header= fread(fidc,4,'*char' )';fclose(fidc);
c1=strcmp(header,'DSAA');c2=strcmp(header,'DSRB');
sumc=sum([c1 c2]);
if sumc>0
switch c1
    case 1
[matrix,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd6txt(k); %format surfer6 text
    case 0
[matrix,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd7bin(k); %surfer7 binary       
end
else
matrix=0;x=0;y=0;nx=0;ny=0;xmin=0;xmax=0;ymin=0;ymax=0;dx=0;dy=0;     
end
end

function [T,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy]=lodgrd6txt(k)
%%%%%reader, Surfer 6 text grid(*.grd)
surfergrd=fopen(k,'r'); % Open *.grid file
dsaa=fgetl(surfergrd);  % Header
% Get the map dimension [NX: East NY: North];
datasize=str2num(fgetl(surfergrd)); nx=datasize(1); ny=datasize(2);
% Map limits: xmin, xmax, ymin ymax
xcoor=str2num(fgetl(surfergrd)); xmin=xcoor(1); xmax=xcoor(2);
ycoor=str2num(fgetl(surfergrd)); ymin=ycoor(1); ymax=ycoor(2);
% check intervals in x and y direction 
dx=(xmax-xmin)/(nx-1);dx=abs(dx);
dy=(ymax-ymin)/(ny-1);dy=abs(dy);
% data limits
anom=str2num(fgetl(surfergrd)); t0min=anom(1); t0max=anom(2);
% data matrix 
[T,numb] = fscanf(surfergrd, '%f', [nx,ny]);
T=T'; % Traspose matrix
fclose(surfergrd);
% map coordinate matrix
[x,y]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);
end

function [T,x,y,nx,ny,xmin,xmax,ymin,ymax,dx,dy] = lodgrd7bin(filename)
%reader, Surfer 7 Binary grid
fid= fopen(filename);
fread(fid,4,'*char' )';
fread(fid,1,'uint32');fread(fid,1,'uint32');
fread(fid,4,'*char' )';fread(fid,1,'uint32');
ny= fread(fid,1,'uint32'); nx= fread(fid,1,'uint32');
xmin= fread(fid,1,'double'); ymin= fread(fid,1,'double');
dx= fread(fid,1,'double'); dy= fread(fid,1,'double');
fread(fid,1,'double');fread(fid,1,'double');
fread(fid,1,'double');
parm= fread(fid,1,'double');
fread(fid,4,'*char' )';
nn= fread(fid,1,'uint32');
if ny*nx ~= nn/8 ; error('error') ;end
T= nan(nx,ny);
T(1:end) = fread(fid,numel(T),'double');
T=T';
fclose(fid);
T(T==parm) = nan;
xv = xmin + (0:nx-1)*dx;
yv = ymin + (0:ny-1)*dy;
[x,y]=meshgrid(xv,yv);
xmax=xv(end);
ymax=yv(end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% check and controls
function errcode=error_checker(data,dx,dy,nx,ny,modeP)
%%% error message displayer
errcode=0;
if dx~=dy;errcode=1;STR='Message [error]: equal spaced grid dx=dy is required !';end
if any(isnan(data(:)));errcode=2;STR='Message [error]: Blanked grid not supported !';end
if nx==0;errcode=3;STR='Message [error]: File format not supported !';end
if ny==0;errcode=3;STR='Message [error]: File format not supported !';end
if modeP==-1;tgmbx='tgmsg';else tgmbx='tgmsg2';end
if errcode>0;drawnow;set(findobj(gcf,'Tag',tgmbx),'string',STR);end
end

function panel_jump(n1,n2)
warning off
h0=findobj(gcf,'Tag',['tbut' num2str(n1)]);h2=findobj(gcf,'Tag',['tbut' num2str(n2)]);
bt=findobj(gcf,'Tag','BG1');
eventdata.OldValue=h0;
eventdata.NewValue=h2;
eventdata.Source=bt;
eventdata.EventName='SelectionChanged';
set(bt,'SelectedObject',h2)
panelup(bt,eventdata);
end

function get_SHWH(~,~,srclabel,L)
[~,~]=axfinder('ax_11');refresh
set (gcf, 'WindowButtonMotionFcn',{@mouseMove,srclabel});
set (gcf, 'WindowButtonDownFcn', {@mouseLog,srclabel,L});
set(gcf,'Pointer','fullcrosshair');
end

function mouseMove(~,~,srclabel)
C = get (gca, 'CurrentPoint');
title([srclabel ': ' num2str(C(1,1))])
end

function mouseLog(~,~,srclabel,L)
try
C = get (gca, 'CurrentPoint');
hChildren = get(gca,'Children');
Y=get(hChildren(end),'YData');
drawnow
switch srclabel
    case 'Raps-WH'
    set(findall(findall(gcf,'type','uipanel','Title','Filter P. WH'),'style','edit'),...
    'string',num2str(C(1,1)))
    delete(findobj(gcf,'Tag',L))
    line([C(1,1) C(1,1)],[min(Y) max(Y)],'Tag',L,'Color','g','LineWidth',3,...
        'ButtonDownFcn',{@get_SHWH,'Raps-WH','L1'})
    case 'Raps-SH'
    set(findall(findall(gcf,'type','uipanel','Title','Filter P. SH'),'style','edit'),...
    'string',num2str(C(1,1)))
    delete(findobj(gcf,'Tag',L))
    line([C(1,1) C(1,1)],[min(Y) max(Y)],'Tag',L,'Color','r','LineWidth',3,...
        'ButtonDownFcn',{@get_SHWH,'Raps-SH','L2'})
end
catch
end
set (gcf, 'WindowButtonMotionFcn', '');
set (gcf, 'WindowButtonDownFcn', '');
set(gcf,'Pointer','arrow');
refresh
end

function setlinSHWH(src,~,L)
[~,~]=axfinder('ax_11');refresh
try
Lclr='g';Tstr='Raps-WH: ';srcL='Raps-WH';
if strcmp(L,'L2');Lclr='r';Tstr='Raps-SH: ';srcL='Raps-SH';end
X=str2double(get(src,'string'));
hChildren = get(gca,'Children');
Y=get(hChildren(end),'YData');
delete(findobj(gcf,'Tag',L))
line([X X],[min(Y) max(Y)],'Tag',L,'Color',Lclr,'LineWidth',3,...
    'ButtonDownFcn',{@get_SHWH,srcL,L})
title([Tstr num2str(X)])
catch
end
end


