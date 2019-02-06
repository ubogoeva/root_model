close all;
clear all;

%%

display('process started');
delete(gcp);
N = 3;
myCluster=parcluster('local');
myCluster.NumWorkers=N; 
parpool(myCluster,N);

%%
load ('variables_for_model.mat');

%%
%save('variables_for_model.mat', 'center_of_cell', '-append')

%% adjustment for real cell sizes
x=zeros(size(labeled));
for i=1:size(labeled,1)
    for j=1:size(labeled,2)
        if labeled(i,j)==2
            x(i,j)=1;
        end;
    end;
end;
x1=sum(ctranspose(x));

i=1;
z=1;
while z~=0
    z=x1(i);
    i=i+1;
end;

diameter=length(labeled(labeled(i,:)~=1));

k=100/diameter;
cellwall_length=cellwall_length*k;
cell_volume=cell_volume*(k^2);

%%

time=20000;
dt=10;
printTime=0:dt:time;

%%

%%this is true from PIN_graph.png
l=2;
q1PIN1=0.3/l;
q2PIN1=0.5/l;
SPIN1=2;
q3PIN1=0.63/l;
hPIN1=4;

q1PIN2=0.13/l;
q2PIN2=0.8/l;
SPIN2=2;
q3PIN2=0.19/l;
hPIN2=4;

q1PIN3=0.88/l;
q2PIN3=1.09/l;
SPIN3=6;
q3PIN3=1.1/l;
hPIN3=10;

q1PIN7=0.55/l;
q2PIN7=0.69/l;
SPIN7=10;
q3PIN7=0.6/l;
hPIN7=12;

q1PIN4=0.75/l;
q2PIN4=0.75/l;
SPIN4=6;
%%dont run it
%%

Auxin_ini=ones(n,1)*0.01;
Auxin_ini(1)=0;
Auxin_ini(2)=0;
%Auxin_ini=Auxin_end;
%Auxin_ini=zeros(n,1);
%Auxin_ini(30)=1;
%Auxin_ini(223)=1;
%Auxin_ini(368)=1;


KsA=0.01;
KdA=0.005;
%KdA=0;
D=0.08;
Kalpha=1.4;
%Kalpha=0;


PIN1_ini=zeros(n,1);

KsPIN1=1000;
%q1PIN1=0.1;
q1PIN1=0.15;
%q2PIN1=0.4;
q2PIN1=0.25;
SPIN1=2;

KdPIN1=1000;
%q3PIN1=0.4;
q3PIN1=0.315;
hPIN1=4;

PIN2_ini=zeros(n,1);

KsPIN2=1000;
%q1PIN2=0.05;
L2=2;
q1PIN2=0.065*L2;
%q2PIN2=0.8;
q2PIN2=0.4*L2;
SPIN2=2;
KdPIN2=1000;
%q3PIN2=0.15;
q3PIN2=0.095*L2;
hPIN2=4;
%hPIN2=6;

PIN3_ini=zeros(n,1);
%PIN3(:,1)=PIN3_end;

KsPIN3=1000;
%q1PIN3=0.9;
q1PIN3=0.375;
%q2PIN3=1.1;
q2PIN3=0.375;
%SPIN3=8;
SPIN3=6;

KdPIN3=1000;


k0PIN1=zeros(n,n);
k0PIN2=zeros(n,n);
k0PIN3=zeros(n,n);
%%

load('cellwalltype.mat');

stele=[169	174	175	178	179	180	181	182	183	185	186	187	188	189	190	192	193	194	195	196	197	198	199	200	202	203	206	207	209	210	211	214	215	216	217	218	219	220	221	222 ...
    223	224	225	226	228	229	230	231	233	234	235	236	237	238	239	240	241	242	243	244	245	246	247	248	249	250	251	258	259	260	261	263	264	265	266	267	269	270	271	273	274 ...
    275	276	277	278	279	280	281	282	283	284	285	286	287	288	289	290	291	292	293	295	296	297	298	299	300	301	302	303	304	305	306	307	308	311	313	314	315	316	317	318	319 ...
    320	322	323	324	325	327	330	331	332	333	335	338	340	341	342	343	347	348	349	369	371	372	373	374	376 ];

pericycle=[107	125	127	130	132	133	134	135	136	137	138	139	140	141	142	144	145	146	149	150	151	152	153	155	156	157	158	159	160	161	163	164	166	167	168	170	191	201	212 ...
    232	306	326	336	337	339	345	346	351	352	353	354	355	356	357	358	359	360	361	362	363	364	366	368	377	378	382	384	386	387	392	405	406	407	408	409	411 ];

endodermis=[84	87	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	108	110	111	112	113	114	115	116	117	118	119	120	126	128	129	131	143	154	162	171 ...
    204	334	344	367	375	379	381	385	388	389	390	391	394	396	397	398	399	400	401	402	403	404	416	417	419	420	424	427	431	432	435	436	437	438	440	444	446 ];

cortex=[40	41	43	45	46	47	49	50	52	53	54	55	57	58	59	60	61	62	63	64	65	66	67	68	69	72	73	75	76	77	78	80	81	85	86	89	121	147	176	370 ...
    383	393	395	413	415	421	422	423	425	426	428	429	430	433	439	441	442	443	445	447	448	449	450	451	452	453	454	455	456	457	458	461	462	464	465	466	467];

epidermis=[6 7	8	9	10	11	12	13	14	15	16	17	18	20	21	22	23	24	25	26	27	28	29	30	32	33	35	39	42	48	56	83	109	148	184	418	434	460	463	469 ...
    471	472	473	475	476	478	479	480	481	482	483	484	485	487	488	489	490	491	492	493	494	495	496	497	498	499	500	501	502	504	505	506];

QC = [272 309];

initiale= [232  251 288 308 343 366 350 365 310 268 208 213];

epidermis_cortex_upper=[6 460 40 370];
%determine type of cells

cells_with_PIN1=[stele, pericycle, endodermis]; %split PIN1 and PIN2
cells_with_PIN2=[epidermis, cortex];
cells_without_PIN1_PIN2=[];
for i=3:n %first 2 cells are not real cells
    if isempty(intersect(i,[cells_with_PIN1, cells_with_PIN2]))==1
        cells_without_PIN1_PIN2=[cells_without_PIN1_PIN2,i]; %determine cells without PIN1 and PIN2
    end;
end;

cells_with_1=zeros(n,n);
cells_with_2=zeros(n,n);
cells_with_3=zeros(n,n);
for i=3:n
    for j=2:n
        if cellwalltype_in_out(i,j)==1
            cells_with_1(i,j)=cellwall_length(i,j);
        elseif cellwalltype_in_out(i,j)==2
            cells_with_2(i,j)=cellwall_length(i,j);
        elseif cellwalltype_in_out(i,j)==3
            cells_with_3(i,j)=cellwall_length(i,j);
        end;
    end;
end;
sum_cells_with_1=sum(ctranspose(cells_with_1)); %we count cells in upper boundary
sum_cells_with_2=sum(ctranspose(cells_with_2)); %we count cells in down boundary
sum_cells_with_3=sum(ctranspose(cells_with_3)); %inner boundary

for i=stele
    for j=3:n
        if cellwalltype_in_out(i,j)==2
            k0PIN1(i,j)=1*cells_with_2(i,j)/sum_cells_with_2(i);
        else
            k0PIN1(i,j)=0;
        end;
    end;
end;

for i=[pericycle, endodermis]
    for j=3:n
        if cellwalltype_in_out(i,j)==2 
            k0PIN1(i,j)=0.8*cells_with_2(i,j)/sum_cells_with_2(i);
        elseif cellwalltype_in_out(i,j)==3
            k0PIN1(i,j)=0.2*cells_with_3(i,j)/sum_cells_with_3(i);
        else
            k0PIN1(i,j)=0;
        end;
    end;
end;


for i=cortex
    for j=2:n %j=3:n %for closed boundary
        if center_of_cell(i, 2) >= center_of_cell(55, 2)
            if cellwalltype_in_out(i,j)==1
                k0PIN2(i,j)=0.8*cells_with_1(i,j)/sum_cells_with_1(i);
            elseif cellwalltype_in_out(i,j)==3
                k0PIN2(i,j)=0.2*cells_with_3(i,j)/sum_cells_with_3(i);
            else
                k0PIN2(i,j)=0;
            end;
        else
             if cellwalltype_in_out(i,j)==2
                k0PIN2(i,j)=0.8*cells_with_2(i,j)/sum_cells_with_2(i);
            elseif cellwalltype_in_out(i,j)==3
                k0PIN2(i,j)=0.2*cells_with_3(i,j)/sum_cells_with_3(i);
            else
                k0PIN2(i,j)=0;
            end;
        end;
    end;
end;

for i=epidermis %in epidermis we determine PIN2 distribution and polarity (/Volume of cell)
    for j=2:n
        if cellwalltype_in_out(i,j)==1
            k0PIN2(i,j)=0.8*cells_with_1(i,j)/sum_cells_with_1(i);
        elseif cellwalltype_in_out(i,j)==3
            k0PIN2(i,j)=0.2*cells_with_3(i,j)/sum_cells_with_3(i);
        else
            k0PIN2(i,j)=0;
        end;
    end;
end;


cellwall_for_diffusion=cellwall_length;
for i=endodermis
    for j=1:n
        if cellwalltype_in_out(i,j)==3 || cellwalltype_in_out(i,j)==4
        cellwall_for_diffusion(i,j)=0;
        cellwall_for_diffusion(j,i)=0;
        end;
    end;
end;
%sum_cellwall_for_diffusion = sum(cellwall_length)-cellwall_length(1,:)-cellwall_length(2,:); %closed upper boundary
%sum_cellwall_for_diffusion = sum(cellwall_length)-cellwall_length(1,:); %opened upper boundary
sum_cellwall_for_diffusion = sum(ctranspose(cellwall_for_diffusion))-cellwall_length(1,:); %opened upper boundary
%sum_cellwall_for_diffusion = sum(ctranspose(cellwall_for_diffusion))-cellwall_length(1,:)-cellwall_length(2,:); %closed upper boundary


sum_cellwall_for_nonpolarised = sum(ctranspose(cellwall_length))-cellwall_length(1,:); 


for i=cells_with_PIN1
    for j=2:n
        k0PIN3(i,j)=k0PIN1(i,j);
    end;
end;

for i=cells_without_PIN1_PIN2 %%for cells that don't contain PIN1, PIN2
    for j=2:n
        k0PIN3(i,j)= cellwall_length(i,j)/sum_cellwall_for_nonpolarised(i);
    end;
end;


for i=cells_with_PIN2
    for j=2:n
        k0PIN3(i,j)=k0PIN2(i,j);
    end;
end;

cellwall_length_for_alpha=0;
cells_with_alpha=[107, 169, 206, 233, 269, 289, 306]; %%it is a flow of auxin from upper cell
for i=cells_with_alpha
    cellwall_length_for_alpha=cellwall_length_for_alpha+cellwall_length(i,2);
end;

neighbour=cell(n,1);
for i=3:n 
    for j=2:n  %i=3:n for closed boundary
        if cellwall_length(i,j)~=0
            neighbour{i}=[neighbour{i},j];
        end;
    end;
end;

cellswithKsA=[]; %here we determine TAA1-mediated synthesis of auxin
for i=[cortex, epidermis]
    if center_of_cell(i, 2) >= center_of_cell(55, 2)
        cellswithKsA=[cellswithKsA,i];
    end
end
for i=QC
    cellswithKsA=[cellswithKsA,i];
    for j=neighbour{i}
        if cellwalltype_in_out(i,j)==2
        cellswithKsA=[cellswithKsA,j];
        end;
    end;
end;




%%
%for Kalpha=0.8:0.1:2
Kalpha=1;

y=[Auxin_ini;PIN1_ini;PIN2_ini;PIN3_ini];
options=[];
tic;
[T, res] = ode15s(@dif,printTime,y,options, n,...
    Kalpha,cellwall_length_for_alpha,cells_with_alpha,...
    neighbour,cellwall_length,cell_volume,... 
    cells_with_PIN1,cells_with_PIN2,cells_without_PIN1_PIN2,...
    KsA, cellswithKsA, KdA,D,...
    KsPIN1,q1PIN1,q2PIN1,SPIN1,KdPIN1,q3PIN1,hPIN1,k0PIN1,...
    KsPIN2,q1PIN2,q2PIN2,SPIN2,KdPIN2,q3PIN2,hPIN2,k0PIN2, ...
    KsPIN3,q1PIN3,q2PIN3,SPIN3,KdPIN3,k0PIN3);
toc;

%%
save('auxin_counted_Kalfa_001_PIN_adjust_PIN2_double.mat', 'res')
%load('auxin_counted.mat')
%%
%save(['D:\Realistic Root Model\New root PIN in all cells\PINs_alpha_variation\A_maxPINintensity_1_time', num2str(time), '_alpha' num2str(Kalpha), '_Ks', num2str(KsA), '_Kd', num2str(KdA), '_D', num2str(D), '_PIN1_l1_', num2str(l1),...
    %'_PIN2_l2_', num2str(l2), '_PIN3_l3_', num2str(l3),  '_PIN3_l4_', num2str(l4), '_PIN7_l7_', num2str(l7), '.mat'], 'res');


tic;

load('MyColormaps1','mycmap');

map1=colormap(mycmap);

T_pic=size(res,1);
image1=picture_function_color(ctranspose(res(:,1:n)), T_pic,labeled,map1);
image2=picture_function_color(ctranspose(res(:,n+1:2*n)), T_pic, labeled,map1);
image3=picture_function_color(ctranspose(res(:,2*n+1:3*n)), T_pic, labeled,map1);
image4=picture_function_color(ctranspose(res(:,3*n+1:4*n)), T_pic, labeled,map1);

fig=figure ('visible','on');

subplot(1, 4, 1), subimage(image1);
colorbar();
colormap(map1);
caxis([0,max(max(res(:,1:n)))]);
axis off;
title ('Auxin');
subplot(1, 4, 2), subimage(image2);
colorbar();
colormap(map1);
caxis([0,max(max(res(:,n+1:2*n)))]);
axis off;
title ('PIN1');
subplot(1, 4, 3), subimage(image3);
colorbar();
colormap(map1);
caxis([0,max(max(res(:,2*n+1:3*n)))]);
axis off;
title ('PIN2');
subplot(1, 4, 4), subimage(image4);
colorbar();
colormap(map1);
caxis([0,max(max(res(:,3*n+1:4*n)))]);
axis off;
title ('PIN3');


%print(fig,'-dpng','-r650',['A_maxPINintensity_1_time', num2str(time), '_alpha' num2str(Kalpha), '_Ks', num2str(KsA), '_Kd', num2str(KdA), '_D', num2str(D),...
    %'_PIN1_l1_', num2str(l1), '_PIN2_l2_', num2str(l2), '_PIN3_l3_', num2str(l3), '.png']);

toc;

%end;

%%

load('MyColormaps1','mycmap');

map1=colormap(mycmap);

T_pic=size(res,1);
image1=picture_function_color(ctranspose(res(:,1:n)), T_pic,labeled,map1);
image2=picture_function_color(ctranspose(res(:,n+1:2*n)), T_pic, labeled,map1);
image3=picture_function_color(ctranspose(res(:,2*n+1:3*n)), T_pic, labeled,map1);
image4=picture_function_color(ctranspose(res(:,3*n+1:4*n)), T_pic, labeled,map1);
image5=picture_function_color(ctranspose(res(:,4*n+1:5*n)), T_pic, labeled,map1);
image6=picture_function_color(ctranspose(res(:,5*n+1:6*n)), T_pic, labeled,map1);

fig=figure ('visible','on');

subplot(2,3,1), subimage(image1);
colorbar();
colormap(map1);
caxis([0,max(max(res(:,1:n)))]);
axis off;
title ('Auxin');
subplot(2,3,2), subimage(image2);
colorbar();
colormap(map1);
caxis([0,max(max(res(:,n+1:2*n)))]);
axis off;
title ('PIN1');
subplot(2,3,3), subimage(image3);
colorbar();
colormap(map1);
caxis([0,max(max(res(:,2*n+1:3*n)))]);
axis off;
title ('PIN2');
subplot(2,3,4), subimage(image4);
colorbar();
colormap(map1);
caxis([0,max(max(res(:,3*n+1:4*n)))]);
axis off;
title ('PIN3');
subplot(2,3,5), subimage(image5);
colorbar();
colormap(map1);
caxis([0,max(max(res(:,4*n+1:5*n)))]);
axis off;
title ('PIN3');
subplot(2,3,6), subimage(image6);
colorbar();
colormap(map1);
caxis([0,max(max(res(:,5*n+1:6*n)))]);
axis off;
title ('PIN7');


%%
KsA_var=0.01:0.01:0.1;


%%
spmd
    ind=labindex;
KsA_var=0.01:0.01:0.1;
for percentIndx=ind:numlabs:size(KsA_var,2)
        KsA=KsA_var(percentIndx);
for l7=0.7:0.05:1.3
    for l3=0.7:0.05:1.3
        for l4=0.7:0.05:1.3
    
    q1PIN1=0.3/l1;
    q2PIN1=0.5/l1;
    SPIN1=2;
    q3PIN1=0.63/l1;
    hPIN1=4;

    q1PIN2=0.13/l2;
    q2PIN2=0.8/l2;
    SPIN2=2;
    q3PIN2=0.19/l2;
    hPIN2=4;

    q1PIN3=0.88/l3;
    q2PIN3=1.09/l3;
    SPIN3=6;
    q3PIN3=1.1/l3;
    hPIN3=10;

    q1PIN7=0.55/l7;
    q2PIN7=0.69/l7;
    SPIN7=10;
    q3PIN7=0.6/l7;
    hPIN7=12;

    q1PIN3=0.75/l4;
    q2PIN3=0.75/l4;
    SPIN3=6;


Y=[Auxin_ini;PIN1_ini;PIN2_ini;PIN3_ini;PIN3_ini;PIN7_ini];
options=[];
tic;
[T, res] = ode15s(@dif_tissues1_withKsAasTAA1_PIN3,printTime,Y,options, n,...
    Kalpha,cellwall_length_for_alpha,cells_with_alpha,...
    neighbour,cellwall_length,cell_volume,... 
    cells_with_PIN1,cells_with_PIN2,cells_without_PIN1_PIN2,...
    KsA, cellswithKsA, KdA,D,...
    KsPIN1,q1PIN1,q2PIN1,SPIN1,KdPIN1,q3PIN1,hPIN1,k0PIN1,...
    KsPIN2,q1PIN2,q2PIN2,SPIN2,KdPIN2,q3PIN2,hPIN2,k0PIN2, ...
    KsPIN3,q1PIN3,q2PIN3,SPIN3,KdPIN3,q3PIN3,hPIN3,k0PIN3,...
    KsPIN3,q1PIN3,q2PIN3,SPIN3,KdPIN3,k0PIN3,...
    KsPIN7,q1PIN7,q2PIN7,SPIN7,KdPIN7,q3PIN7,hPIN7,k0PIN7);
toc;

%save(['A_maxPINintensity_1_time', num2str(time), '_alpha' num2str(Kalpha), '_Ks', num2str(KsA), '_Kd', num2str(KdA), '_D', num2str(D), '_PIN1', num2str(q1PIN1), '..', num2str(q2PIN1), '..', num2str(SPIN1), '..',...
 %   num2str(q3PIN1), '..', num2str(hPIN1), '_PIN2', num2str(q1PIN2), '..', num2str(q2PIN2), '..', num2str(SPIN2), '..', num2str(q3PIN2), '..', num2str(hPIN2),...
  %  '_PIN3', num2str(q1PIN3), '..', num2str(q2PIN3), '..', num2str(SPIN3), '..', num2str(q3PIN3), '..', num2str(hPIN3), '_PIN3', num2str(q1PIN3), '..', num2str(q2PIN3), '..',...
   % num2str(SPIN3), '_PIN7', num2str(q1PIN7), '..', num2str(q2PIN7), '..', num2str(SPIN7), '..', num2str(q3PIN7), '..', num2str(hPIN7), '.mat'], 'res');

   %save(['D:\Realistic Root Model\New root PIN in all cells\PINs_1\A_maxPINintensity_1_time', num2str(time), '_alpha' num2str(Kalpha), '_Ks', num2str(KsA), '_Kd', num2str(KdA), '_D', num2str(D), '_PIN1_l1_', num2str(l1),...
    %'_PIN2_l2_', num2str(l2), '_PIN3_l3_', num2str(l3),  '_PIN3_l4_', num2str(l4), '_PIN7_l7_', num2str(l7), '.mat'], 'res');

 fun_to_save(['A_maxPINintensity_1_time', num2str(time), '_alpha' num2str(Kalpha), '_Ks', num2str(KsA), '_Kd', num2str(KdA), '_D', num2str(D), '_PIN1_l1_', num2str(l1),...
    '_PIN2_l2_', num2str(l2), '_PIN3_l3_', num2str(l3),  '_PIN3_l4_', num2str(l4), '_PIN7_l7_', num2str(l7), '.mat'], res)
   
 

tic;
T_pic=size(res,1);
image1=picture_function(ctranspose(res(:,1:n)), T_pic,labeled);
image2=picture_function(ctranspose(res(:,n+1:2*n)), T_pic, labeled);
image3=picture_function(ctranspose(res(:,2*n+1:3*n)), T_pic, labeled);
image4=picture_function(ctranspose(res(:,3*n+1:4*n)), T_pic, labeled);
image5=picture_function(ctranspose(res(:,4*n+1:5*n)), T_pic, labeled);
image6=picture_function(ctranspose(res(:,5*n+1:6*n)), T_pic, labeled);

fig=figure ('visible','off');
subplot(1,6,1), subimage(image1);
axis off;
title ('Auxin');
subplot(1,6,2), subimage(image2);
axis off;
title ('PIN1');
subplot(1,6,3), subimage(image3);
axis off;
title ('PIN2');
subplot(1,6,4), subimage(image4);
axis off;
title ('PIN3');
subplot(1,6,5), subimage(image5);
axis off;
title ('PIN3');
subplot(1,6,6), subimage(image6);
axis off;
title ('PIN7');


%print(fig,'-dpng','-r650',['D:\Realistic Root Model\New root PIN in all cells\A_maxPINintensity_1_time', num2str(time), '_alpha' num2str(Kalpha), '_Ks', num2str(KsA), '_Kd', num2str(KdA), '_D', num2str(D), '_PIN1', num2str(q1PIN1), '..', num2str(q2PIN1), '..', num2str(SPIN1), '..',...
 %   num2str(q3PIN1), '..', num2str(hPIN1), '_PIN2', num2str(q1PIN2), '..', num2str(q2PIN2), '..', num2str(SPIN2), '..', num2str(q3PIN2), '..', num2str(hPIN2),...
  %  '_PIN3', num2str(q1PIN3), '..', num2str(q2PIN3), '..', num2str(SPIN3), '..', num2str(q3PIN3), '..', num2str(hPIN3), '_PIN3', num2str(q1PIN3), '..', num2str(q2PIN3), '..',...
   % num2str(SPIN3), '_PIN7', num2str(q1PIN7), '..', num2str(q2PIN7), '..', num2str(SPIN7), '..', num2str(q3PIN7), '..', num2str(hPIN7),'.png']);

print(fig,'-dpng','-r650',['D:\Realistic Root Model\New root PIN in all cells\PINs_1\A_maxPINintensity_1_time', num2str(time), '_alpha' num2str(Kalpha), '_Ks', num2str(KsA), '_Kd', num2str(KdA), '_D', num2str(D),...
    '_PIN1_l1_', num2str(l1), '_PIN2_l2_', num2str(l2), '_PIN3_l3_', num2str(l3),  '_PIN3_l4_', num2str(l4), '_PIN7_l7_', num2str(l7), '.png']);

toc;


        end;
    end;
end;
end;
end;
%%
tic;
T_pic=2000;
image1=picture_function(ctranspose(res(:,1:n)), T_pic,labeled);
image2=picture_function(ctranspose(res(:,n+1:2*n)), T_pic, labeled);
image3=picture_function(ctranspose(res(:,2*n+1:3*n)), T_pic, labeled);
image4=picture_function(ctranspose(res(:,3*n+1:4*n)), T_pic, labeled);
image5=picture_function(ctranspose(res(:,4*n+1:5*n)), T_pic, labeled);
image6=picture_function(ctranspose(res(:,5*n+1:6*n)), T_pic, labeled);

subplot(1,6,1), subimage(image1);
axis off;
title ('Auxin');
subplot(1,6,2), subimage(image2);
axis off;
title ('PIN1');
subplot(1,6,3), subimage(image3);
axis off;
title ('PIN2');
subplot(1,6,4), subimage(image4);
axis off;
title ('PIN3');
subplot(1,6,5), subimage(image5);
axis off;
title ('PIN3');
subplot(1,6,6), subimage(image6);
axis off;
title ('PIN7');

toc;

%%
common_auxin2=zeros(size(printTime));

for i=1:size(printTime,2)
    for j=3:n
        common_auxin2(i)=common_auxin2(i)+res(i,j)*cell_volume(j);
    end;
end;
