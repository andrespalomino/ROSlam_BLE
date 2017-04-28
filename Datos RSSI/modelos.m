% node8263

v_d_real = [1 2 3 4]*0.6;
dist_real = ones(11,4);
for i = 1:4
    dist_real(:,i) = dist_real(:,i)*i*0.6;
end  

prom_node_8263 = [sum(-66 -68 -70 -71 -71 -72)/6;
        sum(-78 -78 -85 -86 -86 -87 -87 -87 -85 -84 -84)/11;
        sum(-82 -83 -81 -82 -84 -83 -87 -85 -84)/9;
        sum(-84 -86 -86 -87 -86 -86 -87 -87)/8]
    
load('node_8263')

figure
set(gca,'FontSize',11); set(gcf,'Color','White');
plot(v_d_real,prom_node_8263,'LineWidth',1)
hold on
plot(dist_real,node8263,'-r','LineWidth',2)
xlabel('Distance (m)')
ylabel('Signal Strenght  (dB)')
grid on


% Modelo Distancia funcion RSSI node8263
% y = 0,0037e-0,073x

n_datos = 100;
v_dist_sim = linspace(-90,-60,n_datos);
for j = 1:length(v_dist_sim)
    dis_modelo(j) = 0.0037*exp(-0.073 * v_dist_sim(j)); 
end


figure
set(gca,'FontSize',11); set(gcf,'Color','White');
plot(prom_node_8263,v_d_real,'LineWidth',1)
hold on
plot(node8263,dist_real,'-r','LineWidth',2)
plot(v_dist_sim,dis_modelo,'LineWidth',1)
ylabel('Distance (m)')
xlabel('Signal Strenght  (dB)')
grid on

%%

% node11392
v_d_real = [1 2 3 4 5]*0.6;
dist_real = ones(10,5);
for i = 1:5
    dist_real(:,i) = dist_real(:,i)*i*0.6;
end  

prom_node_11392 = [-70.71428571	-76.3	-81.77777778	-86	-86.42857143];

load('node_11392')

figure
set(gca,'FontSize',11); set(gcf,'Color','White');
plot(v_d_real,prom_node_11392,'LineWidth',1)
hold on
plot(dist_real,node11392,'-r','LineWidth',2)
xlabel('Distance (m)')
ylabel('Signal Strenght  (dB)')
grid on

% Modelo Distancia funcion RSSI node11392
% y = 0,0008e-0,094x

n_datos = 100;
v_dist_sim = linspace(-90,-60,n_datos);
for j = 1:length(v_dist_sim)
    dis_modelo(j) = 0.0008*exp(-0.094 * v_dist_sim(j)); 
end

figure
set(gca,'FontSize',11); set(gcf,'Color','White');
plot(prom_node_11392,v_d_real,'LineWidth',1)
hold on
plot(node11392,dist_real,'-r','LineWidth',2)
plot(v_dist_sim,dis_modelo,'LineWidth',1)
ylabel('Distance (m)')
xlabel('Signal Strenght  (dB)')
grid on
