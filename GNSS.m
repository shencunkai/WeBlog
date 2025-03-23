data = rinexread("C:\Users\sck03\Desktop\BRDM00DLR_S_20250640000_01D_MN.rnx");
GPS = data.GPS; 
GM = 3.986005*10^14;
N = size(GPS,1);
XsYsZs = zeros(N,3);
for i = 1:N

n0 = sqrt(GM)/(GPS.sqrtA(i))^3;
n = n0 + GPS.Delta_n(i);


t = datetime(GPS.Time(i),'Format','yyyy-MM-dd HH:mm:ss');
UT = hour(t) + (minute(t)/60) + ((second(t))/3600);
juliandate = fix(362.25*year(t)) + fix(30.6001*(month(t)+1))+day(t)+(UT/24)+1720981.5; 
t = fix((juliandate -2444244.5)/7);
t = (juliandate -2444244.5-t*7)*24*3600;
t = t-(GPS.SVClockBias(i)+GPS.SVClockDrift(i)*(t-GPS.Toe(i))+GPS.SVClockDriftRate(i)*(t-GPS.Toe(i))*(t-GPS.Toe(i)));
tk = t - GPS.Toe(i);
if tk > 302400
    tk = tk-604800;
elseif tk <-302400
    tk = tk+604800;
end

M = GPS.M0(i) + n*tk;
E = M;
E0 = 0;
while abs(E0 - E)>=10^-12
    E0 = E;
    E  =  M + GPS.Eccentricity(i)*sin(E);
    
end

V = atan2(sqrt(1 - GPS.Eccentricity(i)^2)*sin(E),(cos(E) - GPS.Eccentricity(i)));
U = V + GPS.omega(i);

delta_u = GPS.Cuc(i)*cos(2*U) + GPS.Cus(i)*sin(2*U);
delta_r = GPS.Crc(i)*cos(2*U) + GPS.Crs(i)*sin(2*U);
delta_i = GPS.Cic(i)*cos(2*U) + GPS.Cis(i)*sin(2*U);

uk = U + delta_u;
rk = ((GPS.sqrtA(i))^2)*(1 - GPS.Eccentricity(i)*E) + delta_r;
ik = GPS.i0(i) + tk*GPS.IDOT(i) + delta_i;

xk = rk*cos(uk);
yk = rk*sin(uk);

We = 7.29211567*10^-5;
L = GPS.OMEGA0(i) + GPS.OMEGA_DOT(i)*tk - We*(tk + GPS.Toe(i));

X = xk*cos(L) - yk*cos(ik)*sin(L);
Y = xk*sin(L) + yk*cos(ik)*cos(L);
Zk = yk*sin(ik);

XsYsZs(i,1:3) = [X Y Zk];

end

