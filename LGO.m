%% A LGO Special Case: Apogee is Maximal or Minimal Latitude

omE = 7.2921159e-5; % [s^-1]
mu = 398600.4; % [km^3/s^2] gravitational parameter 

% Data: rA, i

% phi = i; % [rad] angle from equator to parallel of projection of apogee
% rA % [km] radius apogee 
% i < 90 deg
% i = 63.4; % [deg] inclination for absence of perigee regression, i.e. minimum 
%                 characteristic velocity for station keeping

% om = +-90, 135

%Vom = omE * rA * cos(phi); % [km/s]
%V_A = sqrt((2*mu*rP)/(rA*(rA+rP))); %[km/s] = Vom

% rP = (omE^2*rA^4*cos(i)^2)/(2*mu-omE^2*rA^3*cos(i)^2) ; %[km]
% a = (rA+rP)/2;
a = (mu*rA)/(2*mu-omE^2*rA^3*cos(i*pi/180)^2); %[km]
% e = (eA-rP)/(rA+rP);
e = 1-(omE^2*rA^3*cos(i*pi/180)^2)/mu;

%% B LGO General Case

% Data:  rA, i, phi

alpha = acos(tan(phi)/tan(i));
gamma = acos(sqrt(1-sin(i)^2+tan(phi)^2*cos(i)^2));

a = (mu*rA*(1-sin(i)^2+tan(phi)^2*cos(i)^2))/(2*mu*(1-sin(i)^2+tan(phi)^2*cos(i)^2)+...
    -omE^2*rA^3*cos(phi)^2); % [km]
e = 1-(omE^2*rA^3*cos(phi)^2)/(mu*(1-sin(i)^2+tan(phi)^2*cos(i)^2));

uA = asin(sin(phi)/sin(i)) ; % argument of latitude of the apogee
om = asin(sin(phi)/sin(i)) - pi; 

%% C Geosynchronous LGO

% Data: i, (a, e from general or special case)
RE = 6378 ; %[km] Earth Equatorial Radius 
J2 = 0.0010827;
Tnod = 2*pi*a^(3/2)/sqrt(mu)*(1-((3/2)*J2*RE^2)/(a^2*(1-e^2)^2)*...
    ((10*cos(i)^2-2)/4-(3*cos(i)^2-1)/4*(1-e^2)^(1/2))); % Nodal period
dOM = -3*pi*J2*RE^2/(a^2*(1-e^2)^2)*cos(i); %dOM angular shift of the ascending node per one revolution

K = 1-((3/2)*J2*RE^2)/(a^2*(1-e^2)^2) * ...
((10*cos(i)^2-2)/4-(3*cos(i)^2-1)/4*(1-e^2)^(1/2)+ ...
(sqrt(mu)*cos(i))/(omE*a^(7/2)*(1-e^2)^2));

m_n = (omE*(a^(3/2)/sqrt(mu)*K))^(-1); % m/n repetition factor

% Tef = (2*pi+(m/n)*dOM)/omE; % 
Tef = (2*pi+(m_n)*dOM)/omE; % 

% Ttr = m*Tnod; % Track repetition period
% Ttr = n*Tef; % Efficient period

%% PLOT FIG6
clear
clc
close all
Rmean = 6371; % [km] Earth mean radius
mu = 398600.4; % [km^3/s^2] gravitational parameter 
omE = 7.2921159e-5; % [s^-1]
H = [0:1:60000]; % [km] altitude
r = H + Rmean; %[km] geocentric radius from Earth baricenter

Vc = sqrt(mu./r); % [km/s] circular velocity (first cosmic velocity) 

Hp = [500,5000,15000]; % [km] altitude at perigee
% Hp2 = 5000;
% Hp3 = 15000;
% rp1 = Hp1-Rmean;
rp = Hp+Rmean;
rA = r; % apogee radius is argument (variable) for V at apogee

i = [0,63.4]; % [deg] inclination

plot(H,Vc,'linewidth',2)
ax = gca;
hold(ax, 'on'); % Call hold FIRST 
% Then set colororder
colororder(ax,"gem12")
grid on
hold on
set(gcf,'NumberTitle', 'off', 'Name', 'Fig. 6  On geometric interpretation of locally geostationary orbit design')
title("Orbital Velocity vs Altitude")
xlabel("\boldmath $H_A [km]$",'Interpreter','latex')
ylabel("\boldmath $V [km/s]$",'Interpreter','latex')
ylim([0,8])
ax = gca;
ax.XAxis.Exponent = 0;

Va = zeros(length(rA),length(rp));
for k = 1:size(rp,2)
    % velocity at apogee, function of rA, fixed rp
    Va(:,k) = sqrt((2*mu*rp(k))./(rA.*(rA+rp(k)))); % [km/s] from eq.3
    Hindex(:,k) = Vc>=Va(:,k)';
    plot(H(Hindex(:,k)),Va(Hindex(:,k),k),'linewidth',2)
end
clear Hindex

Vi = zeros(length(rA),length(i));
Hindex(:,1) = and(H>=15000,H<=45000);
Hindex(:,2) = H>=19000;
for j = 1:size(i,2)
    % velocity of the end point of the geocentric radius vector, fuction of rA, 
    % fixed inclination
    Vi(:,j) = omE*rA*cos(i(j)*pi/180); % [km/s] from eq.2 where phi = i
    Hi = H';
    plot(Hi(Hindex(:,j)),Vi(Hindex(:,j),j),'linewidth',2)
end

legend('V_c(H_c)','V_A(H_A/H_P=500km)','V_A(H_A/H_P=5000km)','V_A(H_A/H_P=15000km)','V_i(H_A/i=0)','V_i(H_A/i=63.4\circ)')

fig = gcf;
fig.Position = [360.2, 288.2, 660.8, 664.8];

% Special LGO case: find one LGO orbit for a m/n class 
% rA will be computed for each m/n ratio, in the special LGO case
clear
clc
omE = 7.2921159e-5; % [s^-1]
mu = 398600.4; % [km^3/s^2] gravitational parameter 
Rmean = 6371; % [km] Earth mean radius
RE = 6378 ; %[km] Earth Equatorial Radius 
J2 = 0.0010827;
% Data: m/n, i -> plot Ha

i = [0,63.4]*pi/180; % [rad]
H = [0:1:60000]; % [km] altitude
r = H + Rmean; %[km] radius from Earth baricenter
rA = r;
yyaxis right
ylabel('\boldmath $m/n$','Interpreter','latex','FontWeight','bold')
ylim([0,8])
for j = 1:length(i)
    for z = 1:length(H)
        a = (mu*rA(z))/(2*mu-omE^2*rA(z)^3*cos(i(j))^2); %[km]
        e = 1-(omE^2*rA(z)^3*cos(i(j))^2)/mu;

        K = 1-((3/2)*J2*RE^2)/(a^2*(2-e^2)^2) * ...
            ((10*cos(i(j))^2-2)/4-(3*cos(i(j))^2-1)/4*(1-e^2)^(1/2)+ ...
            (sqrt(mu)*cos(i(j)))/(omE*a^(7/2)*(1-e^2)^2));

        m_n(z,j) = (omE*(a^(3/2)/sqrt(mu)*K))^(-1); % m/n repetition factor
    end
    rP(:,j) = (omE^2*rA.^4*cos(i(j))^2)./(2*mu-omE^2*rA.^3*cos(i(j))^2) ; %[km]
    index = and(and(m_n(:,j)>0,imag(m_n(:,j))==0),rA'>=rP(:,j));
    plot(H(index),m_n(index,j),'LineStyle','-.')
end
title("Orbital Velocity and Repetition Factor vs Altitude")
legend('$V_c(H_c)$','$V_A(H_A/H_P=500km)$','$V_A(H_A/H_P=5000km)$','$V_A(H_A/H_P=15000km)$','$V_i(H_A/i=0^\circ)$','$V_i(H_A/i=63.4^\circ)$','$m/n(H_A/i=0)$','$m/n(H_A/i=63.4^\circ)$','Interpreter','latex')
colororder('gem12')

% Adding specific orbits: Velocity at Apocenter
% Special and general LGO
% clear
clc
Rmean = 6371; % [km] Earth mean radius
mu = 398600.4; % [km^3/s^2] gravitational parameter 
omE = 7.2921159e-5; % [s^-1]
RE = 6378 ; %[km] Earth Equatorial Radius 
J2 = 0.0010827;

H = [0:1:6e4]'; % [km] altitude (at apogee)
rA = H+Rmean; % [km] orbital geocentric radius
% mnindex = and(m_n>0,imag(m_n)==0);



i = [0,63.4]*pi/180; % [rad] inclination
% i = [63.4]*pi/180; % [rad] inclination
m_n_0 = [1,1.5,2]; % Fixed Repetition factor
% phi = 53*pi/180;%[30,50]*pi/180; % [rad] latitude (meridian) of the apogee sub-satellite point
phi = [0,30]*pi/180;
% phi = i;%[30,50]*pi/180; % [rad] latitude (meridian) of the apogee sub-satellite point

% From om
om = [135,90]*pi/180;
Z = 1; % Flag to use om as input [1], to use phi as input [0]
% for l=1:length(i)
%     for o=1:length(om)
%         phi(1,o+length(om)*(l-1)) = asin(sin(om(o)+pi).*sin(i(l))); % [rad]
%     end
% end

clear index
colors = colororder('gem12');
% if phi ~= i
    h = get(gca, 'Children');
    delete(h(1:2));
%     figure
%     grid on
%     hold on
%     % set(gcf,'NumberTitle', 'off', 'Name', 'Fig. 6  On geometric interpretation of locally geostationary orbit design')
%     title("Orbital Velocity vs Altitude")
%     xlabel("\it{H [km]}","fontweight",'bold')
%     ylabel("\it{V [km/s]}","fontweight",'bold')
%     ylim([0,8])
%     xlim([0,6e4])
%     ax = gca;
%     ax.XAxis.Exponent = 0;
% end

J = 0;

for k = 1:length(i)
    % Special LGO
    % velocity at apogee, function of rA, fixed rp
    % verify that Vi(Ha|i) is a good evaluation for LGO identification
    if isequal(phi,i)
        rP(:,k) = (omE^2*rA.^4*cos(i(k))^2)./(2*mu-omE^2*rA.^3*cos(i(k))^2) ; %[km]
        
        a = (rA+rP(:,k))/2; % [km] semimajor axis
        e = (rA-rP(:,k))./(rA+rP(:,k)); % eccentiricy

        K = 1-((3/2)*J2*RE^2)./(a.^2.*(2-e.^2).^2) .* ...
            ((10*cos(i(k))^2-2)./4-(3*cos(i(k))^2-1)./4.*(1-e.^2).^(1/2)+ ...
            (sqrt(mu)*cos(i(k)))./(omE*a.^(7/2).*(1-e.^2).^2));

        m_n(:,k) = (omE*(a.^(3/2)./sqrt(mu).*K)).^(-1); % m/n repetition factor

        Va(:,k) = sqrt((2*mu*rP(:,k))./(rA.*(rA+rP(:,k)))); % [km/s] from eq.3
        mnindex = and(and(m_n(:,k)>0,imag(m_n(:,k))==0),rA>=rP(:,k));
        M_N = m_n(mnindex,k);
        % M_N = m_n(mnindex(:,k),k);
        df = M_N-m_n_0;
        [MIN,index(:,k)] = min(abs(df),[],1);
        INDEX = index(MIN<1e-2,k);
        m_n_0_ = m_n_0(MIN<1e-2);
        [Hindex,~] = find(m_n(:,k)==M_N(INDEX)');
        % Hindex(:,k) = Va(:,k) == Va(:,k); 
        % plot(H(Hindex(:,k)),Va(Hindex(:,k),k),'linewidth',2) % verify it coincides with Vi(Ha|i)
        yyaxis left
        valString = char(strjoin(string(m_n_0_), ', '));
        labelStr_dots = ['$V_A(H_A/m/n=',valString,';~i=', num2str(i(k)*180/pi),'^\circ)$'];
        labelStr_m_n = ['$m/n(H_A/i=',num2str(i(k)*180/pi),'^\circ)$'];
        labelStr_Va = ['$V_A(H_A/i=', num2str(i(k)*180/pi),'^\circ)$'];
        p=plot(H(Hindex),Va(Hindex,k),'o','DisplayName', labelStr_dots);
        colororder("gem12")
        p.MarkerFaceColor = p.Color;
        yyaxis right
        color_id = colors(k,:);
        plot(H(mnindex),M_N,'-.','DisplayName', labelStr_m_n,'Color',color_id);
        plot(H(mnindex),Va(mnindex,k),'--','Linewidth',1.5,'DisplayName', labelStr_Va,'Color',color_id);
        for B = 1:size(Hindex,1)
            J = J+1;
            DATA(J).H = H(Hindex(B));
            DATA(J).Va = Va(Hindex(B),k);
            DATA(J).m_n = m_n_0_(B);
            DATA(J).i = i(k)*180/pi;
            DATA(J).a = a(Hindex(B));
            DATA(J).e = e(Hindex(B));
        end
    else 
        % adding general case
        if Z == 1
            % clear phi
            for o=1:length(om)
                phi(1,o) = asin(sin(om(o)+pi).*sin(i(k))); % [rad]
            end
        end
        for j = 1:length(phi)
            if abs(phi(j)) > abs(i(k)) % invalid i < phi
                continue
            elseif j>1
                if phi(j)==phi(j-1) % avoid repetitions
                    continue
                end
            end
            gamma(k,j) = acos(sqrt(1-sin(i(k))^2+tan(phi(j))^2*cos(i(k))^2))*180/pi; % [deg] angle between the vector Va  and frame velocity vector Vom
            cosg(k,j) = cos(gamma(k,j)*pi/180);
          
            rP(:,k,j) = (omE^2*rA.^4*cos(phi(j))^2)./(2*mu*cosg(k,j).^2-omE^2*rA.^3*cos(phi(j))^2) ; %[km]
            
            a = (rA+rP(:,k,j))/2; % [km] semimajor axis
            e = (rA-rP(:,k,j))./(rA+rP(:,k,j)); % eccentiricy

            K = 1-((3/2)*J2*RE^2)./(a.^2.*(2-e.^2).^2) .* ...
            ((10*cos(i(k))^2-2)./4-(3*cos(i(k))^2-1)./4.*(1-e.^2).^(1/2)+ ...
            (sqrt(mu)*cos(i(k)))./(omE*a.^(7/2).*(1-e.^2).^2));
            
            m_n(:,k,j) = (omE*(a.^(3/2)./sqrt(mu).*K)).^(-1); % m/n repetition factor

            Va(:,k,j) = sqrt((2*mu*rP(:,k,j))./(rA.*(rA+rP(:,k,j)))); % [km/s] from eq.3
            ind = and(and(m_n(:,k,j)>0,imag(m_n(:,k,j))==0),rA>=rP(:,k,j));
            M_N = m_n(ind,k,j);
            df = M_N-m_n_0;
            [MIN,index(:,k,j)] = min(abs(df),[],1);
            INDEX = index(MIN<1e-2,k,j);
            m_n_0_ = m_n_0(MIN<1e-2);
            [Hindex,~] = find(m_n(:,k,j)==M_N(INDEX)');
            color_id = colors(j+length(phi)*(k-1),:);
            % color_id = mod(i-1, size(colors, 1)) + 1;
            yyaxis left
            valString = char(strjoin(string(m_n_0_), ', ')); 
            labelStr_dots = ['$V_A(H_A/m/n=',valString,';~i=', num2str(i(k)*180/pi),'^\circ',';~\varphi=', num2str(phi(j)*180/pi,'%.1f'),'^\circ)$'];
            labelStr_Va = ['$V_A(H_A/i=', num2str(i(k)*180/pi),'^\circ;~\varphi=', num2str(phi(j)*180/pi,'%.1f'),'^\circ)$'];
            labelStr_m_n = ['$m/n(H_A/i=',num2str(i(k)*180/pi),'^\circ;~\varphi=', num2str(phi(j)*180/pi,'%.1f'),'^\circ)$'];
            p=plot(H(Hindex),Va(Hindex,k,j),'o','DisplayName', labelStr_dots,'Color',color_id);
            % colororder("gem12")
            p.MarkerFaceColor = p.Color;
            plot(H(ind),Va(ind,k,j),'--','Linewidth',1.5,'DisplayName', labelStr_Va,'Color',color_id);
            % colororder("gem12")
            yyaxis right
            plot(H(ind),m_n(ind,k,j),'-.','DisplayName', labelStr_m_n,'Color',color_id);
            % colororder("gem12")
            % legend('Interpreter', 'latex')
            for B = 1:size(Hindex,1)
                J = J+1;
                DATA(J).H = H(Hindex(B));
                DATA(J).Va = Va(Hindex(B),k,j);
                DATA(J).m_n = m_n_0_(B);
                DATA(J).i = i(k)*180/pi;
                DATA(J).phi = phi(j)*180/pi;
                DATA(J).a = a(Hindex(B));
                DATA(J).e = e(Hindex(B));
            end
        end
    end
end
% colororder("gem12")
% p.MarkerFaceColor = p.Color;
% Here the equation (29) is solved in rA with a GRID METHOD
% in alternative a Numerical method can be used to solve for rA
% e.g. Newton method 
% fsolve starting rA = 20000 km

%% Plotting extra orbits
clc

Ha = [41274,46895]; %[km]
Hp = [10321,26675]; %[km]

rA = Ha+Rmean;
rP = Hp+Rmean;
for j = 1:length(Ha)
    Va(j) = sqrt((2*mu*rP(j))/(rA(j)*(rA(j)+rP(j))));
    % valString = char(strjoin(string(Ha), ', ')); 
    labelStr_dots = ['$V_A(H_A=',num2str(Ha(j)),'km ;~H_P=',num2str(Hp(j)),'km)$'];
    p=plot(Ha(j),Va(j),'o','DisplayName',labelStr_dots,'Color',colors(6+j,:));
    p.MarkerFaceColor = p.Color;
end
% colororder("gem12")

%% Relate ω argument of perigee with φ

om = 135*pi/180;
i = 63.4*pi/180;
% om = 180/pi*asin((sin(phi)./sin(i)))-pi; % [deg]
phi = 180/pi*asin(sin(om+pi).*sin(i)); % [deg]

%% i as a variable (x) instead of rA (H)

%% phi as variable (x)

%% Table 1 from 170 to 190 deg true anomaly diapason - satellite dwelling in the apogee area
% using Kepler's equation from Curtis Appendix D.2

% given e, a, true anomalies
clear
clc
mu = 398600.4; % [km^3/s^2] gravitational parameter 

e = [0.42, 0.62, 0.26, 0.48]';
a = [42168,32219,42156,32169]'; %[km]

th = [170,190]; %[deg] true anomaly diapason

for k = 1:length(e)
    T(k) = 2*pi*sqrt(a(k)^3/mu); % [s] orbital period
    % E = 2*atan(sqrt((1-e)/(1+e))*tan(th/2*pi/180)); % [rad] eccentric anomaly
    % M = E - e*sin(E); % [rad] Kepler's equation - mean anomaly
    for j = 1:length(th)
        E(j) = 2*atan(sqrt((1-e(k))/(1+e(k)))*tan(th(j)/2*pi/180)); % [rad] eccentric anomaly
        if E(j)<0
            E(j) = 2*pi+E(j);
        end
        M(j) = E(j) - e(k)*sin(E(j)); % [rad] Kepler's equation - mean anomaly
    end
    t(k,:) = M*T(k)/(2*pi); % [s] time since pericenter
    dt(k,:) = diff(t(k,:))/60; % [min] Duration of satellite visibility zone over apogee
end

Tundra = (dt(1)-dt(3))/dt(3)*100; % [%] LGO improvement
Raindrop = (dt(2)-dt(4))/dt(4)*100; % [%] LGO improvement 

%{
% Newton method to solve Kepler's equation (NOT REQUIRED)
error = 1.e-8;
%...Select a starting value for E:
if M < pi
    E = M + e/2; % [rad] eccendtric anomaly
else
    E = M - e/2;
end
%...Iterate on Equation 3.14 until E is determined to within ...the error
%tolerance:
ratio = 1;
while abs(ratio) > error
    ratio = (E - e*sin(E) - M)/(1 - e*cos(E));
    E = E - ratio;
end
%}

 