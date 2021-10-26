%% AAE6102 Assignment
% Author:       JIANG Bailun
% Date:         18/10/2021
%% Calculate user position and clock bias based on ephemeris and pseudorange

clear all

% Input parameters
approx = [-2694685.473;...  % Initial user position and clock bias
          -4293642.366;...
          3857878.924;...
          0];
dx_threshold = 1e-4;        % Threshold to stop the iteration

% Constants
global c
global wedot
global mu
global F

c = 299792458;              % Speed of light (m/s)
wedot = 7.2921151467e-5;    % WGS 84 value of earth's rotation rate (r/s)
mu = 3.986005e+14;          % WGS 84 value of earth's universal gravitation constant (m^3/s^2)
F = -4.442807633e-10;       % Relativistic correction term constant

% Load data from files
load rcvr.dat
load eph.dat

% Sort data since the orders in two files are different
eph = sortrows(eph,2);
rcvr = sortrows(rcvr,2);

% Calculate satellite position and clock bias
sv_num = size(eph,1);       % Satellite number
for i = 1:sv_num
    sv_data(i,:) = get_sv_pos(eph(i,:),rcvr(i,3));
end

% Correct the pseudorange
sv_data(:,4) = rcvr(:,3) + c .* sv_data(:,4);
% sv_data(:,4) = rcvr(:,3);

% Use least square to calculate user position and clock bias
dx = ones(1,4);             % Initialize dx to start while loop
i = 0;                      % Number of iterations
disp('---------------------------------------------------------------------------------------')
disp('                  Iteration starts...')
disp('Initial user position in WGS 84 XYZ coordinate (meters) and user clock bias (seconds): ')
disp(num2str(approx(1:3)))
disp(num2str(approx(4)/c))

while norm(dx(1:3)) > dx_threshold
    dx = estimate_user_pos(sv_data,approx);
    approx = approx + dx;   % Update approximation
    i = i+1;
    D = ['Iteration ',num2str(i),' result:'];
    disp(D)
    disp(num2str(approx(1:3)))
    disp(num2str(approx(4)/c))
end

% Print the results
disp('Iteration ends with reached threshold')
disp('---------------------------------------------------------------------------------------')


%% Calculate satellite position and clock bias based on emphemeris and time
function [x] = get_sv_pos(eph,pr)

    % Input parameters:
    rcvr_tow    = eph(1);	% Receiver time of week (s)
    svid        = eph(2);	% Satellite PRN number (1-32)
    toc         = eph(3);	% Reference time of clock parameters (s)
    toe         = eph(4);	% Reference time of ephemeris parameters (s)
    af0         = eph(5);   % Clock correction coefficient - group delay (s)
    af1         = eph(6);	% Clock correction coefficient (s/s)
    af2         = eph(7);	% Clock correction coefficient (s/s/s)
    ura         = eph(8);	% User range accuracy (m)
    e           = eph(9);	% Eccentricity (-)
    sqrta       = eph(10);	% Square root of semi-major axis a (m^1/2)
    dn          = eph(11);	% Mean motion correction (r/s)
    m0          = eph(12);	% Mean anomaly at reference time (r)
    w           = eph(13);	% Argument of perigee (r)
    omg0        = eph(14);	% Right ascension (r)
    i0          = eph(15);  % Inclination angle at reference time (r)
    odot        = eph(16);  % Rate of right ascension (r/s)
    idot        = eph(17);  % Rate of inclination angle (r/s)
    cus         = eph(18);  % Argument of latitude correction, sine (r)
    cuc         = eph(19);  % Argument of latitude correction, cosine (r)
    cis         = eph(20);  % Inclination correction, sine (r)
    cic         = eph(21);  % Inclination correction, cosine (r)
    crs         = eph(22);  % Radius correction, sine (r)
    crc         = eph(23);  % Radius correction, cosine (r)
    iod         = eph(24);  % Issue of data number

    global c
    global wedot
    global mu
    global F
    
    % Calculation follows Table 20-IV in ICD file
    a = sqrta^2;            % Semi-major axis
    n0 = sqrt(mu/a^3);      % Computed mean motion (r/s)
    t = rcvr_tow - pr/c;    % Satellite signal transmition time
    tk = t - toe;           % Time from ephemeris reference epoch
    
    % Account for beginning or end of week crossovers
    if tk > 302400
        tk = tk - 604800;
    elseif tk < -30240
        tk = tk + 604800;
    end
    
    n = n0+dn;              % Corrected mean motion
    mk = m0+n*tk;           % Mean anomaly
    
    % Solve eccentric anomaly
    syms ex
    eqn = ex - e*sin(ex) == mk;
    ek = vpasolve(eqn);
    ek = double(ek);
    clear ex eqn
    
    % True anomaly
    vk = atan2((sqrt(1-e^2)*sin(ek)/(1-e*cos(ek))),((cos(ek)-e)/(1-e*cos(ek))));
    
    % Eccentric anomaly
    phik = vk + w;          % Argument of latitude
    
    % Second harmonic perturbations
    duk = cus*sin(2*phik) + cuc*cos(2*phik);    % Argument of Latitude Correction
    drk = crs*sin(2*phik) + crc*cos(2*phik);    % Radius Correction
    dik = cis*sin(2*phik) + cic*cos(2*phik);    % Inclination Correction
    
    uk = phik + duk;                            % Corrected argument of latitude
    rk = a*(1 - e*cos(ek)) + drk;                % Corrected Radius
    ik = i0 + dik + idot*tk;                    % Corrected Inclination
    xkp = rk*cos(uk);                           % X position in orbital plane
    ykp = rk*sin(uk);                           % Y position in orbital plane
    omgk = omg0 + (odot - wedot)*tk -wedot*toe; % Corrected longitude of ascending node
    
    xk = xkp*cos(omgk) - ykp*cos(ik)*sin(omgk); % Earth-fixed X coordinates
    yk = xkp*sin(omgk) + ykp*cos(ik)*cos(omgk); % Earth-fixed Y coordinates
    zk = ykp*sin(ik);                           % Earth-fixed Z coordinates
    
    % Calculate satellite clock bias (s)
    dtsv = af0 + af1*(t - toc) + af2*(t - toc)^2 + F*e*sqrta*sin(ek);
    
    x = [xk yk zk dtsv];
end

%% Estimate user position and clock bias based on pseudorange and satellite position
function [dx] = estimate_user_pos(sv_data,approx)
    
    % Input parameters
    sv_x = sv_data(:,1);    % Satellite X position
    sv_y = sv_data(:,2);    % Satellite Y position
    sv_z = sv_data(:,3);    % Satellite Z position
    pr = sv_data(:,4);      % Measured pseudorange
    x0 = approx(1);         % Initial X user position
    y0 = approx(2);         % Initial Y user position
    z0 = approx(3);         % Initial Z user position
    b0 = approx(4);         % Initial user clock bias (m)
    
    % Calculate delta rho
    rho = pr;               % Current (measuren) pseudorange
    r = sqrt((sv_x - x0).^2 + (sv_y - y0).^2 + (sv_z - z0).^2);
                            % Geometric range between user and satellite
    rho_hat = r - b0;       % Last (approximated) pseudorange
    delta_rho = rho_hat - rho;
    
    % Calculate H matrix
    ax = (sv_x - x0)./r;
    ay = (sv_y - y0)./r;
    az = (sv_z - z0)./r;
    H = [ax ay az ones(size(sv_data,1),1)];
    
    % Calculate dx
    dx = inv(H'*H)*H'*delta_rho;

end
