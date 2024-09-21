% ---------------------------------------------
% ----- INFORMATIONS -----
%   Author          : louis tomczyk
%   Institution     : Telecom Paris
%   Email           : louis.tomczyk@telecom-paris.fr
%   Version         : 1.1.0
%   Date            : 2024-09-19
%   License         : cc-by-nc-sa
%                       CAN:    modify - distribute
%                       CANNOT: commercial use
%                       MUST:   share alike - include license
%
% ----- CHANGELOG -----
%   2023-09-23  (1.0.0) creation
%   2024-09-19  (1.0.1) better display (axis lim)
% 
% ----- INPUTS -----
% ----- EVOLUTIONS -----
% - define by spectra instead by time pulses
% - display spectra evolution
%
% ----- BIBLIOGRAPHY -----
%   Functions           :
%   Author              :
%   Author contact      :
%   Date                :
%   Title of program    :
%   Code version        :
%   Type                : 
%   Web Address         : 
% ----------------------------------------------
%%


% ----------------------------------------------
% MAINTENANCE
% ----------------------------------------------
clear
close all
clc

% ----------------------------------------------
% PARAMETERS TO CHANGE BY USER
% ----------------------------------------------
% what spectrum you want?
% choose between:
% - "rect"  - rectangular spectrum
% - "gauss" - gaussian spectrum
what    =  "gauss";

% if you chose "gauss", how much energy percentage
% you wanna use for pulse duration estimation?
energy  = 0.999;

% spectral characteristics
lambda  = 1550e-9;      % wavelength        [m]
dnu     = 1e8;          % linewdith         [Hz]

% set the distance over which you want the pulses
% to propagate, in terms of coherence lengths
N_lc    = 1;

% ----------------------------------------------
% PHYSISCAL PARAMETERS
% ----------------------------------------------
c   = 299752458;        % light celerity    [m/s]
nu  = c/lambda;         % carrier frequency [Hz]

% ----------------------------------------------
% TIME AXIS
% ----------------------------------------------
% there will be 3 plots:
% - one with the two sources independent
% - one with the power of the sum of the sources
% - one which is a zoom inside the enveloppe in
%       order to see the phase
% for this last plot we need to zoom in the pulse
% so we draw from the principal time axis a second
% time axes
% (1) creation of the principal axis that spreads
%     over 2*N_PERIODS_E-1 pseudo period of the
%     sinc(.).
% (2) setting the INDEXES in the principal time
%     array such as we zoom inside the mean lobe
N           = 4e5;      % number of points
n_periods_e = 5;        % number of periods to display

if strcmpi(what,'rect') == 1
    Te      = 2/dnu;            % Enveloppe's period  [s]
    tau_c   = Te;
    dt      = n_periods_e*Te/N; % [s]
    n_points_per_period_e = Te/dt;
    lobe_number = 2;
    indexes     = lobe_number*n_points_per_period_e/2+1:...
                 (lobe_number+1)*n_points_per_period_e/2;
elseif strcmpi(what,'gauss') == 1

    % small exercise: f(t=tau/2)=alpha*f(t=0)
    Te          = sqrt(-2*log(1-energy))/pi/dnu;
    n           = 1024;
    indexes     = N/2-int32(N/n):N/2+int32(N/n);
end

t = linspace(-n_periods_e*Te/2,n_periods_e*Te/2,N);

if strcmpi(what,'mono') ~= 1
    tbis        = t(indexes);
end

l_c     = c*Te;             % [m]   coherence length
dist_max= N_lc*l_c;         % [m]   travelling distance of the pulses
FPS     = 30;               % [1/s] Frames Per Second
Nsteps  = 60;               %       dist_max = Nsteps * dz
dz      = dist_max/Nsteps;  % [m]
dt      = dz/c;             % [s]


figure('units','normalized','outerposition',[0 0 1 1])
for step = 1:Nsteps

    dist= step*dz; % [m] propagation distance
    tg  = step*dt; % [s] equivalent time delay

    % time signals
    if strcmpi(what,'mono') == 1
        x1      = cos(2*pi*nu*t);
        x2      = cos(2*pi*nu*(t-tg));
    elseif strcmpi(what,'rect') == 1
        x1      = sinc(dnu*t).*cos(2*pi*nu*t);
        x2      = sinc(dnu*(t-tg)).*cos(2*pi*nu*(t-tg));
    elseif strcmpi(what,'gauss') == 1
        x1      = exp(-pi*dnu^2*t.^2).*cos(2*pi*nu*t);
        x2      = exp(-pi*dnu^2*(t-tg).^2).*cos(2*pi*nu*(t-tg));
    end

    x   = x1+x2;
    y   = abs(x).^2;

    clf
    subplot(2,2,1)
        hold on
        plot(t*1e9,x1)
        plot(t*1e9,x2)
        xlabel("time [ns]")
        ylabel('x1, x2')
        xlim([t(1)*1e9,t(end)*1e9])
        title(sprintf("propagation distance = %.2e [m]",dist))

    if strcmpi(what,'rect') == 1 || strcmpi(what,'gauss') == 1
        subplot(2,2,[3,4])
            hold on
            plot(tbis*1e9,x1(indexes))
            plot(tbis*1e9,x2(indexes))

            imin    = round(length(tbis)/3);
            idt      = round(length(tbis)/100);
            xlim([tbis(imin)*1e9,tbis(imin+idt)*1e9])
            if strcmpi(what,'gauss')
                ylim([-1.1,1.1])
            else
                ylim([-1,1]/5)
            end
            xlabel("time [ns]")
            ylabel('x1, x2')
    end

    subplot(2,2,2)
        plot(t*1e9,y)
        axis([t(1)*1e9,t(end)*1e9,-0.1,4.1])
        xlabel("time [ns]")
        ylabel('|x1+x2|^2')
        title(sprintf("equivalent time delay = %.2e [mu s]",tg*1e6))
    pause(1/FPS)

end