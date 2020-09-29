function [ epsilon, tmax, hmax, treturn ] = projectile( V )
%projectile  For "projectile problem," calculate time to max height and
%            max height and plot perturbation approx vs true solution.
% 
%  >> [ epsilon, tmax, hmax, treturn ] = projectile( V )
% 
%  in:  V = initial velocity of projectile (in miles per hour)
% 
% out: epsilon = key dimensionless param ( = V^2 / g*R )
%         tmax = time to max height (in reduced units)
%         hmax = max height (" " ")
%      treturn = time to return to earth (should be 2*tmax)
% 
% also output is a plot (in the current figure window) of the true
% solution vs the approximate solution corresponding to a constant
% gravitational field (epsilon=0)

% ~/matlab/courses/AppliedMath/projectile.m
% 
% C.G. (10-03)
% C.G. (10-07) minor change: removed "IE" from return arg for ode45()

% parameters and intializations
V=100000;
g = 78500 ;                        % grav const (English units)
R =  4000 ;                        % earth radius ( " " )

epsilon = V^2 / (g*R) ;            % key dimensionless param

tspan = [ 0 3 ] ;                  % interval of time integration

y0 = [ 0 ; 1 ] ;                   % ICs: y0(1)=h(0)=0, y0(2)=h'(0)=1

% set some options for the ode solver: smaller error tol (anticipate
% errors around 1e-8), output # integration steps (etc, for information),
% handle/pointer to "events" function

options = odeset( ...
    'AbsTol', 1e-18,  ...
    'RelTol', 1e-11,  ...
    'Stats',  'on',   ...
    'Events', @events ) ;

% call the ode45() solver

[ T, Y, TE, YE ] = ode45( @odefun, tspan, y0, options, epsilon ) ;

% calculate return values

tmax = TE(1) ;                     % time to max ht (1st event)

hmax = YE(1,1) ;                   % value of y(1)=h at max ht

treturn = TE(2) ;                  % time to return to earth (2nd event)

% plot the computed solution vs the exact solution of the simplified
% model that assumes a constant graviational field

Ysimp =  T-.5*T.^2 ;               % solution of simplified model

plot( T, Y(:,1), T, Ysimp, '--' )  % plot both in current figure window

% labels, titles, dressing up, ...

legend( 'True Model', 'Simplified Model' )

set( gca, 'FontSize', 12 )

xlabel( 't' )
ylabel( 'h' )

title_str = [ 'V = ', num2str(V), ',  \epsilon = ', num2str(epsilon), ...
        ',  tmax = ', num2str(tmax),  ',  hmax = ', num2str(hmax) ] ;

title( title_str, 'FontSize', 14 )

axis( [ T(1) T(end) 0 hmax ] )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dy = odefun( t, y, epsilon )

dy = [ y(2) ;  -1 / ( 1+epsilon*y(1) )^2 ] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ value, isterminal, direction ] = events( t, y, epsilon )

% event(1) = max ht ( y(2)=0 )
% event(2) = return to earth ( y(1)=0 with decreasing y(1) )

value = [ y(2) ; y(1) ] ;

isterminal = [ 0 ; 1 ] ;           % terminate upon return to earth

direction = [ 0 ; -1 ] ;           % cat h=0 only when decreasing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% representative output:
% 
%  V = 100 :    epsilon = 3.184713375796178e-05
%               tmax    = 1.00002123182822
%               hmax    = 0.50000796191022
%               treturn = 2.00004246365641
% 
%  V = 1000 :   epsilon = 0.00318471337580
%               tmax    = 1.00212720661000
%               hmax    = 0.50079744813928
%               treturn = 2.00425441309215
% 
%  V = 10000 :  epsilon = 0.31847133757962
%               tmax    = 1.26184686114108
%               hmax    = 0.59469696969648
%               treturn = 2.52369372228136