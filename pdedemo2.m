%% Helmholtz's Equation on a Unit Disk with a Square Hole
% This example shows how to solve a Helmholtz equation using the |solvepde|
% function in Partial Differential Toolbox(TM).
%
% The Helmholtz equation, an elliptic equation that is the time-independent
% form of the wave equation, is
%
% $-\Delta u-k^2u = 0$.
%
% Solving this equation allows us to compute the waves reflected by a
% square object illuminated by incident waves that are coming from the
% left.

%       Copyright 1994-2015 The MathWorks, Inc.

%% Problem Definition
% The following variables define our problem:
%
% * |g|: A geometry specification function. For more information, see the
% code for |scatterg.m| and the documentation section
% <matlab:helpview(fullfile(docroot,'toolbox','pde','helptargets.map'),'pde_geometry_fcn');
% Create Geometry Using a Geometry Function>.
% * |k|, |c|, |a|, |f|: The coefficients and inhomogeneous term.
g = @hasab;
k = 60;
c = 1;
a = -k^2;
f = 0;

disp('data loaded');

%% Create PDE Model
% Create a PDE Model with a single dependent variable.
numberOfPDE = 1;
model = createpde(numberOfPDE);

disp('model created');

%%
% Convert the geometry and append it to the pde model.
geometryFromEdges(model,g);

disp('geometry converted');

%% Specify PDE Coefficients
specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',f);

disp('pde specified');

%% Boundary Conditions
% Plot the geometry and display the edge labels for use in the boundary
% condition definition.
figure; 
pdegplot(model,'EdgeLabels','on'); 
axis equal
title 'Geometry With Edge Labels Displayed';
ylim([-5,5])

disp('boundary defined');

%%
% Apply the boundary conditions.
bOuter = applyBoundaryCondition(model,'neumann','Edge',(1:4),'g',0,'q',-1i*k);
innerBCFunc = @(loc,state)-exp(-1i*k*loc.x);
bInner = applyBoundaryCondition(model,'dirichlet','Edge',(5:8),'u',innerBCFunc);

disp('boundary applied');

%% Create Mesh
generateMesh(model,'Hmax',2*pi/k/10);
figure
pdemesh(model); 
axis equal

disp('mesh created');

%% Solve for Complex Amplitude
% The real part of the vector |u| stores an approximation to a real-valued solution of the
% Helmholtz equation.
result = solvepde(model);
u = result.NodalSolution;

disp('complex amplitude');

%% Plot FEM Solution
figure
% pdeplot(model,'XYData',real(u),'ZData',real(u),'Mesh','off');
pdeplot(model,'XYData',real(u),'Mesh','off');
colormap(jet)

disp('fem solution plotted');

%% Animate Solution to Wave Equation
% Using the solution to the Helmholtz equation, construct an animation showing
% the corresponding solution to the time-dependent wave equation.
figure
m = 10;
h = newplot; 
hf = h.Parent; 
axis tight
ax = gca;
ax.DataAspectRatio = [1 1 1];
axis off
maxu = max(abs(u));
for j = 1:m
    uu = real(exp(-j*2*pi/m*sqrt(-1))*u);
    pdeplot(model,'XYData',uu,'ColorBar','off','Mesh','off');
    colormap(jet)
    caxis([-maxu maxu]);
    axis tight
    ax = gca;
    ax.DataAspectRatio = [1 1 1]; 
    axis off
    M(j) = getframe(hf);
end
%%
% To play the movie 10 times, use the |movie(hf,M,10)| command.
