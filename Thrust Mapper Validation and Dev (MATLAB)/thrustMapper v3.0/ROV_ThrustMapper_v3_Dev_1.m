%% ROV Triton: Thrust Mapper Understanding and Development
%  thrustMapper v3.0
%  Development Driver 1
%  2021.05.23

%  Playing around with the nullspace basis vectors

close all; clear; clc; tic;
PLOT_NULLBASIS = false;

tm = thrustMapper([0,0,0]);

%% TEST CASES
testCase = 1;
switch(testCase)
    case 1
        F_des = 120*[1;0;0]; %[N]
        thrustList = tm.limitedMap(F_des,[0;0;0]);
    case 2
        F_des = 100*[0;0;1]; %[N]
        thrustList = tm.limitedMap(F_des,[0;0;0]);
    case 3
        F_des = 50*[1;1;0]; %[N]
        thrustList = tm.limitedMap(F_des,[0;0;0]);
    case 4
        M_des = 10*[1;0;0]; %[N
        thrustList = tm.limitedMap([0;0;0],M_des);
end
fprintf('%.4f | Mapping complete.\n',toc);

%% ADD NULL SPACE COMPONENT
b1 = tm.nullBasis(:,1);
b2 = tm.nullBasis(:,2);
a1 = 30; %Coeff. on b1
a2 = 0;  %Coeff. on b2
thrustList_NullMod = thrustList + a1*b1 + a2*b2;

%% PLOT
tm.plotSetup();

% Plot just the basis vectors of the nullspace
if(PLOT_NULLBASIS)
    figure;
    subplot(1,2,1);
    tm.plotThrusters(50*b1,'Nullspace Basis Vector b_1 (Scaled by 50)');
    view([-120,25]);
    subplot(1,2,2);
    tm.plotThrusters(50*b2,'Nullspace Basis Vector b_2 (Scaled by 50)');
    view([-120,25]);
end

% Plot the thrustList w/ and w/o the added null component
figure;
subplot(1,2,1);
tm.plotThrusters(thrustList,'Result of Pseudo-Inverse, T');
gcaExpandable();
subplot(1,2,2);
tm.plotThrusters(thrustList_NullMod,sprintf('Adding null component, T + %db_1 + %db_2',a1,a2));
gcaExpandable();

%figure;
%tm.plotThrusters([0;0;thrustList_NullMod(3:4); 0;0;0;0]);
fprintf('%.4f | Program complete.\n',toc);