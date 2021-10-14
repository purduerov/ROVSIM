%% ROV Triton: Thrust Mapper Understanding and Development
%  thrustMapper v3.0
%  Development Driver 2a
%  2021.05.23-26

% Testing thrustMapper.nullMap()

close all; clear; clc; tic;

tm = thrustMapper([0,0,0]);

%% TEST CASES
testCase = 1;
switch(testCase)
    case 1
        F_des = 120*[1;0;0]; %[N]
        M_des = [0;0;0]; %[Nm]
    case 2
        F_des = 100*[0;0;1]; %[N]
        M_des = [0;0;0]; %[Nm]
    case 3
        F_des = 50*[1;1;0]; %[N]
        M_des = [0;0;0]; %[Nm]
    case 4
        F_des = [0;0;0]; %[N]
        M_des = 10*[1;0;0]; %[N]
    case 5
        F_des = 45*[0;1;0]; %[N]
        M_des = [0;0;0]; %[Nm]
end
thrustList1 = tm.limitedMap(F_des,[0;0;0]);
[thrustList2, numExceedsLimit] = tm.nullMap(F_des,[0;0;0]);
fprintf('%.4f | Mapping complete.\n',toc);

%% PLOT
tm.plotSetup();

% Plot the thrustList w/ and w/o the added null component
figure;
subplot(1,2,1);
tm.plotThrusters(thrustList1,'limitedMap()',F_des,M_des);
gcaExpandable();
subplot(1,2,2);
tm.plotThrusters(thrustList2,'nullMap()',F_des,M_des);
gcaExpandable();

fprintf('%.4f | Program complete.\n',toc);