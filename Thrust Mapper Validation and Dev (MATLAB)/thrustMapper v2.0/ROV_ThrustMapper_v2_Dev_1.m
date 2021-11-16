%% ROV Triton: Thrust Mapper Understanding and Development
%  thrustMapper v2.0
%  Development Driver 1
%  2021.05.09

close all; clear; clc;

tm = thrustMapper([0,0,0]);
%tm.plotSetup();

%% INDIVIDUAL TEST CASES
testCase = 1;
switch(testCase)
    case 1
        F_des = 120*[1;0;0]; %[N]
        thrustList1 = tm.dumbMap(F_des,[0;0;0]);
        thrustList2 = tm.limitedMap(F_des,[0;0;0]);
    case 2
        F_des = 100*[0;0;1]; %[N]
        thrustList1 = tm.dumbMap(F_des,[0;0;0]);
        thrustList2 = tm.limitedMap(F_des,[0;0;0]);
    case 3
        F_des = 50*[1;1;0]; %[N]
        thrustList1 = tm.dumbMap(F_des,[0;0;0]);
        thrustList2 = tm.limitedMap(F_des,[0;0;0]);
    case 4
        M_des = 10*[1;0;0]; %[N
        thrustList1 = tm.dumbMap([0;0;0],M_des);
        thrustList2 = tm.limitedMap([0;0;0],M_des);
end

figure;
subplot(1,2,1);
tm.plotThrusters(thrustList1,'No Limiting');
gcaExpandable();
subplot(1,2,2);
tm.plotThrusters(thrustList2,'With Limiting');
gcaExpandable();