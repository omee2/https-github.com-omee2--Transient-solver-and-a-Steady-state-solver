
close all
clc
clear all

%plot foe the n number of second_queation4
figure(1)  % For n = 20 
second_question4(20, [0 1 0.8 2.4])
figure(2)  % For n = 40 
second_question4(40, [0 1 0.8 2.4])
figure(3)  % For n = 80
second_question4(80,[0 1 0.8 2.4])
figure(4)  % For n = 160
second_question4(160,'tight')
figure(5)  % For n = 101 i.e 100 points, and dx = 0.01 - Optimum stability
second_question4(100,[0 1 0.8 2.4])