
 
Dorsal = logspace(1,3.5,10); % in AUs
Kd = 4000;
timeVector = [1:10:12*60]; % in minutes
k1 = 0.01;
k2 = 0.001;
 
figure
hold on
for D = Dorsal
    pOn1 = OneRevStep(k1*DlBinding(D,Kd),k2,timeVector);
    plot(timeVector,pOn1)
end
hold off
 
 
 
function [On_t,Off_t] = OneRevStep(k1,k2,timeVector)
    On_t = (k1 - exp((-k1-k2).*timeVector) .* k1) ./ (k1 + k2);
    Off_t = 1 - On_t;
end
 
function Occupancy = DlBinding(Dl,Kd)
    Occupancy = (Dl./Kd) ./ (1 + (Dl./Kd));
end

