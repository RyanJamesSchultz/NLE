% Plot the progession of Mlrg for both a constant injection rate and modified rate.
% Used to make Figure 3.
clear;

% Define some constants.
N=1e4;
Ns=1e2;
m1=0.0;
m2=Inf;
b=1.0;
c=2;

% Derive some more variables.
T1=0:N;
v1=c*ones(size(T1));
V1=cumsum(v1);
Mmod=m1+log10(T1)/b;

% Tranform to the new time-axis.
v2=exp(T1*(6/N));
V2=cumsum(v2);
v2=v2*(V1(end)/V2(end)); v2(v2>c)=c; V2=cumsum(v2); v2=v2*(V1(end)/V2(end));
V2=cumsum(v2);
T2=interp1(cumsum(v2),T1,cumsum(v1),'linear');

% Preallocate.
Mlrg=zeros([N+1 Ns]);

% Loop to get multiple Mlrg sequences.
for i=1:Ns
    % Get the largest event sequence.
    M=GR_MFD_Rand(m1,m2, log10(N),b, [1 N]);
    Mlrg(2:end,i)=OrderStatistic(M,N,'none');
    Mlrg(1,i)=m1;
end

% Plot.
figure(3); clf;
plot(T1,Mlrg,'-r'); hold on;
plot(T2,Mlrg,'-m');
plot(T1,mean(Mlrg,2),'-k');
plot(T2,mean(Mlrg,2),'-k');
plot(T1,Mmod,':k');
plot(T1,v1,'-b');
plot(T1,v2,'-c');
xlabel('Time'); ylabel('Magnitude');
ylim([m1 m2]); xlim([min(T1) max(T1)])

% Extra injection-related plots.
figure(4); clf;
subplot(131);
plot(T1,v1); hold on;
plot(T2,v2);
subplot(132);
plot(T1,cumsum(v1)); hold on;
plot(T2,cumsum(v2));
subplot(133);
plot(T1,T2);

