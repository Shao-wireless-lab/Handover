clear all;close all;clc;

%% Create the matrix of TTT combination
% 16 different TTT values for LTE-to-VLC and VLC-to-LTE
TTT_LV = [0.0 0.04 0.08 0.16 0.32 0.64 1.28 2.56 5.12];
TTT_VL = TTT_LV;
N = length(TTT_LV)*length(TTT_VL);
n = sqrt(N);
% Create the matrix 16x16, each column corresponds to one possible value of
% TTT_LV and each row corresponds to one possible value of TTT_VL
for i = 1:N
        [q,r] = quorem(sym(i),sym(n));
        q=double(q);r=double(r);
        if r~=0
            q = q+1;
        else
            r = n;
        end
        TTT_comb(i,:) = [TTT_LV(q) TTT_VL(r)];
end

%% Create reward matrix for limiting the possible actions at each state (TTT combination)
% We have 9 possible actions:
% 
% Increase TTT_LV by a single level: (i+1)
% Decrease TTT_LV by a single level: (i-1)
% Increase TTT_VL by a single level: (i-n)
% Decrease TTT_VL by a single level: (i+n)
% Increase both TTT_LV and TTT_VL by a single level: (i+n+1)
% Decrease both TTT_LV and TTT_VL by a single level: (i-n-1)
% Increase TTT_LV and decrease TTT_VL by a single level: (i-n+1)
% Increase TTT_VL and decrease TTT_LV by a single level: (i+n-1)
% No change to the current values of TTT_LV and TTT_VL: i

% 24 equal-length time periods
count=1;
for t = [12 24 48 96]
reward=ones(N,N,t);
for i=1:N
    for j=1:N
        for k=1:t
            if j~=i+1  && j~=i-1  && j~=i-n && j~=i+n && j~=i+n+1 && j~=i-n-1 && j~=i-n+1 && j~=i+n-1 && j~=i
              reward(i,j,k)=-Inf;
            end
        end    
    end
end
for i=1:n:N
    for j=1:i+n
        for k=1:t
            if j==i+n-1 || j==i-1 || j==i-n-1
                reward(i,j,k)=-Inf;
                reward(j,i,k)=-Inf;
            end
        end
    end
end
for i=1:N
    for j=1:N
        for k=1:t
            if reward(i,j,k)>0
              reward(i,j,k)=Func_of_Cal_reward_two_AP(TTT_comb(j,1),TTT_comb(j,2),k,t);
            end
        end    
    end
end
filename = sprintf('t%dsmall.mat',t);
save(filename);
disp(t)
%% Q-learning algorithm

% Initialize the Q-table with random values ~ N(0,1)
% Set learnning rate to 1 and discount factor to 0.9
% The maximum number of episodes is set to 50

q = randn(size(reward));
gamma = 0.9;
alpha = 1;
maxItr = 10000;
epsilon_initial = 0.2;
% cs -> current state
% ns -> next state
% 
% Repeat until Convergence OR Maximum Iterations
for i=1:maxItr
    epsilon = epsilon_initial/(1+i/500);
    
    % Starting from start position    
    cs=73;
    
    % Repeat for t times
    for k=1:t
        
    % possible actions for the chosen state
    n_actions = find(reward(cs,:,k)>=0);
    % choose an action at random with probability epsilon and set it as the
    % next state
    if rand(1)<epsilon
        ns = n_actions(randi(length(n_actions)));
    else 
        ns = n_actions(find(q(cs,n_actions,k)==max(q(cs,n_actions,k))));
        if length(ns)>1
            ns = ns(randi(length(ns)));
        end
    end
        if k<t
            % find all the possible actions for the selected state
            n_actions = find(reward(ns,:,k+1)>=0);
            
            % find the maximum q-value i.e, next state with best action
            max_q = 0;
            for j=1:length(n_actions)
                max_q = max(max_q,q(ns,n_actions(j),k+1));
            end
            
            % Update q-values as per Bellman's equation
            q(cs,ns,k)=reward(cs,ns,k)+gamma*max_q;
            %fprintf('q(%d,%d)=%d\n',cs,ns,q(cs,ns));
        else
            q(cs,ns,k)=reward(cs,ns,k);
        end 
            % Set current state as next state
            cs=ns;
    Throughput(k) = Func_of_Cal_reward_two_AP(TTT_comb(cs,1),TTT_comb(cs,2),k,t);
    end
    
Throughput_Average(i,count) = mean(Throughput);
       
end
count=count+1;
end

figure
hold on
plot(1:maxItr,Throughput_Average(:,1),'r');
plot(1:maxItr,Throughput_Average(:,2),'b');
plot(1:maxItr,Throughput_Average(:,3),'g');
plot(1:maxItr,Throughput_Average(:,4),'k');
hold off
legend('t = 12','t = 24','t = 48','t = 96');
xlabel('Episode index');
ylabel('Average throughput (Mbps)');
grid on
box on


% cs =241;
% for k=1:t
%     ns = find(q(cs,:,k)==max(q(cs,:,k)));
%     if length(ns)>1
%         ns = ns(randi(length(ns)));
%     end
%     Can_TTT_LV(k) = TTT_comb(ns,1);
%     Can_TTT_VL(k) = TTT_comb(ns,2);
%     cs = ns;
% end
% figure
% hold on
% plot(1:t,Can_TTT_LV,'r');
% plot(1:t,Can_TTT_VL,'b');
% hold off
