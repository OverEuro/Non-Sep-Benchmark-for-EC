clc,clear
% function [DE_gbest,DE_gbestval,DE_fitcount,DE_fit_cut,DE_get_flag]  = newjDE(fname,...
%     VTR,Max_FES,D,XRmin,XRmax,Lbound,Ubound,NP,Max_Gen,F,CR,strategy,fun,M); 
ftest = @f11;
Max_FES = 3 * 10^5;
D = 30;
Lbound = -30;
Ubound = 30;
XRmin = Lbound * ones(1,D);
XRmax = Ubound * ones(1,D);
NP = 30;
DE_fit_cut = Max_FES;
t1=0.1;
t2=0.1;

F=0.5*ones(NP,1);
CR=0.9*ones(NP,1);
%-----Initialize population and some arrays-------------------------------
pop=repmat(XRmin,NP,1)+(repmat(XRmax,NP,1)-repmat(XRmin,NP,1)).*rand(NP,D);

popold    = zeros(size(pop));     % toggle population
val       = zeros(1,NP);          % create and reset the "cost array"
nfeval    = 0;                    % number of function evaluations

%------Evaluate the best member after initialization----------------------

ibest   = 1;                      % start with first population member
val(1)  = ftest(pop(ibest,1:D)); 
DE_gbestval = val(1);                 % best objective function value so far
nfeval  = nfeval + 1;
for i=2:NP                        % check the remaining members
    val(i) = ftest(pop(i,1:D)); 
    if (val(i) < DE_gbestval)           % if member is better
        ibest   = i;                 % save its location
        DE_gbestval = val(i);
    end   
end
DE_gbest = pop(ibest,1:D);         % best member of current iteration
bestvalit = DE_gbestval;              % best value of current iteration

%------DE-Minimization---------------------------------------------
%------popold is the population which has to compete. It is--------
%------static through one iteration. pop is the newly--------------
%------emerging population.----------------------------------------

pm1 = zeros(NP,D);              % initialize population matrix 1
pm2 = zeros(NP,D);              % initialize population matrix 2
pm3 = zeros(NP,D);              % initialize population matrix 3
bm  = zeros(NP,D);              % initialize bestmember  matrix
ui  = zeros(NP,D);              % intermediate population of perturbed vectors
mui = zeros(NP,D);              % mask for intermediate population
mpo = zeros(NP,D);              % mask for old population
rot = (0:1:NP-1);               % rotating index array (size NP)
rt  = zeros(NP);                % another rotating index array
a1  = zeros(NP);                % index array
a2  = zeros(NP);                % index array
a3  = zeros(NP);                % index array

ind = zeros(2);

iter = 0;
pop(:,D+1)=F;
pop(:,D+2)=CR;
while nfeval < Max_FES
    popold = pop(:,1:D);                   % save the old population
    
    ind = randperm(2);              % index pointer array
    
    a1  = randperm(NP);             % shuffle locations of vectors
    rt = rem(rot+ind(1),NP);        % rotate indices by ind(1) positions
    a2  = a1(rt+1);                 % rotate vector locations
    rt = rem(rot+ind(2),NP);
    a3  = a2(rt+1);                               
    
    pm1 = popold(a1,:);             % shuffled population 1
    pm2 = popold(a2,:);             % shuffled population 2
    pm3 = popold(a3,:);             % shuffled population 3
    
    for i=1:NP                      % population filled with the best member
        bm(i,1:D) = DE_gbest(1:D);          % of the last iteration
    end
    
    for i=1:NP
        if rand<t1,
            %             F(i)=0.1+rand*0.9;
            ui(i,D+1)=0.1+rand*0.9;    
        else
            ui(i,D+1)=pop(i,D+1);
        end
        if rand<t2,
            %             CR(i)=rand;
            ui(i,D+2)=rand;
        else
            ui(i,D+2)=pop(i,D+2);
        end
    end
    F=ui(:,D+1);
    CR=ui(:,D+2);
    mui = rand(NP,D) < repmat(CR,1,D);          % all random numbers < CR are 1, 0 otherwise
    
    dd=ceil(D*rand(NP,1));
    for kk=1:NP
        mui(kk,dd(kk))=1;
    end
    mpo = mui < 0.5;                % inverse mask to mui

    ui(:,1:D) = pm3 + repmat(F,1,D).*(pm1 - pm2);       % differential variation
    ui(:,1:D) = popold.*mpo + ui(:,1:D).*mui;     % crossover
        
    %-----Select which vectors are allowed to enter the new population------------
    for i=1:NP
        outbind=find(ui(i,1:D) < Lbound);
        if size(outbind,2)~=0
            ui(i,outbind)=XRmin(outbind);
        end            
        outbind=find(ui(i,1:D) > Ubound);
        if size(outbind,2)~=0
            ui(i,outbind)=XRmax(outbind);
        end
    end
    
    for i=1:NP
        tempval = ftest(ui(i,1:D));   % check cost of competitor
        nfeval  = nfeval + 1;
        if (tempval <= val(i))  % if competitor is better than value in "cost array"
            pop(i,:) = ui(i,:);  % replace old vector with new one (for new iteration)
            val(i)   = tempval;  % save value in "cost array"

            
            %----we update DE_gbestval only in case of success to save time-----------
            if (tempval < DE_gbestval)     % if competitor better than the best one ever
                DE_gbestval = tempval;      % new best value
                DE_gbest = ui(i,:);      % new best parameter vector ever
            end
        end
    end %---end for imember=1:NP
    
    iter = iter + 1;
    Res(iter) = DE_gbestval;
    fprintf('FEs=%d Fmin=%g\n',nfeval,DE_gbestval);

end %---end while ((iter < Max_Gen) ...
