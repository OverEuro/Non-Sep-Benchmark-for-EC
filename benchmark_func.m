%% Shifted non-fully separable benchmark
%% by Yijun Yang
% Dimension=2 for ploting, if you want to benchmark your own algorithm,
% please change the length of "S"(shifted vector),e.g. the line 22 "
% S(1:2,5) ---> S(1:1000,5), if the dimension = 1000".
%--------------------------------------------------------------------------


function fitness = benchmark_func(x, func_num)

	fitness = feval(['f' num2str(func_num)], x);

end

%-----------------------------Base functions-------------------------------

%--------------------------------------------------------------------------
% no.5 Brown Function OK [-1,4]
%--------------------------------------------------------------------------
function fitness = brown(x)
load S.mat
optimum = S(1:2,5);
x = x - optimum;

x1 = x(1:end-1).^2;
x2 = x(2:end).^2;
x3 = x1.^(x2 + 1);
x4 = x2.^(x1 + 1);
fitness = sum(x3) + sum(x4);
end

%--------------------------------------------------------------------------
% No.6 Mishra11 Function OK [-10,10]
%--------------------------------------------------------------------------
function fitness = Mishra(x)
load S.mat
optimum = S(1:2,6);
x = x - optimum;
[D , ~] = size(x);

fitness = (sum(abs(x))/D + (prod(abs(x).^(1/D)))^(1/D))^2;
end


%--------------------------------------------------------------------------
% No.8 Schwefel04 Function OK
%--------------------------------------------------------------------------
function fitness = schwefel04(x)
load S.mat
optimum = S(1:2,8);
x = x - optimum;
x = x + 1;

fitness = sum((x-1).^2 + (x.^2-x(1)).^2);
end


%--------------------------------------------------------------------------
% No.1 Schwefel's Problem 1.2 OK
%--------------------------------------------------------------------------
function fitness = schwefel(x)
load S.mat
optimum = S(1:2,1);
x = x - optimum;
[D , ~] = size(x);

x = repmat(x',D,1);
x = tril(x);
x = sum(x,2);
fitness = sum(x.^2);
end

%--------------------------------------------------------------------------
% No.2 Rosenbrock's Function 不需要改
%--------------------------------------------------------------------------
function fitness = rosenbrock(x)
load S.mat
optimum = S(1:2,2);
x = x - optimum;
x=x+1;
[D , ~] = size(x);

fitness = sum(100.*(x(1:D-1,:).^2-x(2:D, :)).^2+(x(1:D-1, :)-1).^2);
end

%--------------------------------------------------------------------------
%No.3 Griewank Function 不需要改
%--------------------------------------------------------------------------
function fitness = griewank(x)
load S.mat
optimum = S(1:2,3);
x = x - optimum;
[D , ~] = size(x);

fitness=1+1./4000.*sum(x.^2)-prod(cos(x./sqrt(1:D)'));

end


%--------------------------------------------------------------------------
%No.4 Sharp Ridge Function (f(1,...,1)=0) OK
%--------------------------------------------------------------------------
function fitness = sharp(x)
load S.mat
optimum = S(1:2,4);
x = x - optimum;

fitness = x(1)^2 + 100 * sqrt(sum(x(2:end).^2));
end

%--------------------------------------------------------------------------
% No.7 Multimodal Function f(0,...,0)=0; OK [-10,10]
%--------------------------------------------------------------------------
function fitness = Multimodal(x)
[D , ~] = size(x);
load S.mat
optimum = S(1:2,7);
x = x - optimum;

fitness = sum(abs(x)) * prod(abs(x).^(1/D));
end

%--------------------------------------------------------------------------
% No.9 Salomon Function f(0,...,0)=0; No.110 in the papaer
%--------------------------------------------------------------------------
function fitness = salomon(x)
load S.mat
optimum = S(1:2,9);
x = x - optimum;
fitness = sum(x.^2);
fitness = 1 - cos(2*pi*sqrt(fitness)) + 0.1 * sqrt(fitness);
end

%--------------------------------------------------------------------------
% No.10 Zacharov Function f(0,...,0)=0; OK [-5,10]
%--------------------------------------------------------------------------
function fitness = zacharov(x)
load S.mat
optimum = S(1:2,10);
x = x - optimum;
[D , ~] = size(x);
% fitness = sum(x.^4);
w = 1:1:D;
w = w';
fitness = sum(x.^2) + (0.5 * sum(x .* w))^2 + (0.5 * sum(x .* w))^4;
end

%--------------------------------------------------------------------------
% No.11 Schwefel 2.22 Function f(0,...,0)=0; No.124 in the papaer
%--------------------------------------------------------------------------
function fitness = schwefel2(x)
[D , ~] = size(x);
load S.mat
optimum = S(1:2,11);
x = x - optimum;
fitness = abs(x);
fitness = sum(fitness) + prod((fitness).^(1/D));
end

%--------------------------------------------------------------------------
% No.12 Schaffer f6 Function f(0,...,0)=0; No.136 in the papaer
%--------------------------------------------------------------------------
function fitness = schafferf6(x)
load S.mat
optimum = S(1:2,12);
x = x - optimum;
[D,~]=size(x);
fitness = 0;
for i = 1:D-1
    fitness = fitness + (0.5 + ((sin(sqrt(x(i)^2 + x(i+1)^2)))^2-0.5)/((1 + 0.001 * (x(i)^2 + x(i+1)^2))^2));
end
end

%--------------------------------------------------------------------------
% No.13 zerosum Function f(0,...,0)=0; OK [-10,10]
%--------------------------------------------------------------------------
function fitness = zerosum(x)
load S.mat
optimum = S(1:2,13);
x = x - optimum;
% fitness = sum(round(x.^2));
if sum(x) == 0
    fitness = 0;
else
    fitness = 1 + (10000*abs(sum(x)))^0.5;
end
end

%--------------------------------------------------------------------------
% No.14 Whitley function f(1,...,1)=0; OK [-10.24,10,24]
%--------------------------------------------------------------------------
function fitness = whitley(x)
load S.mat
optimum = S(1:2,14);
x = x - optimum;
x=x + 1;
[D,~]=size(x);

x1 = repmat(x,1,D);
x2 = repmat(x',D,1);
m = (100*(x1.^2 - x2).^2 + (1-x2).^2).^2/4000 - cos(100*(x1.^2 - x2).^2 + (1-x2).^2) + 1;
fitness = sum(sum(m,2));
end

%--------------------------------------------------------------------------
% No.15 Happy Cat Function f(-1,...,-1)=0;
%--------------------------------------------------------------------------
function fitness = happycat(x)
load S.mat
optimum = S(1:2,15);
x = x - optimum;
x=x - 1;
[D,~]=size(x);
fitness = sum(x.^2);
fitness = abs(fitness - D)^0.25 + (0.5*fitness + sum(x))/D + 0.5;
end
%% Benchmark functions
%
%----------f1 schwefel------------------------
function fitness = f1(x)
fitness=schwefel(x);
end
%

%------------f2 rosenbrock-----------------
function fitness = f2(x)
fitness=rosenbrock(x);
end
%

%------------f3 griewank---------------------
function fitness = f3(x)

fitness = griewank(x);
end
%

%------------f4 sharp ridge---------------------
function fitness = f4(x)
fitness = sharp(x);
end

%------------f5 brown---------------------
function fitness = f5(x)

fitness = brown(x);
end

%-----------f6 Mishra 11 ----------------------
function fitness = f6(x)

fitness = Mishra(x);
end

%-----------f7 Multimodal---------------------------
function fitness = f7(x)

fitness = Multimodal(x);
end

%-----------f8 schwefel 04-----------------
function fitness = f8(x)

fitness = schwefel04(x);
end

%-----------f9 salomon function-----------------
function fitness = f9(x)
fitness = salomon(x);
end

%-----------f10 zacharov function--------
function fitness = f10(x)

fitness = zacharov(x);
end

%-----------f11 Schwefel 2.22 function--------
function fitness = f11(x)
fitness = schwefel2(x);
end

%-----------f12 Schaffer f6 function----------
function fitness = f12(x)
fitness = schafferf6(x);
end

%-----------f13 zerosum function---------------
function fitness = f13(x)

fitness = zerosum(x);
end

%-----------f14 whitley function------
function fitness = f14(x)
fitness = whitley(x);
end

%-----------f15 Happy Cat function------------
function fitness = f15(x)
fitness = happycat(x);
end
