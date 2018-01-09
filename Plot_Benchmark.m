% plot benchmarkfunction
clc,clear
func_num = 12;
if (func_num == 1 || func_num == 2 || func_num == 4 || func_num == 8 || func_num == 9 || func_num == 11 || func_num == 12 || func_num == 15)
    rang_l = -100;
    rang_r = 100;
end
if (func_num == 3)
    rang_l = -600;
    rang_r = 600;
end
if (func_num == 5)
    rang_l = -1;
    rang_r = 4;
end
if (func_num == 10)
    rang_l = -5;
    rang_r = 10;
end
if (func_num == 6 || func_num == 7 || func_num == 13 || func_num == 14)
    rang_l = -10;
    rang_r = 10;
end
d = (rang_r - rang_l)/200;
x = rang_l:d:rang_r;
y = x;
for i = 1:201
    for j = 1:201
        vector = [x(i);y(j)];
        z(i,j) = benchmark_func(vector,func_num);
    end
end

surf(x,y,z)
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
shading interp