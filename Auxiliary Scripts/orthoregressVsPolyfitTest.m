clc, clear, close all
%Seems like orthogonalRegress is heavily skewed by outliers

x = [0.149689613200379;0.449068839601138;0.748448066001896;1.04782729240265;1.34720651880341;1.64658574520417;1.94596497160493;2.24534419800569;2.54472342440645;2.84410265080720;3.14348187720796;3.44286110360872;3.74224033000948;4.04161955641024;4.34099878281100;4.64037800921176;4.93975723561251;5.23913646201327;5.53851568841403;5.83789491481479;6.13727414121555;6.43665336761631;6.73603259401706;7.03541182041782;7.33479104681858;7.63417027321934;7.93354949962010;8.23292872602086;8.53230795242162;8.83168717882237;9.13106640522313;9.43044563162389;9.72982485802465;10.0292040844254;10.3285833108262;10.6279625372269];
y = [96.8504687144251;-20.9745639700476;-53.6556753799132;-84.6749047093412;-115.995287992002;-198.842910466521;-227.836688940846;-278.383910884915;-299.494550682032;-413.251592862583;-372.099654024158;-419.134732451003;-305.304186088239;-409.920138404274;-421.034683497639;-989.061467530842;-329.827947503054;-387.800059867476;-519.523747769666;-824.721027744531;4450.07725631618;-728.386157492186;7821.40288735142;-480.201038383340;-633.331105536557;-835.972296742431;-637.030361450349;-866.934932324612;-722.811184288278;-745.668291130545;-1404.07147760073;-621.440818238919;-529.321250839922;-1731.66697531493;-1165.14918345045;-660.769607893353];

figure, hold on

scatter(x,y)

orthoFit = orthogonalRegress(x,y);

xlin = linspace(min(x),max(x));
ylin1 = xlin*orthoFit(1) + orthoFit(2);

plot(xlin,ylin1)

pfit = polyfit(x,y,1);
ylin2 = xlin*pfit(1) + pfit(2);
plot(xlin,ylin2)

legend({'Data','Ortho Regress','polyfit'})

ylim([-5000 10000])
%%
% I would like to use pca to do the linear regression but I don't know how

p = pca([x y])
%ylin3 = xlin*

plot(xlin,ylin3)