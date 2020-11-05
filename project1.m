%could do input checking for all the values inputted, but we can assume
%they are mechanical engineers and just write in the documentation about
%the allowed values. 

%documentation should discuss the inputs and what they mean, include the
%image from teh project guidelines.  Discuss the allowed values for the
%inputs.  Discuss the output and how to interpret it.  
tic
clear;
alpha2 = 5; %degrees
p1x = 10; %inches
p1y = 15; %inches
p2x = 4; %inches
p2y = 15; %inches
space = 18; %inches, sidelength of the containment square area

%max number of solutions to be checked, default settings yields 1.5e9 possible solutions
%ERROR: if num_solutions is larger than your PC memory capacity the program
% will through an error
num_solutions = 20; %number of desired solutions possible solutions to be saved., if -1 it will equal solutions to check
solutions_to_check = 1e7; % if value is -1, it checks all possible solutions

step = 1; %degree, step size of each iteration
%part 1, max assumption of variables 180 default for all values
theta_max = 180; %degrees
phi_max = 180; %degrees
beta2_max = 180; %degrees

%part 2, max assumption of variables
sigma_max = 180; %degrees
gamma2_max = 180; %degrees
psi_max = 180; %degrees

p21 = sqrt((p2x-p1x)^2 + (p2y-p1y)^2); %magnitude of vector P21

%calculation of delta2 depends on which quadrant p2 is with respect to p1
%p2 always has to be to the left of p1 (quadrant 2,3)
if (p2x-p1x) <= 0
    if (p2y-p1y) <=0
        delta2 = atand((p2y-p1y)/(p2x-p1x)) + 180;
    else
        delta2 = -1* atand((p2x-p1x)/(p2y-p1y)) + 90;
    end
else
   %error code thrown to the user and the program is terminated
   error("Error 1: The location of P1 and P2 that you chose are not valid. P1 needs to be right of P2")
end

%part 1 
data_pointer = 1;
wz_data = zeros(theta_max*beta2_max * phi_max, 11); %preallocate max size to improve runtime

for theta = 0:step:theta_max %iterate of all values of theta, phi and beta2
   for phi = 0:step:phi_max
      for beta2 = 0:step:beta2_max
          a = cosd(theta) * (cosd(beta2) - 1)+ sind(theta) * sind(beta2); %equations from notes on generation
          b = cosd(phi) * (cosd(alpha2) - 1) + sind(phi) * sind(alpha2);
          c = p21*cosd(delta2);
          d = sind(theta)*(cosd(beta2) -1) + cosd(theta) * sind(beta2);
          e = sind(phi)*(cosd(alpha2)-1) + cosd(phi) * sind(alpha2);
          f = p21*sind(delta2);
          w = (c*e - b*f)/(a*e-b*d);
          z = (a*f-c*d)/(a*e-b*d);
          
          %initial filtering of wz values
          a1x = p1x - z*cosd(phi);
          a1y = p1y - z*sind(phi);
          otwox = a1x - w*cosd(theta);
          otwoy = a1y - w*sind(theta);
          a2x = otwox + w*cosd(theta+beta2);
          a2y = otwoy + w*sind(theta+beta2);
          
          %force w z to be positive
          if a2x > 0 && a2x < space && a2y > 0 && a2y < space && w > 0 && z > 0
              if a1x > 0 && a1x < space && a1y > 0 && a1y < space
                 if otwoy > 0 && otwoy < space && otwox > 0 && otwox < space
                     temp = [w,z,a1x,a1y,a2x,a2y,otwox,otwoy,theta, phi, beta2]; % wz solutions array output
                     wz_data(data_pointer,:) = temp; 
                     data_pointer = data_pointer +1;
                 end
              end
          end
          
      end
   end
end

%wz_data solution size reduction
wz_data = wz_data(1:(data_pointer-1), :);

%part 2  !!!it is worth noting that if the max angular conditions are the
%same for part 1 and part 2 the solutions generated will be the same, could
%speed the program up if that is assumed to be true. and just combine 
data_pointer = 1;
us_data = zeros(sigma_max*gamma2_max*psi_max, 11);

for sigma = 0:step:sigma_max
   for gamma2 = 0:step:gamma2_max
      for psi = 0:step:psi_max
          a = cosd(sigma) * (cosd(gamma2) - 1)+ sind(sigma) * sind(gamma2);
          b = cosd(psi) * (cosd(alpha2) - 1) + sind(psi) * sind(alpha2);
          c = p21*cosd(delta2);
          d = sind(sigma)*(cosd(gamma2) -1) + cosd(sigma) * sind(gamma2);
          e = sind(psi)*(cosd(alpha2)-1) + cosd(psi) * sind(alpha2);
          f = p21*sind(delta2);
          u = (c*e - b*f)/(a*e-b*d);
          s = (a*f-c*d)/(a*e-b*d);
          
          %initial filtering of wz values
          b1x = p1x - s*cosd(psi);
          b1y = p1y - s*sind(psi);
          ofourx = b1x - u*cosd(sigma);
          ofoury = b1y - u*sind(sigma);
          b2x = ofourx + u*cosd(sigma+gamma2);
          b2y = ofoury + u*sind(sigma+gamma2);
          
          %force u s to be positive
          if b2x > 0 && b2x < space && b2y > 0 && b2y < space && u > 0 && s > 0
              if b1x > 0 && b1x < space && b1y > 0 && b1y < space
                 if ofourx > 0 && ofourx < space && ofoury > 0 && ofoury < space
                     temp = [u,s,b1x,b1y,b2x,b2y,ofourx,ofoury,sigma,gamma2,psi]; %data file colums output
                     us_data(data_pointer,:) = temp;
                     data_pointer = data_pointer +1;
                 end
              end
          end
          
      end
   end
end

%us_data size reduction
us_data = us_data(1:(data_pointer-1), :);

possible_solutions = length(us_data(:,11)) * length(wz_data(:,11));
if solutions_to_check == -1 || solutions_to_check > possible_solutions
   solutions_to_check = possible_solutions; 
end
if num_solutions == -1 || num_solutions > solutions_to_check
    num_solutions = solutions_to_check;
end

solutions = zeros(num_solutions,25 ); %preallocate memory for solution matrix and it vastly improves runtime, i cut unused size down later
min_solution = 15000* ones(1,25); %preallocate min solution
 
%in order to get varied solutions, a probablility factor is created and
%check against a random probablility
prob = num_solutions/solutions_to_check *1; %the integer is just a weight factor: ensures that the array is filled 

%in order to get varied solutions, set step size to step across possible
%wz/us solutions, so the whole solution set of 1e9 can be checked
dus = sqrt(possible_solutions/solutions_to_check);
dus = floor(dus);
dwz = dus;

%combonation of the two solution sets, need all criteria to be met
data_pointer = 1;
solutions_checked = 0;
for us = 1:dus:length(us_data(:,8)) %tells up number of us solutions
    for wz = 1:dwz:length(wz_data(:,8)) %tells up number of wz solutions
        solutions_checked = solutions_checked +1; 
        %%%!!!! require that b1 must be +x (to the right) of a1 and that v1
        % and g1 are non zero
        % solution index comment on end of line
        w = wz_data(wz,1); %1
        z = wz_data(wz,2);%2
        u = us_data(us,1);%3
        s = us_data(us,2);%4
        otwox = wz_data(wz,7);%5
        otwoy = wz_data(wz,8);%6
        ofourx = us_data(us,7);%7
        ofoury = us_data(us,8);%8
        a1x = wz_data(wz,3);%9
        a1y = wz_data(wz,4);%10
        b1x = us_data(us,3);%11
        b1y = us_data(us,4);%12
        a2x = wz_data(wz,5);%13
        a2y = wz_data(wz,6);%14
        b2x = us_data(us,5);%15
        b2y = us_data(us,6);%16
        g1 =  sqrt((ofourx - otwox)^2 + (ofourx - otwox)^2);%17
        v1 = sqrt((a1x - b1x)^2 + (a1y - b1y)^2);%18
        length_total = w+v1+u;%19
        theta = wz_data(wz,9);%20
        phi = wz_data(wz,10);%21
        beta2 = wz_data(wz,11);%22
        sigma = us_data(us,9);%23
        psi = us_data(us,10);%24
        gamma2 = us_data(us,11);%25
        
        temp = [w,z,u,s,otwox,otwoy,ofourx,ofoury,a1x,a1y,b1x,b1y,a2x,a2y,b2x,b2y,g1,v1,length_total,theta,phi,beta2,sigma,psi,gamma2];
        
        if b1x - a1x > 0 && g1 > 0 && v1 > 0
            if data_pointer <= num_solutions && rand <= prob %only save so many solutions
                solutions(data_pointer, :) = temp;
                data_pointer = data_pointer +1;
            end
            %will check all possible solutions for minimum. 
            if length_total < min_solution(1,19) %minimum length solution checking and tracking
                min_solution = temp;
            end
        end
        
    end
end

%solution matrix size reduction
solutions = solutions(1:(data_pointer-1), :);
toc