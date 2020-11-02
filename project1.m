%could do input checking for all the values inputted, but we can assume
%they are mechanical engineers and just write in the documentation about
%the allowed values. 

%documentation should discuss the inputs and what they mean, include the
%image from teh project guidelines.  Discuss the allowed values for the
%inputs.  Discuss the output and how to interpret it.  

clear;
alpha2 = 5; %degrees
p1x = 9; %inches
p1y = 17.5; %inches
p2x = 4; %inches
p2y = 16; %inches

space = 18; %inches, sidelength of the containment square area
step = 1; %degree, step size of each iteration

%part 1, max assumption of variables
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
                     wz_data(data_pointer,:) = [w,z,a1x,a1y,a2x,a2y,otwox,otwoy]; % wz solutions array output
                     data_pointer = data_pointer +1;
                 end
              end
          end
          
      end
   end
end

%part 2  !!!it is worth noting that if the max angular conditions are the
%same for part 1 and part 2 the solutions generated will be the same, could
%speed the program up if that is assumed to be true. and just combine 
data_pointer = 1;
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
                     us_data(data_pointer,:) = [u,s,b1x,b1y,b2x,b2y,ofourx,ofoury]; %data file colums output
                     data_pointer = data_pointer +1;
                 end
              end
          end
          
      end
   end
end

%combonation of the two solution sets, need all criteria to be met
data_pointer = 1;
for us = 1:length(us_data(:,8)) %tells up number of us solutions
    for wz = 1:length(wz_data(:,8)) %tells up number of wz solutions
        %%%!!!! require that b1 must be +x (to the right) of a1 and that v1
        % and g1 are non zero
        g1 =  sqrt((us_data(us,7) - wz_data(wz,7))^2 + (us_data(us,8) - wz_data(wz,8))^2);
        if us_data(us,3) - wz_data(wz,3) > 0 && g1 > 0
            solutions(data_pointer, :,:) = [us_data,wz_data];
        end
    end
end