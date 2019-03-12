% Illuminant Estimate Angular Error

T = [0,1,0];
E = [1,0,0];

ThetaInDegrees = atan2d(norm(cross(T,E)),dot(T,E));