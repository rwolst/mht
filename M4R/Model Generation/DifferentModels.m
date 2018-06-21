%MODELS

%%%MATSUDA EMP DATA (VAR_{5}(1))
x = 0.1;
A(:,:,1) = [0.2 0 0.3 0 0.3; 0.3 -0.2 x 0 0; 0.2 x 0.3 0 0;0.2 0.3 0 0.3 0;0.2 0 0.2 0.2 0.2];
%(2,5) and (3,4) missing

%%%Chen and Walden (VAR_{3}(2))
A(:,:,2) = [0.267 0.4 0; 0.4 0.179 0; 0 0.4 0.179];
A(:,:,1) = [-0.25 -0.2 0; 0.2 -0.111 0; 0 0.2 -0.111};
%(1,3) missing

%%%Schelter et al (VAR_{5}(4))
A(:,:,4) = [0.6 0 0 0 0 ; 0 0.5 0 0.6 0; 0 0 0.8 0 0; 0 0 0 0.5 0; 0 0 -0.2 0 0.7];
A(:,:,3) = [0 0.65 0 0 0; 0 -0.3 0 0 0; 0 0 -0.7 0 0; 0 0 0.9 0 0.4; 0 0 0 0 -0.5];
A(:,:,2) = [0 0 0 0 0 ; 0 0 0 0 0 ; 0 0 0 0 -0.1; 0 0 0 0 0; 0 0 0 0 0];
A(:,:,1) = [0 0 0 0 0; 0 0 -0.3 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
%(1,3), (1,4), (1,5), (2,5) missing

%%%Stationary VAR_{2}(2)
A(:,:,2) = [1/2 1; 0 1/2];
A(:,:,1) = [0 0; 1/64 0];

%%%Stationary AR(2)
A(:,:,2) = 1/2;
A(:,:,1) = 1/10;

%%%Stationary VAR3(1)
A = [0.5 0.5 0; 0.5 0.25 0; 0 0 0.25];

%%%Stationary VAR4(1)
A = [0.1 0 0.2 0; 0 -0.1 0.3 0.5; -0.1 0 0.1 0; 0 0.2 0 -0.4];