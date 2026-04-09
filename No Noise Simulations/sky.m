classdef sky
    methods (Static = true)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % % % % % % % % % % % % Simulation Code % % % % % % % % % % % % %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [X,Y,Rho,Phi] = MakeGrid(H,V,dx,dy)
            % This function creates the Cartesian and Polar coordinate
            % grids for plotting

            % H - number of horizontal pixels
            % V - number of vertical points
            % dx, dy - pixel sizes in the H and V directions
            
            x = linspace(-(H/2),(H/2)-1,H).*dx; % defines x spacing
            y = linspace(-(V/2),(V/2)-1,V).*dy; % define y spacing
            x = x + dx/2 ; % centre number line
            y = y + dy/2 ; % centre number line
            [X,Y] = meshgrid(x,y)   ; % creates cartesian coordinates            
            [Phi,Rho]=cart2pol(X,Y) ; % creates polar coordinates
        end

        function U = HG(X,Y,N,M,weights,w0)

            % This function computes a superposition of HG mode at the plane z=0
            % X,Y are matrices, e.g. x = -1:0.1:1, y = -1:0.1:1, [X,Y] = meshgrid(x,y)
            % "weights" is a weight vector for the coefficients in the superposition
            % E.g. generate mode U = 5*HG_0^-5 + I*HG_0^3 - HG_1^4 + 2*HG_3^4 then,
            % p = [0,0,1,3], l = [-5,3,4,4], weights = [5,I,-1,2]
            
            U = zeros(size(X)); %initialise electric field
             
            for i = 1:length(weights)
                U = U + weights(i).*(1./w0).*(sqrt(2.^(1-N(i)-M(i))./(pi.*factorial(N(i)).*factorial(M(i)))).*sky.Hermite(N(i),sqrt(2).*X./w0).*sl.Hermite(M(i),sqrt(2).*Y./w0).*exp(-(X.^2+Y.^2)./w0.^2));
            end
        
        end

        function y = Hermite(n,x)
            % Computes Hermite polynomials
            if n == 0
                y = ones(size(x));
            elseif n == 1
                y = 2.*x;
            else
                y = 2.*x.*sky.Hermite(n-1,x) - 2.*(n-1).*sky.Hermite(n-2,x);
            end

        end

        function U = LG(R,Phi,P,L,weights,w0)

            % This function computes a superposition of LG mode at the plane z=0
            % R,Phi are coordinate matrices
            % "weights" is a weight vector for the coefficients in the superposition
            % E.g. generate mode U = 5*LG_0^-5 + I*LG_0^3 - LG_1^4 + 2*LG_3^4 then,
            % x = -1:0.01:1, y = -1:0.01:1, [X,Y] = meshgrid(x,y); [Phi,R] = cart2pol(X,Y);
            % P = [0,0,1,3], L = [-5,3,4,4], weights = [5,I,-1,2]; w0 = 1;
            % U = LG(R,Phi,P,L,weights,w0);
            
            U = zeros(size(R)); % initialise field e
            
            for i = 1:length(weights)
                U = U + weights(i).*(sqrt(2*factorial(P(i))/(pi*factorial(P(i)+abs(L(i))))).*(1/w0).*(sqrt(2).*R./w0).^(abs(L(i))) .* exp(-R.^2./w0^2).*sky.Laguerre(P(i),abs(L(i)),2.*R.^2./(w0.^2)).* exp(1i.*L(i).*Phi));
            end
            
        end  

        function y = Laguerre(p,l,x)
            % computes associated Laguerre polynomials 
            if p == 0
                y = ones(size(x));
            elseif p == 1
                y = 1 + l - x;
            else
                y = ((2*p+l-1-x)./p).*sky.Laguerre(p-1,l,x) - ((p+l-1)./p).*sky.Laguerre(p-2,l,x);
            end
        end

        function [u2] = AngSpecProp(u1,Lx,Ly,lambda,z)

            % u1 is the original field
            % L is length of original plane
            % lambda is the wavelength
            % z is the propagation distance
            
            % Create coordinate grids            
            Mx = size(u1,2); % get size of original field
            My= size(u1,1); % get size of original field            
            dx = Lx/Mx; % spacing of the starting spatial grid
            dy = Ly/My; % spacing of the starting spatial grid
            fx = -1/(2*dx):1/Lx:1/(2*dx)-1/Lx;
            fy = -1/(2*dy):1/Ly:1/(2*dy)-1/Ly;% frequency coordinates
            [Fx,Fy] = meshgrid(fx,fy);
            
            % Create propogation phase in freqeucnyy domain
            k = 2*pi/lambda;
            H = exp(1i*k*z*sqrt(1-(lambda^2)*(Fx.^2 + Fy.^2)));
            H = fftshift(H); % shift the propagation phase for Matlab's fft function
            
            % Propagate field            
            U1 = fft2(fftshift(u1)); % Fourier transform original function into frequency domain
            U2 = U1.*H; % add propagation phase
            u2 = ifftshift(ifft2(U2)); % Shift function back from frequency domain to spatial domain
        end

        function [s0,s1,s2,s3] = GetStokesRL(UR,UL)
    
            UH = (1/sqrt(2))*(UR + UL);
            UV = (1/sqrt(2))*(UR - UL);
            UD = (1/sqrt(2))*(UR - 1i.*UL);
            UA = (1/sqrt(2))*(UR + 1i.*UL);
            
            H = abs(UH).^2;
            V = abs(UV).^2;
            Tot = max(max(H+V));
            H = H/Tot;
            V = V/Tot;
            
            R = abs(UR).^2;
            R = R/Tot;
            L = abs(UL).^2;
            L = L/Tot;
            D = abs(UD).^2;
            D = D/Tot;
            A = abs(UA).^2;
            A = A/Tot;
            
            s0 = H + V;
            s1 = H - V;
            s2 = D - A;
            s3 = R - L;
        
        end

        function [s0,s1,s2,s3] = GetStokesHV(UH,UV)
            
            UR = (1/sqrt(2))*(UH + 1i*UV);
            UL = (1/sqrt(2))*(UH - 1i*UV);
            UD = (1/sqrt(2))*(UH + UV);
            UA = (1/sqrt(2))*(UH - UV);
            
            H = abs(UH).^2;
            V = abs(UV).^2;
            Tot = max(max(H+V));
            H = H/Tot;
            V = V/Tot;
            
            R = abs(UR).^2;
            R = R/Tot;
            L = abs(UL).^2;
            L = L/Tot;
            D = abs(UD).^2;
            D = D/Tot;
            A = abs(UA).^2;
            A = A/Tot;
            
            s0 = H + V;
            s1 = H - V;
            s2 = D - A;
            s3 = R - L;
        
        end

        function [Vcharge, singularities] = VortexCharge(Field, pixcrop, BeamSizeInPix)
            N = size(Field);
          
            xcrop = (N(1)/2-pixcrop) : (N(1)/2+pixcrop);
            ycrop = (N(2)/2-pixcrop) : (N(2)/2+pixcrop);
        
            xcrop = round(xcrop);
            ycrop = round(ycrop);
            
            Fieldcopy = Field(xcrop, ycrop);
            Fieldcopy(sqrt(xcrop.^2 + ycrop.^2) > BeamSizeInPix ) = 0;
        
            dGx = angle(ShiftHolo(Fieldcopy, 1, 0).*ShiftHolo(conj(Fieldcopy), 0,0));% x gradient;
            dGy = angle(ShiftHolo(Fieldcopy, 0, 1).*ShiftHolo(conj(Fieldcopy),0, 0));%x gradient;
        
            Circulation = round((ShiftHolo(dGx - dGy, 0, 0) - ShiftHolo(dGx, 0, 1) +ShiftHolo(dGy, 1,0))./ (2*pi) );% compute circulation
            Circulation( 1:3, :) = 0;
            Circulation( :, 1:3) = 0;
        
            singularities = Circulation;
            Vcharge = sum(sum(Circulation));

            function Temp = ShiftHolo(Holo, HDispl, VDispl)
                n = size(Holo);
                Temp = zeros(n);
                
                X = 1:n(1);
                Y = 1:n(2);
                 
                XPrim = X-HDispl;
                YPrim = Y-VDispl;
                
                Xclean = XPrim(XPrim<=n(1) & XPrim>=1);
                Yclean = YPrim(YPrim<=n(2) & YPrim>=1);
                Temp(Xclean, Yclean)=Holo(Xclean+HDispl, Yclean+VDispl);
            end

        end

        function [SkyNum, s1, s2, s3,singularities] = SkyNumLine(S1,S2,S3,S0,pixcrop,BeamSizeInPix, thresh)
            N = size(S1);
            xcrop = (N(1)/2-pixcrop) : (N(1)/2+pixcrop);
            ycrop = (N(2)/2-pixcrop) : (N(2)/2+pixcrop);
            xcrop = round(xcrop);
            ycrop = round(ycrop);
            s1 = S1(xcrop, ycrop);
            s2 = S2(xcrop, ycrop);
            s3 = S3(xcrop, ycrop);
            s0 = S0(xcrop, ycrop);
            Norm = sqrt(s1.^2 + s2.^2 + s3.^2);
            s1 = s1./Norm;
            s2 = s2./Norm;
            s3 = s3./Norm;
            Field = s1+1i.*s2;
            Field(isinf(Field)|isnan(Field)) = 0;
            s3(isinf(s3)|isnan(s3)) = 0;
            [~, singularities] = sky.VortexCharge(exp(1i.*angle(Field)), pixcrop, BeamSizeInPix);
            singularities(isinf(singularities)|isnan(singularities)) = 0;
            singularities(s0 < thresh*max(s0(:))) = 0;
            Inner = sum(sum(singularities.*s3));
            OuterPath = ones(size(s1)-2);
            OuterPath = abs(padarray(OuterPath,[1 1])-1);
            s3Inf = sum(sum(OuterPath.*s3))./sum(sum(OuterPath == 1));
            Outer=sum(sum(singularities)).*s3Inf;
            SkyNum = (1/2).*(Inner - Outer);
        end
                
        function [n ,sigz] = SkyNumSurf(S1,S2,S3)

            %Local Normalization of Stokes parameters
            a = sqrt(S1.^2 + S2.^2 + S3.^2);
            S1 = S1./a;
            S2 = S2./a;
            S3 = S3./a;
            
            %Gradients of Stokes parameters
            [dS1_dx, dS1_dy] = gradient(S1);
            [dS2_dx, dS2_dy] = gradient(S2);
            [dS3_dx, dS3_dy] = gradient(S3);
                
            %Calculating z-component of Skyrmion field
            t1 = S1.*dS2_dx.*dS3_dy;
            t2 = S2.*dS1_dy.*dS3_dx;
            t3 = S3.*dS1_dx.*dS2_dy;
            t4 = S1.*dS3_dx.*dS2_dy;
            t5 = S2.*dS1_dx.*dS3_dy;
            t6 = S3.*dS2_dx.*dS1_dy;
            sigz = t1 + t2 + t3 - t4 - t5 - t6; %Skyrmion Field, sigz
            
            %Removing numerical artifacts
            sigz(isnan(sigz))=0;
            sigz(isinf(sigz))=0;
            
            %Skyrme number calculation
            n = -(1/(4*pi))*sum(sum(sigz));
                
        end
        
    end
end
