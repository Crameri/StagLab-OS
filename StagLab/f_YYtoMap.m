
%%                                                 YINYANG GRID TO MAP 2.1
%
%                                                             Paul Tackley
%                                 modified by Fabio Crameri, December 2015
%
%% NOTES
% this version does not require the theta and phi arguments, and expects
% ayy to be 3-dimensional ayy(itheta,iphi,iblock)

function [amap,theta,phi] = f_YYtoMap(ayy)
nth=size(ayy,1);
nph=nth*3;

%% load yin grid directly into map

nthmap=2*nth; nphmap=(4*nph)/3;

amap = zeros(nthmap,nphmap);

theta = zeros(nthmap,1);
phi = zeros(nphmap,1);
for i=1:nthmap
    theta1=pi*(i-0.5)/nthmap;   % range 0 to pi
    %****fc***************************** make sure to range from 0 to pi
    if i==1
        theta1=0;
    elseif i==nthmap
        theta1=pi;
    end
    %****fc*****************************
    for j=1:nphmap
        phi1=2*pi*(j-0.5)/nphmap-pi; % range -pi to +pi
        %****fc***************************** make sure to range from -pi to +pi
        if j==1
            phi1=-pi;
        elseif j==nphmap
            phi1=pi;
        end
        %****fc*****************************
        theta2=acos(sin(theta1)*sin(phi1));
        phi2=atan2(cos(theta1),-sin(theta1)*cos(phi1));
        
        ongrid1 = (theta1>0.25*pi) && (theta1<0.75*pi) && (abs(phi1)<0.75*pi); %on Yin
        ongrid2 = (theta2>0.25*pi) && (theta2<0.75*pi) && (abs(phi2)<0.75*pi); %on Yang
        
        if ongrid1 && ongrid2
            if abs(phi1)<abs(phi2)
                ongrid2 = 0;
            else
                ongrid1 = 0;
            end
        end
        
        if ongrid1
            % in yin grid
            itheta=round(0.5+nth*(theta1-0.25*pi)/(0.5*pi)); itheta=min(itheta,nth);
            iphi=round(0.5+nph*(phi1+0.75*pi)/(1.5*pi))+1; iphi=min(iphi,nph);
            amap(i,j) = ayy(itheta,iphi,1);
        elseif ongrid2
            % in yang grid
            itheta=round(0.5+nth*(theta2-0.25*pi)/(0.5*pi)); itheta=min(itheta,nth);
            iphi=round(0.5+nph*(phi2+0.75*pi)/(1.5*pi))+1; iphi=min(iphi,nph);
            amap(i,j) = ayy(itheta,iphi,2);
        end
        
        theta(i) = theta1;
        phi(j)   = phi1;
    end
end
