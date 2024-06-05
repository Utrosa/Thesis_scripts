function [PLh, PLc, PLe]=PLextension(PL)
PLh=round(PL/2)+[0 -rem(PL,2)]; % pad length for temporally extending a QPP
PLc=round(PL/2)+(1:PL); % range of interest in an extended QPP, matches PLh
PLe=PL+sum(PLh); % length of an extended QPP (derivable but saved/read)