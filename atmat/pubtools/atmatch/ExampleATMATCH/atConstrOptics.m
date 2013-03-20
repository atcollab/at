function constraint=atConstrOptics(r,dpp,opticstoconstrain)
% output a constraint structure using atfuntofitoptics function
% 
% ring is an AT lattice
% dpp is the DP/P 
% 
% opticstoconstrain={'betx',refpts,valminbetx,valmaxbetx,w,...
%                    'bety',refpts,valminbety,valmaxbety,w,...
%                    'alfx',refpts,valminalfx,valmaxalfx,w,...
%                    'alfy', ........
%                    'dispx', ........
%                    'dispy', ........
%                    'mux', ........
%                    'muy', ........
%                   }
% % parameters to fit may be (atlinopt output variables):
%            'betx',
%            'bety',
%            'alfx',
%            'alfy',
%            'dispx',
%            'dispy',
%            'disppx',
%            'disppy',
%            'gammax',
%            'gammay',
%            'mux', % phase advance at element entrance /2/pi
%            'muy', % phase advance at element entrance /2/pi
%            'xco',  % closed orbit
%            'xpco', 
%            'xco',
%            'xpco',
%            'Qx',  % fractional part tune
%            'Qy',  % fractional part tune
%            'chromX',
%            'chromY',
%            'Mij', (M44 i,j component)
%            'SPos'
%            
%
% the output is a structure with fields of a constriant for atmatch
% 
% usage example: 
% 
% vincoli={'betx',qfmindx(2),17.3,17.3,1,...
%          'bety',qfmindx(2),0.58,0.58,1,...
%          'dispx',1,0,0,1,...
%          'Qx',length(RING),0.35,0.35,1};
%  constraint{i}=atConstrOptics(lattice,0,vincoli);
% 
% 


[~,min,max,w]=atfuntofitoptics(r,dpp,opticstoconstrain(1:5:end),... name
                          opticstoconstrain(2:5:end),... refpts
                          opticstoconstrain(3:5:end),... min
                          opticstoconstrain(4:5:end),... max
                          opticstoconstrain(5:5:end)... w
                          );
                     
constraint.Fun   =@(r)atfuntofitoptics(r,dpp,opticstoconstrain(1:5:end),... name
                          opticstoconstrain(2:5:end),... refpts
                          opticstoconstrain(3:5:end),... min
                          opticstoconstrain(4:5:end),... max
                          opticstoconstrain(5:5:end)... w
                          );
constraint.Min   =min;
constraint.Max   =max;
constraint.Weight=w;




return