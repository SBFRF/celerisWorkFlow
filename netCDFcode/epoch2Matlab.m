function dnum = epoch2Matlab( epoch )

% Usage:
%
%    dnum = epoch2Matlab( epoch )
%
%   convert a UNIX epoch time as used in the Argus system into a Matlab
%   datenum. You can convert the datenum to whatever you want.
%
%   WARNING: DANGER: UNIX epoch times are all GMT based, and the datenum
%   returned from this routine is GMT. 
%

% get reference time
unixEpoch = datenum( '1-jan-1970 00:00:00' );

% how much later than reference time is input?
offset = epoch / (24*3600);

% add and return
dnum = unixEpoch + offset;

%
% Copyright by Oregon State University, 2002
% Developed through collaborative effort of the Argus Users Group
% For official use by the Argus Users Group or other licensed activities.
%
% $Id: epoch2Matlab.m,v 1.11 2004/03/25 16:55:25 stanley Exp $
%
% $Log: epoch2Matlab.m,v $
% Revision 1.11  2004/03/25 16:55:25  stanley
% auto insert keywords
%
%
%key time 
%comment  Converts epoch time to Matlab datenum 
%
