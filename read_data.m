% function data = read_data(expidx, subjidx, condidx)
%
% INPUT
%  expidx  : experiment number (1, 2, or 3; see paper)
%  subjidx : subject index (in range 1-30, 1-85, and 1-20 in experiments 1, 2, 3, respectively)
%  condidx : condition index, as follows:
%             exp1 -- 1=larger/smaller task, 2=same/different task
%             exp2 -- 1=50ms ISI, 2=300ms ISI, 3=2000ms ISI
%             exp3 -- 1=estimation task
% 
% This file is part of the code published with the paper "Recent is more: 
% a negative time-order effect in non-symbolic numerical judgment" by 
% R. van den Berg, M. Lindskog, L. Poom, and A. Winman (JEP:HPP, 2017).
%
% For questions, bug reports, etc, please email nronaldvdberg@gmail.com

function data = read_data(expidx,subjidx,taskidx)

load VandenBergEtAl2017_data.mat;
data = all_data{expidx}{subjidx}{taskidx};