function [parab] = genParab(score_neig, dispar_int_val, params)
% Calculate paraoblas 

parab = subPixInterp(score_neig, params); % Sub-pixel interpolation of current scores

% Shift parabaols by disparaity values so all parabolas are centered around zero.
aNew = parab.a;
bNew = parab.b-2*parab.a.*dispar_int_val;
cNew = parab.a.*dispar_int_val.^2-parab.b.*dispar_int_val+parab.c;
parab.a = aNew;
parab.b = bNew;
parab.c = cNew;

end