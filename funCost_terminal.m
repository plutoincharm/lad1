function Ct = funCost_terminal(x)


global XT  r L MI g nx ny lam Q R Qf 


Ct = 0.5*(((XT-x)')*Q*(XT-x));