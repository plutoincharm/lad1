function Cs = funStagecost(x,u)


global XT  r L MI g nx ny lam Q R Qf 


Cs= 0.5*(((XT-x)')*Q*(XT-x)) + 0.5*(u')*R*u;