function G=build_fishery

%1 = gill net, 2 = spear gun, 3 = hand line, 4 = ring net, 5 = seine net
G=zeros(5);
G(1,4)=0.43;
G(3,4)=0.72;
G(3,5)=1.01;
G(4,1)=0.43;
G(4,3)=0.72;
G(4,5)=0.36;
G(5,3)=1.01;
G(5,4)=0.36;

end


