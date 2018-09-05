function R=get_reward(Snew,BR)

RSnew=BR(Snew==1,:);
R=length(find(sum(RSnew,1))); % get benefit for the state

end