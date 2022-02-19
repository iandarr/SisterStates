function [chosen_row,chosen_col]=choose_reaction(rates_all,rand_number)
nrow=size(rates_all,1);
ncol=size(rates_all,2);

rates_vector=reshape(rates_all,1,numel(rates_all));
rates_vector_cumsum_norm=cumsum(rates_vector)/sum(rates_vector);

vector_index=find(rand_number<rates_vector_cumsum_norm,1,'first');
chosen_col=ceil(vector_index/nrow);
chosen_row=vector_index - (chosen_col-1)*nrow;

if rates_all(chosen_row,chosen_col)==0,error("theres a problem with the choose_reaction function"),end

end