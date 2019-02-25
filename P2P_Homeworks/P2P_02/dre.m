function f = dre(s)

s = vpa(s);

sym_vars = symvar(s);

for i=1:length(sym_vars)
   eval(['syms ' char(sym_vars(i))]) 
end

for i=1:size(s,1)
    for j = 1:size(s,2)
        k = char(sym(s(i,j)));
        % k([1:6]) = [];
        
        l = strfind(k,'e-');
        
        count =0;
        
        while ~isempty(l)
            count=count+1;
            
            if length(k) >= l(1)+3
                val = k(l(1)+2:l(1)+3);
            else
                val = k(l(1)+2:l(1)+2);
            end
            
            if numel(val) == 1
                
                if str2double(val)>4
                    k = strrep(k,k(l(1):l(1)+1),'*0*');
                end
                
            else
                
                if (val(2) == '*' && str2double(val(1))>4) || str2double(val) > 4
                    k = strrep(k,k(l(1):l(1)+1),'*0*');
                end
                
            end
            
            l = strfind(k,'e-');
            
        end
        
        % f = simplify(sym(k));
%         f = simplify(sym(k));
        f(i,j) = simplify(vpa(eval(k)));
        
    end
end
end