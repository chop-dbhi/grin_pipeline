function divfunc(a, b)
	if(a == 0) then 
		return 0.0
	else 
    	return string.format("%.9f", a / b)
    end
end