function fluid_property_at_T_ref = property_interpolating_function(T_ref,imp_fluid_temperature,imp_fluid_property)
% Interpolate the fluid property to the appropriate temperature

[~, loc] = min(abs(imp_fluid_temperature - T_ref));

% If the temperature is out of range return an error 
if T_ref > max(imp_fluid_temperature) || T_ref < min(imp_fluid_temperature)
    error('temperture is out of range')
    
% Do the interpolation    
elseif imp_fluid_temperature(loc) >= T_ref
    
    scaling = (T_ref - imp_fluid_temperature(loc-1))/(imp_fluid_temperature(loc) - imp_fluid_temperature(loc-1));
    fluid_property_at_T_ref = imp_fluid_property(loc-1) + scaling*(imp_fluid_property(loc) - imp_fluid_property(loc-1));
    
else
    scaling = (T_ref - imp_fluid_temperature(loc))/(imp_fluid_temperature(loc+1) - imp_fluid_temperature(loc));
    fluid_property_at_T_ref = imp_fluid_property(loc) + scaling*(imp_fluid_property(loc+1) - imp_fluid_property(loc));
    
end   


end