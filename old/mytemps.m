function mytemps(uipanel_handle)
c = uicontrol(uipanel_handle,'Style','popupmenu');
%c.Position = position;
c.String = {'Celsius','Kelvin','Fahrenheit'};
c.Callback = @selection;

    function selection(src,event)
        val = c.Value;
        str = c.String;
        str{val};
        disp(['Selection: ' str{val}]);
    end

end