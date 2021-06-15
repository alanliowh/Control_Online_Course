% MakeWSP;

function [wsp,wsp_file_name] = MakeWSP(wind_no,t)
switch wind_no
    case 0
        wsp_file_name = 'WindFiles/Exercise_wind.hh';
    case 1
        wsp_file_name = 'WindFiles/step.hh';
    case 2
        wsp_file_name = 'WindFiles/Kaimal_8ms.hh';
    case 3
        wsp_file_name = 'WindFiles/Kaimal_12ms.hh';
    case 4
        wsp_file_name = 'WindFiles/Kaimal_15ms.hh';
    case 5
        wsp_file_name = 'WindFiles/Kaimal_18ms.hh';
    case 6
        wsp_file_name = 'WindFiles/Part_4_wind.hh';
    case 7
        wsp_file_name = 'WindFiles/Part_5_wind.hh';
    otherwise
        error('Please give me an appropriate wind profile!!!');
end

wsp_data = load(wsp_file_name);

wsp  = interp1(wsp_data(:,1),wsp_data(:,2),t);
