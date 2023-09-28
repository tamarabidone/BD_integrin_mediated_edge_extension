function FormatedNumber = SimFormat(Num)
    
    FormatedNumber = sprintf('%08.3f',Num);
    FormatedNumber(FormatedNumber == '.') = 'd';

end