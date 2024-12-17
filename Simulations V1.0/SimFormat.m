function FormatedNumber = SimFormat(Num)
    
    FormatedNumber = sprintf('%08.4f',Num);
    FormatedNumber(FormatedNumber == '.') = 'd';

end